mod analytic;
mod block;
mod edge;
mod march;

use numpy::ndarray::Array2;
use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const BAND: usize = 5;
const STRIDE: usize = 3 * BAND + 1;

#[derive(Clone)]
struct Viscosity {
    kind: u8,
    parameters: Vec<f64>,
}

impl Viscosity {
    fn evaluate(&self, t: f64, out_of_range: &mut bool) -> Option<f64> {
        if !t.is_finite() || t <= 0.0 {
            return None;
        }
        let p = &self.parameters;
        if self.kind == 0 {
            Some(p[0] * (t / p[1]).powf(1.5) * (p[1] + p[2]) / (t + p[2]))
        } else {
            *out_of_range |= t < p[3] || t > p[4];
            Some(p[0] * t.powf(1.5) / (t + p[1] * 10_f64.powf(-p[2] / t)))
        }
    }

    fn evaluate_with_derivative(&self, t: f64, out: &mut bool) -> Option<(f64, f64)> {
        // Reuse the value expression verbatim, including range/failure handling.
        let mu = self.evaluate(t, out)?;
        let p = &self.parameters;
        let logarithmic_derivative = if self.kind == 0 {
            1.5 / t - 1.0 / (t + p[2])
        } else {
            let term = p[1] * 10_f64.powf(-p[2] / t);
            let denominator = t + term;
            let derivative = 1.0 + term * std::f64::consts::LN_10 * p[2] / t.powi(2);
            1.5 / t - derivative / denominator
        };
        Some((mu, mu * logarithmic_derivative))
    }
}

#[pyclass]
struct LaminarKernel {
    analytic: bool,
    n: usize,
    pr: f64,
    viscosity: Viscosity,
    h: Vec<f64>,
    width: Vec<f64>,
    weights: Vec<[f64; 3]>,
}

struct Edge {
    xi: f64,
    beta: f64,
    heating: f64,
    temperature: f64,
    wall_ratio: Option<f64>,
    coefficient: f64,
}

impl LaminarKernel {
    fn residual(
        &self,
        state: &[f64],
        lagged: &[f64],
        edge: &Edge,
        result: &mut [f64],
        c: &[f64],
        conductivity: &[f64],
    ) {
        let n = self.n;
        let f = &state[..n];
        let u = &state[n..2 * n];
        let g = &state[2 * n..];
        result[0] = f[0];
        result[n] = u[0];
        result[2 * n - 1] = u[n - 1] - 1.0;
        result[3 * n - 1] = g[n - 1] - 1.0;
        for i in 1..n {
            result[i] = f[i] - f[i - 1] - self.h[i - 1] / 2.0 * (u[i] + u[i - 1]);
        }
        let scaled_xi = 2.0 * edge.xi;
        for i in 1..n - 1 {
            let w = self.weights[i - 1];
            let du = w[0] * u[i - 1] + w[1] * u[i] + w[2] * u[i + 1];
            let dg = w[0] * g[i - 1] + w[1] * g[i] + w[2] * g[i + 1];
            let convection = f[i] + scaled_xi * (edge.coefficient * f[i] + lagged[i]);
            let flux_right = (c[i + 1] + c[i]) / 2.0 * ((u[i + 1] - u[i]) / self.h[i]);
            let flux_left = (c[i] + c[i - 1]) / 2.0 * ((u[i] - u[i - 1]) / self.h[i - 1]);
            let kg_right =
                (conductivity[i + 1] + conductivity[i]) / 2.0 * ((g[i + 1] - g[i]) / self.h[i]);
            let kg_left =
                (conductivity[i] + conductivity[i - 1]) / 2.0 * ((g[i] - g[i - 1]) / self.h[i - 1]);
            result[n + i] = (flux_right - flux_left) / self.width[i - 1] + convection * du
                - scaled_xi * u[i] * (edge.coefficient * u[i] + lagged[n + i])
                + edge.beta * (g[i] - u[i].powi(2));
            result[2 * n + i] = (kg_right - kg_left) / self.width[i - 1]
                + convection * dg
                + edge.heating * c[i] * du.powi(2)
                - scaled_xi * u[i] * (edge.coefficient * g[i] + lagged[2 * n + i]);
        }
        result[2 * n] = if let Some(wall) = edge.wall_ratio {
            g[0] - wall
        } else {
            let face = (conductivity[0] + conductivity[1]) / 2.0;
            let slope = (u[1] - u[0]) / self.h[0];
            face * (g[1] - g[0]) / self.h[0] + self.h[0] / 2.0 * edge.heating * c[0] * slope.powi(2)
        };
    }

    fn iterate(
        &self,
        mut state: Vec<f64>,
        lagged: &[f64],
        edge: &Edge,
        iterations: usize,
        tolerance: f64,
    ) -> (Option<Vec<f64>>, bool) {
        let n = self.n;
        let m = 3 * n;
        let mut out = false;
        let mut residual = vec![0.0; m];
        let mut perturbed_residual = vec![0.0; m];
        let mut perturbed = vec![0.0; m];
        let mut c = vec![0.0; n];
        let mut perturbed_c = vec![0.0; n];
        let mut conductivity = vec![0.0; n];
        let mut perturbed_conductivity = vec![0.0; n];
        let mut step = vec![0.0; m];
        let mut matrix = BandMatrix::new(m);
        let mut rhs = vec![0.0; m];
        // Python's empty iteration range performs no viscosity evaluations.
        if iterations == 0 {
            return (None, out);
        }
        let Some(mu_edge) = self.viscosity.evaluate(edge.temperature, &mut out) else {
            return (None, out);
        };
        for _ in 0..iterations {
            for i in 0..n {
                let g = state[2 * n + i];
                let Some(mu) = self.viscosity.evaluate(edge.temperature * g, &mut out) else {
                    return (None, out);
                };
                c[i] = mu / mu_edge / g;
                conductivity[i] = c[i] / self.pr;
            }
            self.residual(&state, lagged, edge, &mut residual, &c, &conductivity);
            for i in 0..m {
                step[i] = f64::EPSILON.sqrt()
                    * if state[i].is_nan() {
                        f64::NAN
                    } else {
                        state[i].abs().max(1.0)
                    };
            }
            // Some scalar bands are absent from the block stencil and the LU
            // creates upper fill, so both must be cleared before reassembly.
            matrix.data.fill(0.0);
            for color in 0..3 {
                for variable in 0..3 {
                    perturbed.copy_from_slice(&state);
                    for col in (color..n).step_by(3) {
                        perturbed[variable * n + col] += step[variable * n + col];
                    }
                    let (viscosity, thermal) = if variable == 2 {
                        perturbed_c.copy_from_slice(&c);
                        perturbed_conductivity.copy_from_slice(&conductivity);
                        for col in (color..n).step_by(3) {
                            let g = perturbed[2 * n + col];
                            let Some(mu) = self.viscosity.evaluate(edge.temperature * g, &mut out)
                            else {
                                return (None, out);
                            };
                            perturbed_c[col] = mu / mu_edge / g;
                            perturbed_conductivity[col] = perturbed_c[col] / self.pr;
                        }
                        (&perturbed_c, &perturbed_conductivity)
                    } else {
                        (&c, &conductivity)
                    };
                    self.residual(
                        &perturbed,
                        lagged,
                        edge,
                        &mut perturbed_residual,
                        viscosity,
                        thermal,
                    );
                    for col in (color..n).step_by(3) {
                        for row in col.saturating_sub(1)..=(col + 1).min(n - 1) {
                            for eq in 0..3 {
                                let value = (perturbed_residual[eq * n + row]
                                    - residual[eq * n + row])
                                    / step[variable * n + col];
                                matrix.set(3 * row + eq, 3 * col + variable, value);
                            }
                        }
                    }
                }
            }
            for i in 0..n {
                for k in 0..3 {
                    rhs[3 * i + k] = -residual[k * n + i];
                }
            }
            if !matrix.solve(&mut rhs) || rhs.iter().any(|v| !v.is_finite()) {
                return (None, out);
            }
            let mut du: f64 = 0.0;
            let mut dg: f64 = 0.0;
            for i in 0..n {
                du = du.max(rhs[3 * i + 1].abs());
                let v = (rhs[3 * i + 2] / state[2 * n + i]).abs();
                if v.is_nan() {
                    dg = f64::NAN;
                } else if !dg.is_nan() {
                    dg = dg.max(v);
                }
            }
            // Python's built-in min/max retain the first argument on unordered comparisons.
            let d1 = 0.5 / if du < 1e-30 { 1e-30 } else { du };
            let d2 = 0.7 / if dg < 1e-30 { 1e-30 } else { dg };
            let mut damping = 1.0;
            if d1 < damping {
                damping = d1;
            }
            if d2 < damping {
                damping = d2;
            }
            if damping == 0.0 {
                return (None, out);
            }
            let size = if dg > du { dg } else { du };
            for i in 0..n {
                for k in 0..3 {
                    state[k * n + i] += damping * rhs[3 * i + k];
                }
            }
            if size < tolerance {
                return (Some(state), out);
            }
        }
        (None, out)
    }
}

// LAPACK-style column bands: the main diagonal is at row 2 * BAND.
// The top BAND rows hold the upper fill introduced by row pivoting.
struct BandMatrix {
    n: usize,
    data: Vec<f64>,
}
impl BandMatrix {
    fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * STRIDE],
        }
    }
    fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[col * STRIDE + 2 * BAND + row - col] = value;
    }
    fn solve(&mut self, b: &mut [f64]) -> bool {
        const DIAG: usize = 2 * BAND;
        let mut last_upper = 0;
        for j in 0..self.n {
            let km = BAND.min(self.n - 1 - j);
            let base = j * STRIDE + DIAG;
            let mut pivot = 0;
            for i in 1..=km {
                if self.data[base + i].abs() > self.data[base + pivot].abs() {
                    pivot = i;
                }
            }
            if self.data[base + pivot] == 0.0 {
                return false;
            }
            last_upper = last_upper.max((j + BAND + pivot).min(self.n - 1));
            if pivot != 0 {
                for col in j..=last_upper {
                    let row = col * STRIDE + DIAG - (col - j);
                    self.data.swap(row, row + pivot);
                }
                b.swap(j, j + pivot);
            }
            let diagonal = self.data[base];
            let mut factors = [0.0; BAND];
            for (factor, &value) in factors[..km]
                .iter_mut()
                .zip(&self.data[base + 1..base + 1 + km])
            {
                *factor = value / diagonal;
            }
            let bj = b[j];
            for (rhs, &factor) in b[j + 1..j + 1 + km].iter_mut().zip(&factors[..km]) {
                *rhs -= factor * bj;
            }
            for col in j + 1..=last_upper {
                let pivot_index = col * STRIDE + DIAG - (col - j);
                let upper = self.data[pivot_index];
                let column = &mut self.data[pivot_index + 1..pivot_index + 1 + km];
                for (value, &factor) in column.iter_mut().zip(&factors[..km]) {
                    *value -= factor * upper;
                }
            }
        }
        for j in (0..self.n).rev() {
            let diagonal = j * STRIDE + DIAG;
            b[j] /= self.data[diagonal];
            let value = b[j];
            let count = (2 * BAND).min(j);
            for (rhs, &upper) in b[j - count..j]
                .iter_mut()
                .zip(&self.data[diagonal - count..diagonal])
            {
                *rhs -= upper * value;
            }
        }
        true
    }
}

#[pymethods]
impl LaminarKernel {
    #[new]
    #[pyo3(signature = (eta, prandtl_number, viscosity_kind, viscosity_parameters, jacobian="analytic"))]
    fn new(
        eta: PyReadonlyArray1<'_, f64>,
        prandtl_number: f64,
        viscosity_kind: i32,
        viscosity_parameters: Vec<f64>,
        jacobian: &str,
    ) -> PyResult<Self> {
        let analytic = match jacobian {
            "analytic" => true,
            "finite-difference" => false,
            _ => {
                return Err(PyValueError::new_err(
                    "jacobian must be analytic or finite-difference",
                ))
            }
        };
        let eta: Vec<f64> = eta.as_array().iter().copied().collect();
        if eta.len() < 3
            || eta
                .windows(2)
                .any(|w| w[1].partial_cmp(&w[0]) != Some(std::cmp::Ordering::Greater))
        {
            return Err(PyValueError::new_err(
                "eta must contain at least three strictly increasing values",
            ));
        }
        if !matches!(
            (viscosity_kind, viscosity_parameters.len()),
            (0, 3) | (1, 5)
        ) {
            return Err(PyValueError::new_err(
                "invalid viscosity kind or parameter count",
            ));
        }
        let h: Vec<_> = eta.windows(2).map(|w| w[1] - w[0]).collect();
        let mut width = Vec::new();
        let mut weights = Vec::new();
        for pair in h.windows(2) {
            let (left, right) = (pair[0], pair[1]);
            width.push((left + right) / 2.0);
            weights.push([
                -right / (left * (left + right)),
                (right - left) / (left * right),
                left / (right * (left + right)),
            ]);
        }
        Ok(Self {
            analytic,
            n: eta.len(),
            pr: prandtl_number,
            viscosity: Viscosity {
                kind: viscosity_kind as u8,
                parameters: viscosity_parameters,
            },
            h,
            width,
            weights,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn newton<'py>(
        &self,
        py: Python<'py>,
        xi: f64,
        beta: f64,
        heating: f64,
        edge_temperature: f64,
        wall_ratio: Option<f64>,
        coefficient: f64,
        lagged: PyReadonlyArray2<'_, f64>,
        guess: PyReadonlyArray2<'_, f64>,
        max_iterations: i64,
        tolerance: f64,
    ) -> PyResult<(Option<Bound<'py, PyArray2<f64>>>, bool)> {
        let lagged = lagged.as_array();
        let guess = guess.as_array();
        if lagged.shape() != [3, self.n] || guess.shape() != [3, self.n] {
            return Err(PyValueError::new_err(
                "lagged and guess must have shape (3, N)",
            ));
        }
        let lagged: Vec<f64> = lagged.iter().copied().collect();
        let guess: Vec<f64> = guess.iter().copied().collect();
        let edge = Edge {
            xi,
            beta,
            heating,
            temperature: edge_temperature,
            wall_ratio,
            coefficient,
        };
        let (state, out) = py.detach(|| {
            let iterate = if self.analytic {
                Self::iterate_analytic
            } else {
                Self::iterate
            };
            iterate(
                self,
                guess,
                &lagged,
                &edge,
                max_iterations.max(0) as usize,
                tolerance,
            )
        });
        let array = state
            .map(|values| {
                Array2::from_shape_vec((3, self.n), values)
                    .map(|a| PyArray2::from_owned_array(py, a))
            })
            .transpose()
            .map_err(|_| PyValueError::new_err("invalid state shape"))?;
        Ok((array, out))
    }
}

#[pymodule]
fn aerophysics_native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<LaminarKernel>()?;
    module.add_class::<edge::NativeEdge>()?;
    module.add_function(wrap_pyfunction!(march::march, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_matrix(dense: &[Vec<f64>], expected: &[f64]) {
        let n = dense.len();
        let mut matrix = BandMatrix::new(n);
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            for (j, &value) in dense[i].iter().enumerate() {
                if i.abs_diff(j) <= BAND {
                    matrix.set(i, j, value);
                }
                rhs[i] += value * expected[j];
            }
        }
        let mut reference = dense.to_vec();
        let mut dense_rhs = rhs.clone();
        for col in 0..n {
            let pivot = (col..n)
                .max_by(|&a, &b| reference[a][col].abs().total_cmp(&reference[b][col].abs()))
                .expect("nonempty pivot range");
            reference.swap(col, pivot);
            dense_rhs.swap(col, pivot);
            for row in col + 1..n {
                let factor = reference[row][col] / reference[col][col];
                let (above, below) = reference.split_at_mut(row);
                for (value, &pivot_value) in
                    below[0][col + 1..].iter_mut().zip(&above[col][col + 1..])
                {
                    *value -= factor * pivot_value;
                }
                dense_rhs[row] -= factor * dense_rhs[col];
            }
        }
        for row in (0..n).rev() {
            for col in row + 1..n {
                dense_rhs[row] -= reference[row][col] * dense_rhs[col];
            }
            dense_rhs[row] /= reference[row][row];
        }
        assert!(matrix.solve(&mut rhs));
        for ((actual, reference), expected) in rhs.iter().zip(&dense_rhs).zip(expected) {
            assert!((actual - reference).abs() < 1e-10);
            assert!((actual - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn banded_random_and_pivoting() {
        let n = 31;
        let mut seed = 1123_u64;
        let mut dense = vec![vec![0.0; n]; n];
        for (i, row) in dense.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                if i.abs_diff(j) <= BAND {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    *value = ((seed >> 32) as f64) / (u32::MAX as f64) - 0.5;
                }
            }
            row[i] += 4.0;
        }
        let expected: Vec<_> = (0..n).map(|i| (i as f64).sin()).collect();
        check_matrix(&dense, &expected);
        dense[0][0] = 0.0;
        dense[1][0] = 2.0;
        check_matrix(&dense, &expected);
        // Pivot row BAND into row zero, carrying its original upper-band
        // entry to the furthest possible fill position (row 0, column 2*BAND).
        let mut full_fill = vec![vec![0.0; 2 * BAND + 2]; 2 * BAND + 2];
        for (i, row) in full_fill.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        full_fill[BAND][0] = 4.0;
        full_fill[BAND][2 * BAND] = 3.0;
        let solution: Vec<_> = (0..full_fill.len()).map(|i| i as f64 + 1.0).collect();
        check_matrix(&full_fill, &solution);
        // Exercise repeated pivots and fill on non-diagonally-dominant bands.
        for _ in 0..20 {
            for (i, row) in dense.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    *value = if i.abs_diff(j) <= BAND {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) as f64) / (u32::MAX as f64) - 0.5
                    } else {
                        0.0
                    };
                }
            }
            check_matrix(&dense, &expected);
        }
        let mut singular = BandMatrix::new(5);
        assert!(!singular.solve(&mut [1.0; 5]));
    }

    #[test]
    fn viscosities() {
        let s = Viscosity {
            kind: 0,
            parameters: vec![1.7894e-5, 288.15, 110.4],
        };
        let k = Viscosity {
            kind: 1,
            parameters: vec![1.488e-6, 122.1, 5.0, 79.0, 1845.0],
        };
        let mut out = false;
        assert!((s.evaluate(288.15, &mut out).expect("positive") - 1.7894e-5).abs() < 1e-18);
        let expected = 1.488e-6 * 300_f64.powf(1.5) / (300.0 + 122.1 * 10_f64.powf(-5.0 / 300.0));
        assert!((k.evaluate(300.0, &mut out).expect("positive") - expected).abs() < 1e-18);
        assert!(!out);
        assert!(k.evaluate(1900.0, &mut out).is_some() && out);
        assert!(k.evaluate(0.0, &mut out).is_none());
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(s.evaluate(bad, &mut out).is_none());
            assert!(k.evaluate(bad, &mut out).is_none());
        }
    }
}

#[cfg(test)]
mod derivative_tests {
    use super::*;
    #[test]
    fn viscosity_derivatives_values_and_failures() {
        for viscosity in [
            Viscosity {
                kind: 0,
                parameters: vec![1.7894e-5, 288.15, 110.4],
            },
            Viscosity {
                kind: 1,
                parameters: vec![1.488e-6, 122.1, 5.0, 79.0, 1845.0],
            },
        ] {
            for t in [60.0, 79.0, 288.15, 500.0, 1845.0, 2000.0] {
                let mut out = false;
                let mut derivative_out = false;
                let mu = viscosity.evaluate(t, &mut out).unwrap();
                let (value, derivative) = viscosity
                    .evaluate_with_derivative(t, &mut derivative_out)
                    .unwrap();
                assert_eq!(mu.to_bits(), value.to_bits());
                assert_eq!(out, derivative_out);
                let h = t * 1e-5;
                let central = (viscosity.evaluate(t + h, &mut out).unwrap()
                    - viscosity.evaluate(t - h, &mut out).unwrap())
                    / (2.0 * h);
                assert!((derivative - central).abs() < 1e-9 * derivative.abs());
            }
            for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                for initial in [false, true] {
                    let mut out = initial;
                    assert!(viscosity.evaluate_with_derivative(bad, &mut out).is_none());
                    assert_eq!(out, initial);
                }
            }
        }
    }
}
