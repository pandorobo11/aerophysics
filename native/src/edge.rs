//! Scalar PPoly evaluation and laminar isentropic edge states (SI units).
use crate::{Edge, Viscosity};
use numpy::ndarray::Array2;
use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub(crate) struct PPoly {
    x: Vec<f64>,
    c: Vec<f64>,
    order: usize,
}

impl PPoly {
    fn new(x: PyReadonlyArray1<'_, f64>, c: PyReadonlyArray2<'_, f64>) -> PyResult<Self> {
        let x: Vec<_> = x.as_array().iter().copied().collect();
        let c = c.as_array();
        if x.len() < 2
            || c.shape()[1] != x.len() - 1
            || c.shape()[0] == 0
            || x.iter().any(|v| !v.is_finite())
            || x.windows(2).any(|w| w[1] <= w[0])
        {
            return Err(PyValueError::new_err(
                "invalid PPoly breakpoints or coefficients",
            ));
        }
        Ok(Self {
            x,
            c: c.iter().copied().collect(),
            order: c.shape()[0],
        })
    }

    fn interval(&self, s: f64) -> usize {
        self.x
            .partition_point(|&v| v <= s)
            .saturating_sub(1)
            .min(self.x.len() - 2)
    }

    fn evaluate(&self, s: f64) -> f64 {
        let i = self.interval(s);
        let t = s - self.x[i];
        let mut result = 0.0;
        let mut z = 1.0;
        // Match scipy.interpolate._ppoly.evaluate_poly1, including sum order.
        for kp in 0..self.order {
            result += self.c[(self.order - 1 - kp) * (self.x.len() - 1) + i] * z;
            z *= t;
        }
        result
    }
}

pub(crate) struct State {
    pub equation: Edge,
    pub velocity: f64,
    pub pressure: f64,
}

impl State {
    pub fn values(&self) -> [f64; 7] {
        let e = &self.equation;
        [
            e.xi,
            self.velocity,
            e.beta,
            e.temperature,
            self.pressure,
            e.heating,
            e.wall_ratio.unwrap_or(f64::NAN),
        ]
    }
}

#[pyclass]
pub(crate) struct NativeEdge {
    log_pressure: PPoly,
    log_pressure_slope: PPoly,
    wall_temperature: Option<PPoly>,
    panels: Vec<f64>,
    xi_at_panels: Vec<f64>,
    nodes: Vec<f64>,
    weights: Vec<f64>,
    total_temperature: f64,
    total_pressure: f64,
    gamma: f64,
    cp: f64,
    gas_constant: f64,
    viscosity: Viscosity,
}

impl NativeEdge {
    fn thermodynamics(&self, s: f64) -> (f64, f64, f64, f64) {
        let pressure = self.log_pressure.evaluate(s).exp();
        let temperature = self.total_temperature
            * (pressure / self.total_pressure).powf((self.gamma - 1.0) / self.gamma);
        let velocity = (2.0 * self.cp * (self.total_temperature - temperature)).sqrt();
        let density = pressure / (self.gas_constant * temperature);
        (pressure, temperature, velocity, density)
    }

    fn flux(&self, s: f64, out: &mut bool) -> Result<f64, &'static str> {
        let (_, t, u, rho) = self.thermodynamics(s);
        Ok(rho * viscosity(&self.viscosity, t, out)? * u)
    }

    fn panel(&self, s: f64) -> usize {
        self.panels
            .partition_point(|&v| v <= s)
            .saturating_sub(1)
            .min(self.panels.len() - 1)
    }

    fn xi(&self, s: f64, out: &mut bool) -> Result<f64, &'static str> {
        let i = self.panel(s);
        let left = self.panels[i];
        let mut increment = 0.0;
        if s > left {
            let half = (s - left) / 2.0;
            for (&node, &weight) in self.nodes.iter().zip(&self.weights) {
                increment += weight * self.flux(left + half * (node + 1.0), out)?;
            }
            increment *= half;
        }
        Ok(self.xi_at_panels[i] + increment)
    }

    pub fn state(&self, s: f64, out: &mut bool) -> Result<State, &'static str> {
        let (pressure, temperature, velocity, density) = self.thermodynamics(s);
        let viscosity = viscosity(&self.viscosity, temperature, out)?;
        let xi = self.xi(s, out)?;
        let log_velocity_slope =
            -pressure * self.log_pressure_slope.evaluate(s) / (density * velocity.powi(2));
        Ok(State {
            equation: Edge {
                xi,
                beta: 2.0 * xi * log_velocity_slope / (density * viscosity * velocity),
                heating: velocity.powi(2) / (self.cp * temperature),
                temperature,
                wall_ratio: self
                    .wall_temperature
                    .as_ref()
                    .map(|p| p.evaluate(s) / temperature),
                coefficient: 0.0,
            },
            velocity,
            pressure,
        })
    }
}

fn viscosity(model: &Viscosity, t: f64, out: &mut bool) -> Result<f64, &'static str> {
    model.evaluate(t, out).ok_or(if t.is_finite() {
        "temperature must be greater than zero"
    } else {
        "temperature must contain only finite values"
    })
}

#[pymethods]
impl NativeEdge {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        log_pressure: (PyReadonlyArray1<'_, f64>, PyReadonlyArray2<'_, f64>),
        log_pressure_slope: (PyReadonlyArray1<'_, f64>, PyReadonlyArray2<'_, f64>),
        wall_temperature: Option<(PyReadonlyArray1<'_, f64>, PyReadonlyArray2<'_, f64>)>,
        panels: PyReadonlyArray1<'_, f64>,
        xi_at_panels: PyReadonlyArray1<'_, f64>,
        nodes: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
        gas: (f64, f64, f64, f64, f64),
        viscosity_kind: i32,
        viscosity_parameters: Vec<f64>,
    ) -> PyResult<Self> {
        let panels: Vec<_> = panels.as_array().iter().copied().collect();
        let xi_at_panels: Vec<_> = xi_at_panels.as_array().iter().copied().collect();
        let nodes: Vec<_> = nodes.as_array().iter().copied().collect();
        let weights: Vec<_> = weights.as_array().iter().copied().collect();
        if panels.is_empty()
            || panels.len() != xi_at_panels.len()
            || panels.iter().any(|v| !v.is_finite())
            || panels.windows(2).any(|w| w[1] <= w[0])
            || nodes.len() != 5
            || weights.len() != 5
        {
            return Err(PyValueError::new_err(
                "invalid integration panels or quadrature",
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
        Ok(Self {
            log_pressure: PPoly::new(log_pressure.0, log_pressure.1)?,
            log_pressure_slope: PPoly::new(log_pressure_slope.0, log_pressure_slope.1)?,
            wall_temperature: wall_temperature
                .map(|(x, c)| PPoly::new(x, c))
                .transpose()?,
            panels,
            xi_at_panels,
            nodes,
            weights,
            total_temperature: gas.0,
            total_pressure: gas.1,
            gamma: gas.2,
            cp: gas.3,
            gas_constant: gas.4,
            viscosity: Viscosity {
                kind: viscosity_kind as u8,
                parameters: viscosity_parameters,
            },
        })
    }

    /// Columns: xi, velocity, beta, temperature, pressure, heating, wall_ratio.
    fn states<'py>(
        &self,
        py: Python<'py>,
        distances: Vec<f64>,
    ) -> PyResult<(Bound<'py, PyArray2<f64>>, bool)> {
        let mut out = false;
        let values = py
            .detach(|| {
                distances
                    .iter()
                    .map(|&s| self.state(s, &mut out).map(|v| v.values()))
                    .collect::<Result<Vec<_>, _>>()
            })
            .map_err(PyValueError::new_err)?;
        let array =
            Array2::from_shape_vec((distances.len(), 7), values.into_iter().flatten().collect())
                .unwrap();
        Ok((PyArray2::from_owned_array(py, array), out))
    }

    /// Test hook: evaluate the supplied PPoly coefficients without thermodynamics.
    fn polynomials(&self, distances: Vec<f64>) -> Vec<(f64, f64, Option<f64>)> {
        distances
            .iter()
            .map(|&s| {
                (
                    self.log_pressure.evaluate(s),
                    self.log_pressure_slope.evaluate(s),
                    self.wall_temperature.as_ref().map(|p| p.evaluate(s)),
                )
            })
            .collect()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    #[test]
    fn intervals_and_power_sum() {
        let p = PPoly {
            x: vec![0.0, 1.0, 3.0],
            c: vec![2.0, 4.0, 3.0, 5.0, 7.0, 11.0],
            order: 3,
        };
        for (s, i, value) in [
            (-1.0, 0, 6.0),
            (0.0, 0, 7.0),
            (0.5, 0, 9.0),
            (1.0, 1, 11.0),
            (3.0, 1, 37.0),
            (4.0, 1, 62.0),
        ] {
            assert_eq!(p.interval(s), i);
            assert_eq!(p.evaluate(s), value);
        }
    }

    pub(crate) fn flat_edge() -> NativeEdge {
        let constant = |value| PPoly {
            x: vec![0.0, 1.0],
            c: vec![value],
            order: 1,
        };
        let mut edge = NativeEdge {
            log_pressure: constant(100000_f64.ln()),
            log_pressure_slope: constant(0.0),
            wall_temperature: None,
            panels: vec![0.0, 0.5, 1.0],
            xi_at_panels: vec![0.0; 3],
            nodes: vec![0.0; 5],
            weights: vec![0.4; 5],
            total_temperature: 300.0,
            total_pressure: 100007.0,
            gamma: 1.4,
            cp: 1004.675,
            gas_constant: 287.05,
            viscosity: Viscosity {
                kind: 0,
                parameters: vec![1.7894e-5, 288.15, 110.4],
            },
        };
        let flux = edge.flux(0.0, &mut false).unwrap();
        edge.xi_at_panels = vec![0.0, 0.5 * flux, flux];
        edge
    }

    #[test]
    fn panel_lookup_and_partial_integral() {
        let edge = flat_edge();
        let flux = edge.flux(0.0, &mut false).unwrap();
        for (s, i) in [(-1.0, 0), (0.0, 0), (0.25, 0), (0.5, 1), (1.0, 2), (1.5, 2)] {
            assert_eq!(edge.panel(s), i);
            assert!((edge.xi(s, &mut false).unwrap() - s.max(0.0) * flux).abs() < 1e-18);
        }
    }
}
