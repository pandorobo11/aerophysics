use crate::{block::Blocks, BandMatrix, Edge, LaminarKernel};

impl LaminarKernel {
    pub(crate) fn iterate_analytic(
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
        let mut c = vec![0.0; n];
        let mut conductivity = vec![0.0; n];
        let mut dc = vec![0.0; n];
        let mut blocks = Blocks::new(n);
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
                let Some((mu, derivative)) = self
                    .viscosity
                    .evaluate_with_derivative(edge.temperature * g, &mut out)
                else {
                    return (None, out);
                };
                c[i] = mu / mu_edge / g;
                conductivity[i] = c[i] / self.pr;
                dc[i] = (edge.temperature * derivative / mu_edge - c[i]) / g;
            }
            self.residual(&state, lagged, edge, &mut residual, &c, &conductivity);
            self.analytic_jacobian(&state, lagged, edge, &c, &dc, &mut blocks);
            for i in 0..n {
                for k in 0..3 {
                    rhs[3 * i + k] = -residual[k * n + i];
                }
            }
            if !blocks.solve(&mut rhs, &mut matrix) {
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

impl LaminarKernel {
    #[allow(clippy::too_many_arguments)]
    fn analytic_jacobian(
        &self,
        state: &[f64],
        lagged: &[f64],
        edge: &Edge,
        c: &[f64],
        dc: &[f64],
        blocks: &mut Blocks,
    ) {
        let n = self.n;
        let f = &state[..n];
        let u = &state[n..2 * n];
        let g = &state[2 * n..];
        // Write full blocks (including structural zeros) on every assembly.
        blocks.lower[0] = [[0.0; 3]; 3];
        blocks.diagonal[0] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0; 3]];
        blocks.upper[0] = [[0.0; 3]; 3];
        let outer_half = self.h[n - 2] / 2.0;
        blocks.lower[n - 1] = [[-1.0, -outer_half, 0.0], [0.0; 3], [0.0; 3]];
        blocks.diagonal[n - 1] = [[1.0, -outer_half, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        blocks.upper[n - 1] = [[0.0; 3]; 3];
        let sx = 2.0 * edge.xi;
        for i in 1..n - 1 {
            let w = self.weights[i - 1];
            let du = w[0] * u[i - 1] + w[1] * u[i] + w[2] * u[i + 1];
            let dg = w[0] * g[i - 1] + w[1] * g[i] + w[2] * g[i + 1];
            let convection = f[i] + sx * (edge.coefficient * f[i] + lagged[i]);
            let width = self.width[i - 1];
            let left = (c[i] + c[i - 1]) / 2.0 / self.h[i - 1] / width;
            let right = (c[i + 1] + c[i]) / 2.0 / self.h[i] / width;
            let su_left = (u[i] - u[i - 1]) / self.h[i - 1] / 2.0 / width;
            let su_right = (u[i + 1] - u[i]) / self.h[i] / 2.0 / width;
            let sg_left = (g[i] - g[i - 1]) / self.h[i - 1] / 2.0 / width / self.pr;
            let sg_right = (g[i + 1] - g[i]) / self.h[i] / 2.0 / width / self.pr;
            let half = self.h[i - 1] / 2.0;
            let mut a = [[-1.0, -half, 0.0], [0.0; 3], [0.0; 3]];
            let mut b = [[1.0, -half, 0.0], [0.0; 3], [0.0; 3]];
            let mut upper = [[0.0; 3]; 3];
            a[1][1] = left + convection * w[0];
            b[1][1] = -left - right + convection * w[1]
                - sx * (2.0 * edge.coefficient * u[i] + lagged[n + i])
                - 2.0 * edge.beta * u[i];
            upper[1][1] = right + convection * w[2];
            b[1][0] = (1.0 + sx * edge.coefficient) * du;
            a[1][2] = -dc[i - 1] * su_left;
            b[1][2] = dc[i] * (su_right - su_left) + edge.beta;
            upper[1][2] = dc[i + 1] * su_right;
            b[2][0] = (1.0 + sx * edge.coefficient) * dg;
            let dissipation = 2.0 * edge.heating * c[i] * du;
            a[2][1] = dissipation * w[0];
            b[2][1] = dissipation * w[1] - sx * (edge.coefficient * g[i] + lagged[2 * n + i]);
            upper[2][1] = dissipation * w[2];
            a[2][2] = left / self.pr - dc[i - 1] * sg_left + convection * w[0];
            b[2][2] = (-left - right) / self.pr
                + dc[i] * (sg_right - sg_left)
                + convection * w[1]
                + edge.heating * dc[i] * du.powi(2)
                - sx * u[i] * edge.coefficient;
            upper[2][2] = right / self.pr + dc[i + 1] * sg_right + convection * w[2];
            blocks.lower[i] = a;
            blocks.diagonal[i] = b;
            blocks.upper[i] = upper;
        }
        if edge.wall_ratio.is_some() {
            blocks.diagonal[0][2][2] = 1.0;
        } else {
            let h = self.h[0];
            let slope = (u[1] - u[0]) / h;
            let face = (c[0] / self.pr + c[1] / self.pr) / 2.0;
            let dg = (g[1] - g[0]) / h;
            blocks.diagonal[0][2][1] = -edge.heating * c[0] * slope;
            blocks.upper[0][2][1] = edge.heating * c[0] * slope;
            blocks.diagonal[0][2][2] = dc[0] / self.pr / 2.0 * dg - face / h
                + h / 2.0 * edge.heating * dc[0] * slope.powi(2);
            blocks.upper[0][2][2] = dc[1] / self.pr / 2.0 * dg + face / h;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Viscosity;

    #[test]
    fn reassembly_overwrites_every_block_entry() {
        let kernel = LaminarKernel {
            analytic: true,
            n: 3,
            pr: 0.72,
            viscosity: Viscosity {
                kind: 0,
                parameters: vec![1.7894e-5, 288.15, 110.4],
            },
            h: vec![0.1, 0.2],
            width: vec![0.15],
            weights: vec![[-1.0, 0.0, 1.0]],
        };
        let state = [0.0, 0.01, 0.05, 0.0, 0.2, 1.0, 1.4, 1.2, 1.0];
        for wall_ratio in [None, Some(1.4)] {
            let edge = Edge {
                xi: 0.4,
                beta: 0.2,
                heating: 1.3,
                temperature: 300.0,
                wall_ratio,
                coefficient: 0.8,
            };
            let mut reference = Blocks::new(3);
            let mut poisoned = Blocks::new(3);
            for values in [
                &mut poisoned.lower,
                &mut poisoned.diagonal,
                &mut poisoned.upper,
            ] {
                values.fill([[f64::NAN; 3]; 3]);
            }
            for blocks in [&mut reference, &mut poisoned] {
                kernel.analytic_jacobian(
                    &state,
                    &[0.1; 9],
                    &edge,
                    &[0.9, 0.95, 1.0],
                    &[-0.1; 3],
                    blocks,
                );
            }
            assert_eq!(reference.lower, poisoned.lower);
            assert_eq!(reference.diagonal, poisoned.diagonal);
            assert_eq!(reference.upper, poisoned.upper);
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)] // Independently perturb each dense column.
    fn jacobian_matches_forward_difference() {
        for kind in 0..2 {
            for wall_ratio in [None, Some(1.4)] {
                for perturbation in [0.0, 0.01, -0.02] {
                    let n = 41;
                    let eta: Vec<_> = (0..n).map(|i| 0.03 * (i as f64).powf(1.5)).collect();
                    let h: Vec<_> = eta.windows(2).map(|v| v[1] - v[0]).collect();
                    let kernel = LaminarKernel {
                        analytic: true,
                        n,
                        pr: 0.72,
                        viscosity: Viscosity {
                            kind,
                            parameters: if kind == 0 {
                                vec![1.7894e-5, 288.15, 110.4]
                            } else {
                                vec![1.488e-6, 122.1, 5.0, 79.0, 1845.0]
                            },
                        },
                        width: h.windows(2).map(|v| (v[0] + v[1]) / 2.0).collect(),
                        weights: h
                            .windows(2)
                            .map(|v| {
                                let (l, r) = (v[0], v[1]);
                                [-r / (l * (l + r)), (r - l) / (l * r), l / (r * (l + r))]
                            })
                            .collect(),
                        h,
                    };
                    let edge = Edge {
                        xi: 0.37,
                        beta: 0.23,
                        heating: 1.8,
                        temperature: 290.0,
                        wall_ratio,
                        coefficient: 1.3,
                    };
                    let mut state = vec![0.0; 3 * n];
                    let recovery = 1.0 + kernel.pr.sqrt() * edge.heating / 2.0;
                    let wall = wall_ratio.unwrap_or(recovery);
                    for i in 0..n {
                        let u = (0.5 * eta[i]).tanh();
                        state[n + i] = u + perturbation * (i as f64).sin();
                        state[2 * n + i] = wall
                            + (recovery - wall) * u
                            + (1.0 - recovery) * u * u
                            + perturbation * (i as f64).cos();
                        if i > 0 {
                            state[i] = state[i - 1]
                                + kernel.h[i - 1] / 2.0 * (state[n + i] + state[n + i - 1]);
                        }
                    }
                    let lagged: Vec<_> = state.iter().map(|v| -0.8 * v + 0.03).collect();
                    let mut out = false;
                    let mu_edge = kernel
                        .viscosity
                        .evaluate(edge.temperature, &mut out)
                        .unwrap();
                    let mut c = vec![0.0; n];
                    let mut dc = c.clone();
                    let mut thermal = c.clone();
                    for i in 0..n {
                        let g = state[2 * n + i];
                        let (mu, dmu) = kernel
                            .viscosity
                            .evaluate_with_derivative(edge.temperature * g, &mut out)
                            .unwrap();
                        c[i] = mu / mu_edge / g;
                        dc[i] = (edge.temperature * dmu / mu_edge - c[i]) / g;
                        thermal[i] = c[i] / kernel.pr;
                    }
                    let mut blocks = Blocks::new(n);
                    kernel.analytic_jacobian(&state, &lagged, &edge, &c, &dc, &mut blocks);
                    let mut base = vec![0.0; 3 * n];
                    kernel.residual(&state, &lagged, &edge, &mut base, &c, &thermal);
                    let mut dense = vec![vec![0.0; 3 * n]; 3 * n];
                    for col in 0..3 * n {
                        let node = col / 3;
                        let variable = col % 3;
                        let index = variable * n + node;
                        let step = f64::EPSILON.sqrt() * state[index].abs().max(1.0);
                        let mut trial = state.clone();
                        trial[index] += step;
                        let mut pc = c.clone();
                        let mut pk = thermal.clone();
                        if variable == 2 {
                            let g = trial[index];
                            pc[node] = kernel
                                .viscosity
                                .evaluate(edge.temperature * g, &mut out)
                                .unwrap()
                                / mu_edge
                                / g;
                            pk[node] = pc[node] / kernel.pr;
                        }
                        let mut residual = vec![0.0; 3 * n];
                        kernel.residual(&trial, &lagged, &edge, &mut residual, &pc, &pk);
                        for row in 0..3 * n {
                            dense[row][col] = (residual[(row % 3) * n + row / 3]
                                - base[(row % 3) * n + row / 3])
                                / step;
                        }
                    }
                    for (row, entries) in dense.iter().enumerate() {
                        let scale = entries.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
                        for (col, &fd) in entries.iter().enumerate() {
                            let (i, j, r, v) = (row / 3, col / 3, row % 3, col % 3);
                            let exact = if i == j {
                                blocks.diagonal[i][r][v]
                            } else if i == j + 1 {
                                blocks.lower[i][r][v]
                            } else if j == i + 1 {
                                blocks.upper[i][r][v]
                            } else {
                                0.0
                            };
                            // Forward differences use sqrt(eps): O(step) truncation
                            // plus cancellation in flux differences. 1e-5 of the
                            // row scale allows this while checking every entry,
                            // including all off-stencil zeros and wall terms.
                            assert!((exact-fd).abs() <= 1e-5*scale, "kind={kind}, wall={wall_ratio:?}, row={row}, col={col}, exact={exact}, fd={fd}, scale={scale}");
                        }
                    }
                }
            }
        }
    }
}
