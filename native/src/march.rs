//! The Python laminar march, preserving its step decisions and BDF operation order.
use crate::edge::{NativeEdge, State};
use crate::{Edge, LaminarKernel};
use numpy::ndarray::{Array2, Array3};
use numpy::{PyArray1, PyArray2, PyArray3, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

struct Controls {
    first_step: f64,
    minimum_step: f64,
    growth: f64,
    max_ratio: f64,
    merge: f64,
    relative_step: f64,
    velocity_change: f64,
    shear_change: f64,
    heat_change: f64,
    heat_floor: f64,
    separation_fraction: f64,
    edge_tolerance: f64,
    probe: usize,
    iterations: usize,
    tolerance: f64,
}

struct Station {
    xi: f64,
    state: Vec<f64>,
    slopes: [f64; 2],
}

struct Output {
    status: &'static str,
    candidate: f64,
    positions: Vec<f64>,
    states: Vec<f64>,
    slopes: Vec<f64>,
    edges: Vec<f64>,
    accepted_positions: Vec<f64>,
    accepted_shear: Vec<f64>,
    out_of_range: bool,
}

fn solve(
    kernel: &LaminarKernel,
    edge: &Edge,
    guess: Vec<f64>,
    lagged: &[f64],
    wall_weights: &[f64],
    controls: &Controls,
    out: &mut bool,
) -> Option<Station> {
    let iterate = if kernel.analytic {
        LaminarKernel::iterate_analytic
    } else {
        LaminarKernel::iterate
    };
    let (state, outside) = iterate(
        kernel,
        guess,
        lagged,
        edge,
        controls.iterations,
        controls.tolerance,
    );
    *out |= outside;
    state.map(|state| {
        let mut slopes = [0.0; 2];
        for (k, slope) in slopes.iter_mut().enumerate() {
            for i in 0..3 {
                *slope += wall_weights[i] * state[(k + 1) * kernel.n + i];
            }
        }
        Station {
            xi: edge.xi,
            state,
            slopes,
        }
    })
}

fn starting_guess(kernel: &LaminarKernel, eta: &[f64], edge: &Edge) -> Vec<f64> {
    let n = kernel.n;
    let mut state = vec![0.0; 3 * n];
    let recovery = 1.0 + kernel.pr.sqrt() * edge.heating / 2.0;
    let wall = edge.wall_ratio.unwrap_or(recovery);
    for i in 0..n {
        let u = (0.5 * eta[i]).tanh();
        state[n + i] = u;
        state[2 * n + i] = wall + (recovery - wall) * u + (1.0 - recovery) * u.powi(2);
        if i > 0 {
            state[i] = state[i - 1] + kernel.h[i - 1] * (u + state[n + i - 1]) / 2.0;
        }
    }
    state
}

fn stencil(xi: f64, last: &Station, older: Option<&Station>) -> (f64, Vec<f64>) {
    let first_gap = xi - last.xi;
    if let Some(older) = older {
        let second_gap = xi - older.xi;
        let coefficient = 1.0 / first_gap + 1.0 / second_gap;
        let a = -second_gap / (first_gap * (second_gap - first_gap));
        let b = first_gap / (second_gap * (second_gap - first_gap));
        (
            coefficient,
            last.state
                .iter()
                .zip(&older.state)
                .map(|(&v, &o)| a * v + b * o)
                .collect(),
        )
    } else {
        let coefficient = 1.0 / first_gap;
        (
            coefficient,
            last.state.iter().map(|&v| -coefficient * v).collect(),
        )
    }
}

fn resolved(previous: &Station, station: &Station, edge: &Edge, c: &Controls) -> bool {
    let shear = station.slopes[0];
    if shear <= 0.0 || (shear / previous.slopes[0]).ln().abs() > c.shear_change {
        return false;
    }
    if edge.wall_ratio.is_none() {
        return true;
    }
    let (slope, earlier) = (station.slopes[1], previous.slopes[1]);
    let floor = c.heat_floor * shear * (1.0 + edge.heating / 2.0);
    // Python max keeps its first value on unordered comparisons.
    let mut scale = slope.abs();
    if earlier.abs() > scale {
        scale = earlier.abs();
    }
    if floor > scale {
        scale = floor;
    }
    (slope - earlier).abs() <= c.heat_change * scale
}

fn run(
    kernel: &LaminarKernel,
    edge_at: &NativeEdge,
    targets: &[f64],
    eta: &[f64],
    wall_weights: &[f64],
    c: &Controls,
) -> Result<Output, &'static str> {
    let mut output = Output {
        status: "ok",
        candidate: 0.0,
        positions: vec![],
        states: vec![],
        slopes: vec![],
        edges: vec![],
        accepted_positions: vec![],
        accepted_shear: vec![],
        out_of_range: false,
    };
    let end = targets[targets.len() - 1];
    let start_edge = edge_at.state(0.0, &mut output.out_of_range)?;
    if start_edge.equation.xi != 0.0 {
        output.status = "invalid_start";
        return Ok(output);
    }
    let start = solve(
        kernel,
        &start_edge.equation,
        starting_guess(kernel, eta, &start_edge.equation),
        &vec![0.0; 3 * kernel.n],
        wall_weights,
        c,
        &mut output.out_of_range,
    );
    let Some(mut last) = start else {
        output.status = "leading_edge_failed";
        return Ok(output);
    };
    let mut older = None;
    let mut position = 0.0;
    let mut velocity = start_edge.velocity;
    let mut step = c.first_step * targets[0].min(end);
    let mut last_step = step;
    let minimum_step = c.minimum_step * end;
    output.accepted_positions.push(0.0);
    output.accepted_shear.push(last.slopes[0]);
    let mut latest_edge: State = start_edge;
    for &target in targets {
        while position < target {
            let mut trial = step.min(c.max_ratio * last_step);
            if position > 0.0 {
                trial = trial.min(c.relative_step * position);
            }
            let remaining = target - position;
            if remaining <= c.merge * trial {
                trial = remaining;
            } else if remaining < 2.0 * trial {
                trial = remaining / 2.0;
            }
            let candidate = if trial == remaining {
                target
            } else {
                position + trial
            };
            output.candidate = candidate;
            let mut edge = edge_at.state(candidate, &mut output.out_of_range)?;
            let mut station = None;
            if (edge.velocity / velocity).ln().abs() <= c.velocity_change {
                let (coefficient, lagged) = stencil(edge.equation.xi, &last, older.as_ref());
                edge.equation.coefficient = coefficient;
                let mut guess = last.state.clone();
                if let Some(wall) = edge.equation.wall_ratio {
                    guess[2 * kernel.n] = wall;
                }
                station = solve(
                    kernel,
                    &edge.equation,
                    guess,
                    &lagged,
                    wall_weights,
                    c,
                    &mut output.out_of_range,
                );
            }
            if station
                .as_ref()
                .is_some_and(|s| !resolved(&last, s, &edge.equation, c))
            {
                station = None;
            }
            let Some(station) = station else {
                step = trial / 2.0;
                if step >= minimum_step {
                    continue;
                }
                let peak = output
                    .accepted_shear
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                output.status =
                    if edge.equation.beta < 0.0 && last.slopes[0] < c.separation_fraction * peak {
                        "separated"
                    } else {
                        "not_converged"
                    };
                return Ok(output);
            };
            let n = kernel.n;
            let outer = (station.state[2 * n - c.probe] - 1.0).abs()
                + (station.state[3 * n - c.probe] - 1.0).abs();
            if outer > c.edge_tolerance {
                output.status = "outer_boundary";
                return Ok(output);
            }
            older = Some(last);
            last = station;
            position = candidate;
            velocity = edge.velocity;
            last_step = trial;
            step = trial * c.growth;
            output.accepted_positions.push(position);
            output.accepted_shear.push(last.slopes[0]);
            latest_edge = edge;
        }
        output.positions.push(position);
        output.states.extend_from_slice(&last.state);
        output.slopes.extend_from_slice(&last.slopes);
        output.edges.extend_from_slice(&latest_edge.values());
    }
    Ok(output)
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub(crate) fn march<'py>(
    py: Python<'py>,
    kernel: &LaminarKernel,
    edge: &NativeEdge,
    targets: PyReadonlyArray1<'_, f64>,
    eta: PyReadonlyArray1<'_, f64>,
    wall_weights: PyReadonlyArray1<'_, f64>,
    constants: Vec<f64>,
    probe: usize,
    iterations: i64,
    tolerance: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let targets: Vec<_> = targets.as_array().iter().copied().collect();
    let eta: Vec<_> = eta.as_array().iter().copied().collect();
    let wall_weights: Vec<_> = wall_weights.as_array().iter().copied().collect();
    if targets.is_empty()
        || targets.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || targets.windows(2).any(|w| w[1] <= w[0])
        || eta.len() != kernel.n
        || wall_weights.len() != 3
        || constants.len() != 12
        || probe == 0
        || probe > kernel.n
    {
        return Err(PyValueError::new_err("invalid march arrays or controls"));
    }
    let c = Controls {
        first_step: constants[0],
        minimum_step: constants[1],
        growth: constants[2],
        max_ratio: constants[3],
        merge: constants[4],
        relative_step: constants[5],
        velocity_change: constants[6],
        shear_change: constants[7],
        heat_change: constants[8],
        heat_floor: constants[9],
        separation_fraction: constants[10],
        edge_tolerance: constants[11],
        probe,
        iterations: iterations.max(0) as usize,
        tolerance,
    };
    let result = py
        .detach(|| run(kernel, edge, &targets, &eta, &wall_weights, &c))
        .map_err(PyValueError::new_err)?;
    let count = result.positions.len();
    let dict = PyDict::new(py);
    dict.set_item("status", result.status)?;
    dict.set_item("candidate", result.candidate)?;
    dict.set_item("out_of_range", result.out_of_range)?;
    dict.set_item("positions", PyArray1::from_vec(py, result.positions))?;
    dict.set_item(
        "states",
        PyArray3::from_owned_array(
            py,
            Array3::from_shape_vec((count, 3, kernel.n), result.states).unwrap(),
        ),
    )?;
    dict.set_item(
        "slopes",
        PyArray2::from_owned_array(
            py,
            Array2::from_shape_vec((count, 2), result.slopes).unwrap(),
        ),
    )?;
    dict.set_item(
        "edges",
        PyArray2::from_owned_array(
            py,
            Array2::from_shape_vec((count, 7), result.edges).unwrap(),
        ),
    )?;
    dict.set_item("accepted_positions", result.accepted_positions)?;
    dict.set_item("accepted_shear", result.accepted_shear)?;
    Ok(dict)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn short_flat_plate_blasius() {
        // f uses eta = y sqrt(U/(2 nu x)), so f''(0)/sqrt(2) = 0.332057.
        let n = 801;
        let h = 0.02;
        let eta: Vec<_> = (0..n).map(|i| i as f64 * h).collect();
        let edge = crate::edge::tests::flat_edge();
        let kernel = LaminarKernel {
            analytic: true,
            n,
            pr: 0.72,
            viscosity: crate::Viscosity {
                kind: 0,
                parameters: vec![1.7894e-5, 288.15, 110.4],
            },
            h: vec![h; n - 1],
            width: vec![h; n - 2],
            weights: vec![[-0.5 / h, 0.0, 0.5 / h]; n - 2],
        };
        // Fixture controls only; production receives every control from Python.
        let c = Controls {
            first_step: 1e-6,
            minimum_step: 1e-10,
            growth: 1.25,
            max_ratio: 2.0,
            merge: 1.3,
            relative_step: 0.03,
            velocity_change: 0.000625,
            shear_change: 0.0125,
            heat_change: 0.02,
            heat_floor: 0.1,
            separation_fraction: 0.25,
            edge_tolerance: 1e-4,
            probe: 5,
            iterations: 20,
            tolerance: 1e-10,
        };
        let result = run(
            &kernel,
            &edge,
            &[0.01, 0.02],
            &eta,
            &[-1.5 / h, 2.0 / h, -0.5 / h],
            &c,
        )
        .unwrap();
        assert_eq!(result.status, "ok");
        assert_eq!(result.positions, [0.01, 0.02]);
        for slopes in result.slopes.chunks(2) {
            assert!((slopes[0] / 2_f64.sqrt() - 0.332057).abs() < 1e-4);
        }
    }
}
