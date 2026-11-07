//! Allocation-free block Thomas elimination with pivoted 3x3 factors.
use crate::BandMatrix;

type Matrix = [[f64; 3]; 3];
// Three columns of the normalized upper block, followed by the reduced rhs.
type Augmented = [[f64; 4]; 3];

#[derive(Clone, Copy, Default)]
struct Lu {
    data: Matrix,
    pivots: [usize; 3],
    reciprocal: [f64; 3],
}

impl Lu {
    #[inline]
    fn factor(&mut self, matrix: Matrix) -> bool {
        // The fixed 3x3 factorization keeps the three rows in locals. Pivot
        // comparisons and elimination order match ordinary partial-pivot LU.
        let [mut a, mut b, mut c] = matrix;
        let mut first = 0;
        let mut largest = a[0].abs();
        if b[0].abs() > largest {
            first = 1;
            largest = b[0].abs();
        }
        if c[0].abs() > largest {
            first = 2;
        }
        match first {
            1 => std::mem::swap(&mut a, &mut b),
            2 => std::mem::swap(&mut a, &mut c),
            _ => {}
        }
        if a[0] == 0.0 {
            return false;
        }
        let r0 = a[0].recip();
        b[0] *= r0;
        c[0] *= r0;
        b[1] -= b[0] * a[1];
        b[2] -= b[0] * a[2];
        c[1] -= c[0] * a[1];
        c[2] -= c[0] * a[2];
        let second = if c[1].abs() > b[1].abs() {
            std::mem::swap(&mut b, &mut c);
            2
        } else {
            1
        };
        if b[1] == 0.0 {
            return false;
        }
        let r1 = b[1].recip();
        c[1] *= r1;
        c[2] -= c[1] * b[2];
        if c[2] == 0.0 {
            return false;
        }
        let r2 = c[2].recip();
        self.data = [a, b, c];
        self.pivots = [first, second, 2];
        self.reciprocal = [r0, r1, r2];
        // Zero reciprocals identify infinite pivots; overflowed reciprocals
        // also need fallback even if the original nonzero pivot was finite.
        self.reciprocal.iter().all(|v| v.is_finite() && *v != 0.0)
    }

    #[inline]
    fn solve_four(&self, mut rhs: Augmented) -> Augmented {
        rhs.swap(0, self.pivots[0]);
        rhs.swap(1, self.pivots[1]);
        let [mut a, mut b, mut c] = rhs;
        // Each column has the same factors; all four columns are independent.
        // Fixed rows remove dynamic triangular indexing from this hot path.
        for k in 0..4 {
            b[k] -= self.data[1][0] * a[k];
            c[k] -= self.data[2][0] * a[k];
            c[k] -= self.data[2][1] * b[k];
            c[k] *= self.reciprocal[2];
            b[k] -= self.data[1][2] * c[k];
            b[k] *= self.reciprocal[1];
            a[k] -= self.data[0][1] * b[k];
            a[k] -= self.data[0][2] * c[k];
            a[k] *= self.reciprocal[0];
        }
        [a, b, c]
    }
}

pub(crate) struct Blocks {
    pub lower: Vec<Matrix>,
    pub diagonal: Vec<Matrix>,
    pub upper: Vec<Matrix>,
    factors: Vec<Lu>,
    reduced: Vec<Augmented>,
    work: Vec<[f64; 3]>,
}

impl Blocks {
    pub fn new(n: usize) -> Self {
        Self {
            lower: vec![[[0.0; 3]; 3]; n],
            diagonal: vec![[[0.0; 3]; 3]; n],
            upper: vec![[[0.0; 3]; 3]; n],
            factors: vec![Lu::default(); n],
            reduced: vec![[[0.0; 4]; 3]; n],
            work: vec![[0.0; 3]; n],
        }
    }

    // Normalize each reduced row with its LU. This is algebraically the same
    // Schur complement as L_i = A_i B'_(i-1)^-1, without explicit inverses.
    // Preserve both the input blocks and rhs until success for banded fallback.
    pub fn solve_block(&mut self, rhs: &mut [f64]) -> bool {
        let mut previous = [[0.0; 4]; 3];
        let blocks = self.lower.iter().zip(&self.diagonal).zip(&self.upper);
        let storage = self.factors.iter_mut().zip(&mut self.reduced);
        for (i, ((((lower, diagonal), upper), (factor, reduced)), rhs_row)) in
            blocks.zip(storage).zip(rhs.as_chunks::<3>().0).enumerate()
        {
            // Fixed-size locals keep the sequential recurrence off Vec indexing.
            let lower = *lower;
            let mut diagonal = *diagonal;
            let mut combined = [
                [upper[0][0], upper[0][1], upper[0][2], rhs_row[0]],
                [upper[1][0], upper[1][1], upper[1][2], rhs_row[1]],
                [upper[2][0], upper[2][1], upper[2][2], rhs_row[2]],
            ];
            if i > 0 {
                for r in 0..3 {
                    for (k, previous_row) in previous.iter().enumerate() {
                        combined[r][3] -= lower[r][k] * previous_row[3];
                        for (c, value) in diagonal[r].iter_mut().enumerate() {
                            *value -= lower[r][k] * previous_row[c];
                        }
                    }
                }
            }
            if !factor.factor(diagonal) {
                return false;
            }
            previous = factor.solve_four(combined);
            *reduced = previous;
        }
        let mut next = [0.0; 3];
        for (i, (reduced, solution)) in self.reduced.iter().zip(&mut self.work).rev().enumerate() {
            let mut current = [reduced[0][3], reduced[1][3], reduced[2][3]];
            // The last block has no upper neighbor, even if unused storage is
            // nonzero. All earlier normalized upper entries take part here.
            if i > 0 {
                for r in 0..3 {
                    for (k, value) in next.iter().enumerate() {
                        current[r] -= reduced[r][k] * value;
                    }
                }
            }
            *solution = current;
            next = current;
        }
        // Non-finite elimination/solve arithmetic propagates to a pivot or to
        // this solution. Keep the sole solution scan before publishing rhs.
        if self.work.iter().flatten().any(|v| !v.is_finite()) {
            return false;
        }
        for (row, values) in rhs.as_chunks_mut::<3>().0.iter_mut().zip(&self.work) {
            row.copy_from_slice(values);
        }
        true
    }

    pub fn copy_to_band(&self, band: &mut BandMatrix) {
        band.data.fill(0.0);
        for i in 0..self.diagonal.len() {
            for r in 0..3 {
                for c in 0..3 {
                    band.set(3 * i + r, 3 * i + c, self.diagonal[i][r][c]);
                    if i > 0 {
                        band.set(3 * i + r, 3 * (i - 1) + c, self.lower[i][r][c]);
                    }
                    if i + 1 < self.diagonal.len() {
                        band.set(3 * i + r, 3 * (i + 1) + c, self.upper[i][r][c]);
                    }
                }
            }
        }
    }

    pub fn solve(&mut self, rhs: &mut [f64], fallback: &mut BandMatrix) -> bool {
        if self.solve_block(rhs) {
            return true;
        }
        self.copy_to_band(fallback);
        fallback.solve(rhs) && rhs.iter().all(|v| v.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn four_right_hand_sides_with_two_row_pivots() {
        let matrix = [[0.0, 2.0, 1.0], [1.0, 0.0, 4.0], [3.0, 1.0, 1.0]];
        let expected = [
            [1.0, -2.0, 3.0, 4.0],
            [5.0, 6.0, -7.0, 8.0],
            [9.0, 10.0, 11.0, -12.0],
        ];
        let mut rhs = [[0.0; 4]; 3];
        for (row, values) in rhs.iter_mut().zip(matrix) {
            for (value, solution) in values.iter().zip(expected) {
                for (entry, x) in row.iter_mut().zip(solution) {
                    *entry += value * x;
                }
            }
        }
        let mut lu = Lu::default();
        assert!(lu.factor(matrix));
        assert_eq!(lu.pivots, [2, 2, 2]);
        for (actual, expected) in lu
            .solve_four(rhs)
            .iter()
            .flatten()
            .zip(expected.iter().flatten())
        {
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn reciprocal_overflow_falls_back_to_banded_division() {
        let mut blocks = Blocks::new(1);
        for k in 0..3 {
            blocks.diagonal[0][k][k] = 1e-320;
        }
        let original = [1e-320; 3];
        let mut rhs = original;
        assert!(!blocks.solve_block(&mut rhs));
        assert_eq!(rhs, original);
        assert!(blocks.solve(&mut rhs, &mut BandMatrix::new(3)));
        assert_eq!(rhs, [1.0; 3]);
    }

    #[test]
    fn failed_recurrence_and_back_substitution_preserve_rhs() {
        // Fail at a later diagonal block, after the first row was reduced.
        let mut blocks = Blocks::new(2);
        for k in 0..3 {
            blocks.diagonal[0][k][k] = 1.0;
        }
        let original = [1.0; 6];
        let mut rhs = original;
        assert!(!blocks.solve_block(&mut rhs));
        assert_eq!(rhs, original);
        // Finite inputs and finite pivots overflow in the four-column solve.
        let mut blocks = Blocks::new(1);
        blocks.diagonal[0] = [[1.0, 1e308, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let original = [1e308; 3];
        let mut rhs = original;
        assert!(!blocks.solve_block(&mut rhs));
        assert_eq!(rhs, original);
        // Overflow only in the final block back substitution (no lower block).
        let mut blocks = Blocks::new(2);
        for i in 0..2 {
            for k in 0..3 {
                blocks.diagonal[i][k][k] = 1.0;
            }
        }
        blocks.upper[0][0][0] = 1e308;
        let original = [1.0, 1.0, 1.0, 1e308, 1.0, 1.0];
        let mut rhs = original;
        assert!(!blocks.solve_block(&mut rhs));
        assert_eq!(rhs, original);
    }

    #[test]
    fn random_blocks_and_pivoting() {
        let mut seed = 987654_u64;
        for diagonal_shift in [5.0, 0.0] {
            for _ in 0..20 {
                let n = 19;
                let mut blocks = Blocks::new(n);
                for i in 0..n {
                    for r in 0..3 {
                        for c in 0..3 {
                            for matrix in [
                                &mut blocks.lower[i],
                                &mut blocks.diagonal[i],
                                &mut blocks.upper[i],
                            ] {
                                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                                matrix[r][c] = (seed >> 32) as f64 / u32::MAX as f64 - 0.5;
                            }
                        }
                        blocks.diagonal[i][r][r] += diagonal_shift;
                    }
                }
                let expected: Vec<_> = (0..3 * n).map(|i| (i as f64).sin()).collect();
                let mut rhs = vec![0.0; 3 * n];
                for i in 0..n {
                    for r in 0..3 {
                        for c in 0..3 {
                            rhs[3 * i + r] += blocks.diagonal[i][r][c] * expected[3 * i + c];
                            if i > 0 {
                                rhs[3 * i + r] += blocks.lower[i][r][c] * expected[3 * (i - 1) + c];
                            }
                            if i + 1 < n {
                                rhs[3 * i + r] += blocks.upper[i][r][c] * expected[3 * (i + 1) + c];
                            }
                        }
                    }
                }
                let mut band = BandMatrix::new(3 * n);
                blocks.copy_to_band(&mut band);
                let mut reference = rhs.clone();
                assert!(band.solve(&mut reference));
                assert!(blocks.solve_block(&mut rhs));
                for ((actual, reference), expected) in rhs.iter().zip(reference).zip(expected) {
                    assert!((actual - reference).abs() < 1e-9);
                    assert!((actual - expected).abs() < 1e-9);
                }
            }
        }
    }

    #[test]
    fn singular_first_block_falls_back_without_corrupting_rhs() {
        let mut blocks = Blocks::new(2);
        for k in 0..3 {
            blocks.upper[0][k][k] = 1.0;
            blocks.lower[1][k][k] = 1.0;
            blocks.diagonal[1][k][k] = 2.0;
        }
        let original = [4.0, 5.0, 6.0, 9.0, 12.0, 15.0];
        let mut rhs = original;
        assert!(!blocks.solve_block(&mut rhs));
        assert_eq!(rhs, original);
        assert!(blocks.solve(&mut rhs, &mut BandMatrix::new(6)));
        assert_eq!(rhs, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        blocks.lower[1] = [[0.0; 3]; 3];
        assert!(!blocks.solve(&mut rhs, &mut BandMatrix::new(6)));
    }

    #[test]
    fn nonfinite_intermediate_and_result_fail() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut blocks = Blocks::new(2);
            for i in 0..2 {
                for k in 0..3 {
                    blocks.diagonal[i][k][k] = 1.0;
                }
            }
            blocks.diagonal[0][0][0] = bad;
            assert!(!blocks.solve_block(&mut [1.0; 6]));
        }
        let mut blocks = Blocks::new(2);
        for i in 0..2 {
            for k in 0..3 {
                blocks.diagonal[i][k][k] = 1e-300;
            }
        }
        assert!(!blocks.solve_block(&mut [1e300; 6]));
    }
}
