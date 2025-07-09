use rand::distr::{Distribution, StandardUniform};

use crate::binary_field::{BinaryField, random};
use crate::multilinear_poly::{MultiLinearPoly, eval_013_product};
use std::ops::{Add, Mul};

pub struct QuadraticEvals<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    e0: F,
    e1: F,
    e2: F,
}

impl<F> QuadraticEvals<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(e0: F, e1: F, e2: F) -> Self {
        Self { e0, e1, e2 }
    }
}

pub struct QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    a: F,
    b: F,
    c: F,
}

impl<F> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    // TODO: do we really need this function if we are not performing
    // any operation on a,b,c and simply initializing the struct?
    pub fn new(a: F, b: F, c: F) -> Self {
        Self { a, b, c }
    }

    pub fn eval_quadratic(&self, r: F) -> F {
        self.a * r * r + self.b * r + self.c
    }
}

pub fn quadratic_from_evals<F>(at0: F, at1: F, atx: F, x: Option<F>) -> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    // TODO: Can we unwrap x in a simpler way from Option<F>?
    let x = x.unwrap_or_else(|| F::new(F::ValueType::try_from(3).unwrap()));
    let numerator = atx + at0 + x * (at1 + at0);
    let denominator = x * x + x;
    // TODO: is there a better way to do this than .unwrap()?
    let a = numerator * denominator.inverse().unwrap();
    let b = at1 + at0 + a;
    return QuadraticPoly { a, b, c: at0 };
}

pub fn fold_quadratic<F>(p1: QuadraticPoly<F>, p2: QuadraticPoly<F>, alpha: F) -> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    QuadraticPoly::new(
        p1.a + alpha * p2.a,
        p1.b + alpha * p2.b,
        p1.c + alpha * p2.c,
    )
}

pub fn sumcheck_prover<F>(f: &MultiLinearPoly<F>, claimed_sum: F) -> (Vec<(F, F)>, Vec<F>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
    StandardUniform: Distribution<F::ValueType>,
{
    let num_vars = f.num_vars();
    let mut transcript: Vec<(F, F)> = Vec::with_capacity(num_vars);
    let mut ris: Vec<F> = Vec::with_capacity(num_vars);
    let mut current_poly = f.clone();
    let mut current_sum = claimed_sum;

    for _ in 0..num_vars {
        let (s0, s1) = current_poly.eval_012();
        assert_eq!(s0 + s1, current_sum);

        transcript.push((s0, s1));
        let r_i = random::<F>();
        ris.push(r_i);

        current_poly = current_poly.partial_eval(vec![r_i]);
        current_sum = current_poly.sum();
    }

    (transcript, ris)
}

pub fn sumcheck_verifier<F>(
    transcript: Vec<(F, F)>,
    ris: Vec<F>,
    claimed_sum: F,
    f: &MultiLinearPoly<F>,
) -> bool
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut h = claimed_sum;
    for (i, (g0, g1)) in transcript.iter().enumerate() {
        assert_eq!(*g0 + *g1, h);
        h = *g0 * (F::one() + ris[i]) + (*g1 * ris[i]);
    }

    let f_eval = f.partial_eval(ris).sum();
    f_eval == h
}

pub fn double_sumcheck_prover<F>(
    fp: &MultiLinearPoly<F>,
    gp: &MultiLinearPoly<F>,
) -> (Vec<(F, F, F)>, Vec<F>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
    StandardUniform: Distribution<F::ValueType>,
{
    let num_vars = fp.num_vars();
    let mut transcript: Vec<(F, F, F)> = Vec::with_capacity(num_vars);
    let mut f = fp.clone();
    let mut g = gp.clone();
    let mut ris: Vec<F> = Vec::with_capacity(num_vars);

    for _ in 0..num_vars {
        let (s0, s1, s2) = eval_013_product(&f, &g);
        transcript.push((s0, s1, s2));

        let r_i = random::<F>();
        ris.push(r_i);

        f = f.partial_eval(vec![r_i]);
        g = g.partial_eval(vec![r_i]);
    }

    (transcript, ris)
}

pub fn double_sumcheck_verifier<F>(
    transcript: Vec<(F, F, F)>,
    ris: Vec<F>,
    claimed_sum: F,
    fp: &MultiLinearPoly<F>,
    gp: &MultiLinearPoly<F>,
) -> bool
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut h = claimed_sum;
    for (i, (g0, g1, g2)) in transcript.iter().enumerate() {
        assert_eq!(*g0 + *g1, h);
        let gi = quadratic_from_evals(*g0, *g1, *g2, None);
        h = gi.eval_quadratic(ris[i]);
    }

    // TODO: there has to be a better way than cloning ris twice
    let f_eval = fp.partial_eval(ris.clone()).sum();
    let g_eval = gp.partial_eval(ris.clone()).sum();
    f_eval * g_eval == h
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_field::BinaryElem16;

    // Helper function to create a random polynomial for testing
    fn random_poly(k: usize) -> MultiLinearPoly<BinaryElem16> {
        let evals: Vec<BinaryElem16> = (0..(1 << k)).map(|_| random::<BinaryElem16>()).collect();
        MultiLinearPoly::new(evals)
    }

    #[test]
    fn test_sumcheck_basic() {
        // Test with small polynomial (2 variables)
        let f = random_poly(2);
        let claimed_sum = f.sum();

        // Run prover
        let (transcript, ris) = sumcheck_prover(&f, claimed_sum);

        // Run verifier
        let verified = sumcheck_verifier(transcript, ris, claimed_sum, &f);
        assert!(verified, "Sumcheck verification should pass");
    }

    #[test]
    fn test_double_sumcheck_basic() {
        // Test with small polynomials (2 variables each)
        let f = random_poly(2);
        let g = random_poly(2);

        // Calculate claimed sum: sum of element-wise products
        let claimed_sum: BinaryElem16 = f
            .evals()
            .iter()
            .zip(g.evals().iter())
            .map(|(fi, gi)| *fi * *gi)
            .fold(BinaryElem16::zero(), |acc, x| acc + x);

        // Run prover
        let (transcript, ris) = double_sumcheck_prover(&f, &g);

        // Run verifier
        let verified = double_sumcheck_verifier(transcript, ris, claimed_sum, &f, &g);
        assert!(verified, "Double sumcheck verification should pass");
    }

    #[test]
    fn test_sumcheck_larger() {
        // Test with larger polynomial (4 variables)
        let f = random_poly(4);
        let claimed_sum = f.sum();

        let (transcript, ris) = sumcheck_prover(&f, claimed_sum);
        let verified = sumcheck_verifier(transcript, ris, claimed_sum, &f);
        assert!(verified, "Larger sumcheck verification should pass");
    }

    #[test]
    fn test_sumcheck_single_variable() {
        // Test edge case with 1 variable polynomial
        let f = random_poly(1);
        let claimed_sum = f.sum();

        let (transcript, ris) = sumcheck_prover(&f, claimed_sum);
        let verified = sumcheck_verifier(transcript, ris, claimed_sum, &f);
        assert!(verified, "Single variable sumcheck should pass");
    }

    #[test]
    fn test_double_sumcheck_larger() {
        // Test with larger polynomials (3 variables each)
        let f = random_poly(3);
        let g = random_poly(3);

        // Calculate claimed sum: sum of element-wise products
        let claimed_sum: BinaryElem16 = f
            .evals()
            .iter()
            .zip(g.evals().iter())
            .map(|(fi, gi)| *fi * *gi)
            .fold(BinaryElem16::zero(), |acc, x| acc + x);

        let (transcript, ris) = double_sumcheck_prover(&f, &g);
        let verified = double_sumcheck_verifier(transcript, ris, claimed_sum, &f, &g);
        assert!(verified, "Larger double sumcheck should pass");
    }

    #[test]
    fn test_sumcheck_known_values() {
        // Test with known small values for deterministic behavior
        let evals = vec![
            BinaryElem16::new(1),
            BinaryElem16::new(2),
            BinaryElem16::new(3),
            BinaryElem16::new(4),
        ];
        let f = MultiLinearPoly::new(evals);
        let claimed_sum = f.sum(); // Should be 1+2+3+4 = 10 in binary field

        let (transcript, ris) = sumcheck_prover(&f, claimed_sum);
        let verified = sumcheck_verifier(transcript, ris, claimed_sum, &f);
        assert!(verified, "Known values sumcheck should pass");
    }
}
