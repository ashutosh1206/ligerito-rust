use crate::binary_fft::is_power_of_2;
use crate::binary_field::BinaryField;
use rayon::prelude::*;
use std::ops::{Add, Mul};

#[derive(Debug, Clone)]
pub struct MultiLinearPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    evals: Vec<F>,
    n: usize,
}

impl<F> MultiLinearPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + Send + Sync,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(evals: Vec<F>) -> Self {
        let n: usize = evals.len().ilog2() as usize;
        assert!(is_power_of_2(evals.len()));

        Self { evals, n }
    }

    pub fn sum(&self) -> F {
        self.evals.iter().fold(F::zero(), |acc, &x| acc + x)
    }

    pub fn num_vars(&self) -> usize {
        self.n
    }

    pub fn evals(&self) -> &Vec<F> {
        &self.evals
    }

    pub fn scale_evals(&mut self, alpha: F) {
        for eval in &mut self.evals {
            *eval = *eval * alpha;
        }
    }

    pub fn partial_eval_at_0(&self) -> Self {
        let half = self.evals.len() >> 1;
        MultiLinearPoly::new(self.evals[..half].to_vec())
    }

    pub fn partial_eval_at_1(&self) -> Self {
        let half = self.evals.len() >> 1;
        MultiLinearPoly::new(self.evals[half..].to_vec())
    }

    pub fn eval_012(&self) -> (F, F) {
        let f0 = self.partial_eval_at_0();
        let f1 = self.partial_eval_at_1();
        (f0.sum(), f1.sum())
    }

    pub fn partial_eval(&self, rs: &[F]) -> Self {
        let mut n = self.evals.len() / 2;
        let mut partial_evals: Vec<F> = Vec::with_capacity(n);
        let one = F::one();

        for i in 0..n {
            partial_evals.push((one + rs[0]) * self.evals[i] + rs[0] * self.evals[i + n]);
        }

        for i in 1..rs.len() {
            n /= 2;

            let (left, right) = partial_evals.split_at_mut(n);
            left.par_iter_mut()
                .zip(right.par_iter())
                .for_each(|(l, r)| {
                    *l = (one + rs[i]) * *l + rs[i] * *r;
                });
        }

        MultiLinearPoly::new(partial_evals[0..n].to_vec())
    }
}

pub fn eval_013_product<F>(f: &MultiLinearPoly<F>, g: &MultiLinearPoly<F>) -> (F, F, F)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + Send + Sync,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let f0 = f.partial_eval_at_0();
    let g0 = g.partial_eval_at_0();
    let s0 = f0
        .evals
        .iter()
        .zip(g0.evals.iter())
        .map(|(fi, gi)| *fi * *gi)
        .fold(F::zero(), |acc, x| acc + x);

    let f1 = f.partial_eval_at_1();
    let g1 = g.partial_eval_at_1();
    let s1 = f1
        .evals
        .iter()
        .zip(g1.evals.iter())
        .map(|(fi, gi)| *fi * *gi)
        .fold(F::zero(), |acc, x| acc + x);

    let field_three = F::new(F::ValueType::try_from(3).unwrap());
    let f2_vals = f1
        .evals
        .iter()
        .zip(f0.evals.iter())
        .map(|(f1i, f0i)| (*f1i * field_three) + (F::one() + field_three) * *f0i);
    let g2_vals = g1
        .evals
        .iter()
        .zip(g0.evals.iter())
        .map(|(g1i, g0i)| (*g1i * field_three) + (F::one() + field_three) * *g0i);
    let s2 = f2_vals
        .zip(g2_vals)
        .map(|(f2i, g2i)| f2i * g2i)
        .fold(F::zero(), |acc, x| acc + x);

    (s0, s1, s2)
}
