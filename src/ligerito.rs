use crate::binary_field::BinaryField;
use crate::multilinear_poly::{MultiLinearPoly, eval_013_product};
use crate::sumcheck::QuadraticEvals;
use std::ops::{Add, Mul};

pub struct SumcheckProverInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    f: MultiLinearPoly<F>,
    basis_polys: Vec<MultiLinearPoly<F>>,
    sum: F,
    transcript: Vec<(F, F, F)>,
    to_be_glued: Option<MultiLinearPoly<F>>,
}

impl<F> SumcheckProverInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(f: MultiLinearPoly<F>, b1: MultiLinearPoly<F>, h1: F) -> (Self, QuadraticEvals<F>) {
        let mut transcript: Vec<(F, F, F)> = vec![];
        let (s0, s1, s2) = eval_013_product(&f, &b1);
        assert_eq!(s0 + s1, h1);
        transcript.push((s0, s1, s2));
        let quad_evals = QuadraticEvals::new(s0, s1, s2);

        (
            Self {
                f,
                basis_polys: vec![b1],
                sum: h1,
                transcript,
                to_be_glued: None,
            },
            quad_evals,
        )
    }

    pub fn introduce_new(&mut self, bi: MultiLinearPoly<F>, h: F) -> QuadraticEvals<F> {
        let (s0, s1, s2) = eval_013_product(&self.f, &bi);
        assert_eq!(s0 + s1, h);
        self.transcript.push((s0, s1, s2));
        self.to_be_glued = Some(bi);

        QuadraticEvals::new(s0, s1, s2)
    }

    pub fn glue(&mut self, alpha: F) {
        assert!(self.to_be_glued.is_some());
        let mut poly = self.to_be_glued.take().unwrap();
        poly.scale_evals(alpha);
        self.basis_polys.push(poly);
        // self.to_be_glued.take() already leaves None in its place
        // So this is redundant and a relic of Julia port
        self.to_be_glued = None;
    }

    pub fn eval_01x_product(&self, x: Option<F>) -> (F, F, F) {
        let x = x.unwrap_or_else(|| F::new(F::ValueType::try_from(3).unwrap()));

        let val_at_0 = self.f.partial_eval_at_0();
        let val_at_1 = self.f.partial_eval_at_1();
        let f0 = val_at_0.evals().clone();
        let f1 = val_at_1.evals().clone();
        let f2 = evals_at_x(&f0, &f1, x);

        let bs_val_at_0 = self.basis_polys[0].partial_eval_at_0();
        let bs_val_at_1 = self.basis_polys[0].partial_eval_at_1();
        let mut b0s_sum = bs_val_at_0.evals().clone();
        let mut b1s_sum = bs_val_at_1.evals().clone();
        let mut b2s_sum = evals_at_x(&b0s_sum, &b1s_sum, x);

        for i in 1..self.basis_polys.len() {
            let bi_val_at_0 = self.basis_polys[i].partial_eval_at_0();
            let bi_val_at_1 = self.basis_polys[i].partial_eval_at_1();
            let b0i = bi_val_at_0.evals().clone();
            let b1i = bi_val_at_1.evals().clone();
            let b2i = evals_at_x(&b0i, &b1i, x);

            for i in 0..b0s_sum.len() {
                b0s_sum[i] = b0s_sum[i] + b0i[i];
            }
            for i in 0..b1s_sum.len() {
                b1s_sum[i] = b1s_sum[i] + b1i[i];
            }
            for i in 0..b2s_sum.len() {
                b2s_sum[i] = b2s_sum[i] + b2i[i];
            }
        }
        let s0 = f0
            .iter()
            .zip(b0s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        let s1 = f1
            .iter()
            .zip(b1s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        let s2 = f2
            .iter()
            .zip(b2s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        (s0, s1, s2)
    }

    pub fn fold(&mut self, r: F) -> QuadraticEvals<F> {
        self.f = self.f.partial_eval(vec![r]);
        for i in 0..self.basis_polys.len() {
            self.basis_polys[i] = self.basis_polys[i].partial_eval(vec![r]);
        }
        let (s0, s1, s2) = self.eval_01x_product(None);
        self.transcript.push((s0, s1, s2));
        QuadraticEvals::new(s0, s1, s2)
    }
}

fn evals_at_x<F>(evals_at_0: &Vec<F>, evals_at_1: &Vec<F>, alpha: F) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    evals_at_0
        .iter()
        .zip(evals_at_1.iter())
        .map(|(&val0, &val1)| alpha * val0 + (F::one() + alpha) * val1)
        .collect()
}
