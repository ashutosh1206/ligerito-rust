use crate::binary_field::BinaryField;
use crate::utils::evaluate_lagrange_basis;
use std::ops::{Add, Mul};

pub fn precompute_alpha_powers<F>(alpha: F, n: usize) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut alpha_pows = Vec::with_capacity(n);
    alpha_pows.push(F::one());

    for i in 1..n {
        alpha_pows.push(alpha_pows[i - 1] * alpha);
    }

    alpha_pows
}

pub fn induce_sumcheck_poly<F>(
    n: usize,
    sks_vks: &[F],
    opened_rows: &[&[F]],
    v_challenges: &[F],
    sorted_queries: &[usize],
    alpha: F,
) where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let gr = evaluate_lagrange_basis(v_challenges);

    assert!(opened_rows.iter().all(|row| row.len() == gr.len()));
    assert_eq!(opened_rows.len(), sorted_queries.len());

    let n_rows = opened_rows.len();
    let alpha_pows = precompute_alpha_powers(alpha, n_rows);
    let partial_basis = vec![F::zero(); 1 << n];
}
