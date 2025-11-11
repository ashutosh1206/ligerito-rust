use crate::utils::evaluate_lagrange_basis;
use crate::{binary_field::BinaryField, evaluate_scaled_basis_inplace};
use rayon::prelude::*;
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
    opened_rows: &[Vec<F>],
    v_challenges: &[F],
    sorted_queries: &[usize],
    alpha: F,
) -> (Vec<F>, F)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let gr = evaluate_lagrange_basis(v_challenges);

    assert!(opened_rows.iter().all(|row| row.len() == gr.len()));
    assert_eq!(opened_rows.len(), sorted_queries.len());

    let n_rows = opened_rows.len();
    let alpha_pows = precompute_alpha_powers(alpha, n_rows);
    let mut enforced_sum = F::zero();
    let mut basis_poly = vec![F::zero(); 1 << n];

    let mut local_basis = vec![F::zero(); 1 << n];
    let mut local_sks_x = vec![F::zero(); sks_vks.len()];

    for i in 0..n_rows {
        let row = &opened_rows[i];
        // FIXME: this will panic if sorted_queries[i] > F::ValueType::MAX, handle this condition
        // and return an error in that case
        let query = F::ValueType::try_from(sorted_queries[i]).unwrap();

        let dot = row
            .iter()
            .zip(gr.iter())
            .map(|(&r, &g)| r * g)
            .fold(F::zero(), |acc, x| acc + x);

        let alpha_pow = alpha_pows[i];
        enforced_sum = enforced_sum + (dot * alpha_pow);

        let qf = F::new(query);
        evaluate_scaled_basis_inplace(&mut local_sks_x, &mut local_basis, sks_vks, qf, alpha_pow);

        for j in 0..basis_poly.len() {
            basis_poly[j] = basis_poly[j] + local_basis[j];
        }
    }

    (basis_poly, enforced_sum)
}

pub fn induce_sumcheck_poly_parallel<F>(
    n: usize,
    sks_vks: &[F],
    opened_rows: &[Vec<F>],
    v_challenges: &[F],
    sorted_queries: &[usize],
    alpha: F,
) -> (Vec<F>, F)
where
    F: Copy
        + BinaryField
        + Add<Output = F>
        + Mul<Output = F>
        + PartialEq
        + std::fmt::Debug
        + Copy
        + Sync
        + Send,
    F::ValueType: TryFrom<usize> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let gr = evaluate_lagrange_basis(v_challenges);

    assert!(opened_rows.iter().all(|row| row.len() == gr.len()));
    assert_eq!(opened_rows.len(), sorted_queries.len());

    let n_rows = opened_rows.len();
    let alpha_pows = precompute_alpha_powers(alpha, n_rows);
    let (basis_poly, enforced_sum) = (0..n_rows)
        .into_par_iter()
        .fold(
            || {
                (
                    vec![F::zero(); 1 << n],        // local_basis_poly
                    F::zero(),                      // local_enforced_sum
                    vec![F::zero(); 1 << n],        // local_basis (scratch)
                    vec![F::zero(); sks_vks.len()], // local_sks_x (scratch)
                )
            },
            |(mut local_basis_poly, mut local_enforced_sum, mut local_basis, mut local_sks_x),
             i| {
                let row = &opened_rows[i];
                let query = F::ValueType::try_from(sorted_queries[i]).unwrap();

                let dot = row
                    .iter()
                    .zip(gr.iter())
                    .map(|(&r, &g)| r * g)
                    .fold(F::zero(), |acc, x| acc + x);

                let alpha_pow = alpha_pows[i];
                local_enforced_sum = local_enforced_sum + (dot * alpha_pow);

                let qf = F::new(query);
                evaluate_scaled_basis_inplace(
                    &mut local_sks_x,
                    &mut local_basis,
                    sks_vks,
                    qf,
                    alpha_pow,
                );

                for j in 0..local_basis_poly.len() {
                    local_basis_poly[j] = local_basis_poly[j] + local_basis[j];
                }

                (
                    local_basis_poly,
                    local_enforced_sum,
                    local_basis,
                    local_sks_x,
                )
            },
        )
        .map(|(basis_poly, enforced_sum, _, _)| (basis_poly, enforced_sum))
        .reduce(
            || (vec![F::zero(); 1 << n], F::zero()),
            |(mut basis_poly1, enforced_sum1), (basis_poly2, enforced_sum2)| {
                for j in 0..basis_poly1.len() {
                    basis_poly1[j] = basis_poly1[j] + basis_poly2[j];
                }
                (basis_poly1, enforced_sum1 + enforced_sum2)
            },
        );

    (basis_poly, enforced_sum)
}
