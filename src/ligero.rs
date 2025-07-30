use crate::binary_fft::{eval_sk_at_vks, evaluate_scaled_basis_inplace};
use crate::binary_field::BinaryField;
use crate::data_structures::RecursiveLigeroWitness;
use crate::merkle_tree::build_merkle_tree;
use crate::reed_solomon::ReedSolomonEncoding;
use crate::utils::evaluate_lagrange_basis;
use std::ops::{Add, Mul};

fn poly2flatmat<F>(poly: &[F], m: usize, n: usize, inv_rate: usize) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let m_target = m * inv_rate;
    let mut flat_mat = vec![F::zero(); m_target * n];

    for j in 0..n {
        for i in 0..m {
            flat_mat[j * m_target + i] = poly[j * m + i];
        }
    }

    flat_mat
}

fn encode_cols<F>(
    poly_flat_mat: &mut Vec<F>,
    m_target: usize,
    n: usize,
    rs: &ReedSolomonEncoding<F>,
) where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    for j in 0..n {
        rs.encode_non_systematic(&mut poly_flat_mat[j * m_target..(j + 1) * m_target]);
    }
}

fn extract_row<F: BinaryField + Copy>(
    flat_mat: &[F],
    row_idx: usize,
    m_target: usize,
    n: usize,
) -> Vec<F> {
    let mut row = Vec::with_capacity(n);
    for j in 0..n {
        row.push(flat_mat[j * m_target + row_idx]);
    }
    row
}

fn extract_leaves<F: BinaryField + Copy>(
    flat_mat: &[F],
    m_target: usize,
    n: usize,
) -> Vec<Vec<u8>> {
    (0..m_target)
        .map(|row_idx| {
            let row = extract_row(flat_mat, row_idx, m_target, n);
            // Serialize the row to bytes
            row.iter().flat_map(|elem| elem.to_bytes()).collect()
        })
        .collect()
}

pub fn ligero_commit<F>(
    poly: &[F],
    m: usize,
    n: usize,
    rs: &ReedSolomonEncoding<F>,
) -> RecursiveLigeroWitness<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let inv_rate = 4;
    let mut poly_flat_mat = poly2flatmat(poly, m, n, inv_rate);
    let m_target = m * inv_rate;
    encode_cols(&mut poly_flat_mat, m_target, n, rs);

    let leaves = extract_leaves(&poly_flat_mat, m_target, n);
    let tree = build_merkle_tree(leaves);

    RecursiveLigeroWitness {
        flat_mat: poly_flat_mat,
        tree,
        num_rows: m_target,
        num_cols: n,
    }
}

pub fn ligero_verify<F>(queries: &[F::ValueType], opened_rows: &[&[F]], yr: &[F], challenges: &[F])
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let gr = evaluate_lagrange_basis(challenges);
    let n = yr.len().ilog2() as usize;
    // Why do we need to specify <F> here? Could we do something in eval_sk_at_vks to fix this?
    let sks_vks = eval_sk_at_vks::<F>(n);

    let mut local_basis = vec![F::zero(); 1 << n];
    let mut local_sks_x = vec![F::zero(); sks_vks.len()];

    for i in 0..opened_rows.len() {
        let row = opened_rows[i];
        let query: F::ValueType = queries[i];

        let dot = row
            .iter()
            .zip(gr.iter())
            .map(|(&r, &g)| r * g)
            .fold(F::zero(), |acc, x| acc + x);

        let qf = F::new(query);
        evaluate_scaled_basis_inplace(&mut local_sks_x, &mut local_basis, &sks_vks, qf, F::one());
        let e = yr
            .iter()
            .zip(local_basis.iter())
            .map(|(&y, &l)| y * l)
            .fold(F::zero(), |acc, x| acc + x);

        assert_eq!(e, dot, "Verification failed at index {}", i);
    }
}
