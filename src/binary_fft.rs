use crate::binary_field::{BinaryElem16, BinaryElem32, BinaryField};
use std::{fmt::write, ops::Add};

pub fn fft<F>(v: &mut [F], twiddles: &[F]) {
    // TODO: Implement
}

pub fn ifft<F>(v: &mut [F], twiddles: &[F]) {
    // TODO: Implement
}

pub fn compute_twiddles<F>(log_n: usize, beta: Option<F>) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: From<usize>,
{
    let beta = beta.unwrap_or_else(|| F::zero());

    let mut twiddles = vec![F::zero(); (1 << log_n) - 1];
    let mut layer_size = 1 << (log_n - 1);
    let mut layer = vec![F::zero(); layer_size];
    let mut write_at = layer_size;
    let s_prev_at_root = layer_0(&mut layer, beta, log_n);
    twiddles[(write_at - 1)..].copy_from_slice(&layer);

    for _ in 0..(log_n - 1) {
        write_at >>= 1;
        layer_size = write_at;
    }

    // TODO: Implement the rest
    twiddles
}

fn layer_0<F>(layer: &mut Vec<F>, beta: F, k: usize) -> F
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: From<usize>,
{
    for i in 0usize..(1 << (k - 1)) {
        let mut l0i = beta;
        l0i = l0i + F::new(F::ValueType::from(i << 1));
        layer[i] = l0i;
    }
    F::new(F::ValueType::from(1))
}

fn layer_i<F>(layer: &mut Vec<F>, layer_size: usize, s_prev_at_root: F) {}

// Internal implementation (private)
fn fft_twiddles<F>(v: &mut [F], twiddles: &[F], idx: usize) {
    // TODO: Implement
}

fn ifft_twiddles<F>(v: &mut [F], twiddles: &[F], idx: usize) {
    // TODO: Implement
}

fn fft_mul<F>(v: &mut [F], lambda: F) {
    // TODO: Implement
}

fn ifft_mul<F>(v: &mut [F], lambda: F) {
    // TODO: Implement
}

fn split_half<F>(v: &mut [F]) -> (&mut [F], &mut [F]) {
    let mid = v.len() / 2;
    v.split_at_mut(mid)
}

// Utilities (public)
pub fn eval_sk_at_vks<F>(n: usize) -> Vec<F> {
    // TODO: Implement
    Vec::new()
}

pub fn evaluate_basis<F>(basis_len: usize, sks_vks: &[F], x: F) -> Vec<F> {
    // TODO: Implement
    Vec::new()
}

pub fn compute_pis<F>(pis_len: usize, sks_vks: &[F]) -> Vec<F> {
    // TODO: Implement
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_fft() {}
}
