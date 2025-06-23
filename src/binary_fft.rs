use crate::binary_field::{BinaryElem16, BinaryElem32, BinaryField};
use std::ops::Add;

pub fn fft<F>(v: &mut [F], twiddles: &[F]) {
    // TODO: Implement
}

pub fn ifft<F>(v: &mut [F], twiddles: &[F]) {
    // TODO: Implement
}

pub fn compute_twiddles<F>(log_n: usize, beta: Option<F>) -> Vec<F>
where
    F: Copy + BinaryField,
{
    let beta = beta.unwrap_or_else(|| F::zero());

    let mut twiddles = vec![F::zero(); (1 << log_n) - 1];
    let layer_size = 1 << (log_n - 1);
    let mut layer = vec![F::zero(); layer_size];
    let write_at = layer_size;

    // TODO: Implement the rest
    twiddles
}

fn layer_0<F>(mut layer: Vec<F>, beta: F, k: usize) -> F
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: From<usize>,
{
    for i in 1usize..=(1<<(k-1)) {
        let mut l0i = beta;
        l0i = l0i + F::new(F::ValueType::from((i-1) << 1));
        layer[i-1] = l0i;
    }
    F::new(F::ValueType::from(1))
}

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
