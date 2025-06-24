use crate::binary_field::{BinaryElem16, BinaryElem32, BinaryField};
use std::ops::{Add, Mul};

pub fn fft<F>(v: &mut [F], twiddles: &[F])
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(is_power_of_2(v.len()));

    fft_twiddles(v, twiddles, Some(1));
}

pub fn ifft<F>(v: &mut [F], twiddles: &[F])
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(is_power_of_2(v.len()));

    ifft_twiddles(v, twiddles, Some(1));
}

pub fn compute_twiddles<F>(log_n: usize, beta: Option<F>) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let beta = beta.unwrap_or_else(|| F::zero());

    let mut twiddles = vec![F::zero(); (1 << log_n) - 1];
    let mut layer_size = 1 << (log_n - 1);
    let mut layer = vec![F::zero(); layer_size];
    let mut write_at = layer_size;
    let mut s_prev_at_root = layer_0(&mut layer, beta, log_n);
    twiddles[(write_at - 1)..].copy_from_slice(&layer);

    for _ in 0..(log_n - 1) {
        write_at >>= 1;
        layer_size = write_at;

        s_prev_at_root = layer_i(&mut layer, layer_size, s_prev_at_root);
        let s_inv = s_prev_at_root.inverse().unwrap();

        for i in 0..layer_size {
            twiddles[write_at + i - 1] = s_inv * layer[i];
        }
    }

    twiddles
}

fn layer_0<F>(layer: &mut Vec<F>, beta: F, k: usize) -> F
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    for i in 0usize..(1 << (k - 1)) {
        let mut l0i = beta;
        l0i = l0i + F::new(F::ValueType::try_from(i << 1).unwrap());
        layer[i] = l0i;
    }
    F::new(F::ValueType::try_from(1).unwrap())
}

fn layer_i<F>(layer: &mut Vec<F>, layer_size: usize, s_prev_at_root: F) -> F
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let prev_layer_size = layer_size * 2;
    let s_at_root = compute_s_at_root(layer, s_prev_at_root);

    for (write_idx, read_idx) in (0..prev_layer_size).step_by(2).enumerate() {
        let s_prev = layer[read_idx];
        layer[write_idx] = next_s(s_prev, s_prev_at_root);
    }

    s_at_root
}

fn compute_s_at_root<F>(prev_layer: &Vec<F>, s_prev_at_root: F) -> F
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    next_s(prev_layer[1] + prev_layer[0], s_prev_at_root)
}

fn next_s<F>(s_prev: F, s_prev_at_root: F) -> F
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    s_prev * s_prev + s_prev_at_root * s_prev
}

// Internal implementation (private)
fn fft_twiddles<F>(v: &mut [F], twiddles: &[F], idx: Option<usize>)
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    if v.len() == 1 {
        return;
    }
    let idx = idx.unwrap_or_else(|| 1);

    fft_mul(v, twiddles[idx - 1]);
    let (u, w) = split_half(v);

    fft_twiddles(u, twiddles, Some(idx * 2));
    fft_twiddles(w, twiddles, Some(idx * 2 + 1));
}

fn ifft_twiddles<F>(v: &mut [F], twiddles: &[F], idx: Option<usize>)
where
    F: Copy + BinaryField + Add<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    if v.len() == 1 {
        return;
    }
    let idx = idx.unwrap_or_else(|| 1);

    let (u, w) = split_half(v);

    ifft_twiddles(u, twiddles, Some(idx * 2));
    ifft_twiddles(w, twiddles, Some(idx * 2 + 1));

    ifft_mul(v, twiddles[idx - 1]);
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

fn is_power_of_2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
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
    use crate::binary_field::BinaryElem16;

    #[test]
    fn test_compute_twiddles_matches_julia() {
        let twiddles = compute_twiddles::<BinaryElem16>(10, None);

        // Expected values from Julia: compute_twiddles(BinaryElem16, 10)
        let expected_start = [0x0000, 0x0000, 0x6ba9, 0x0000, 0x6a23];

        println!("Rust twiddles length: {}", twiddles.len());
        println!("First 5 twiddles:");
        for (i, &t) in twiddles.iter().take(5).enumerate() {
            println!("  {}: 0x{:04x}", i + 1, t.value);
        }

        // Check length matches Julia (2^10 - 1 = 1023)
        assert_eq!(twiddles.len(), 1023);

        // Check first few values match
        for (i, &expected) in expected_start.iter().enumerate() {
            assert_eq!(twiddles[i].value, expected, "Mismatch at index {}", i);
        }
    }
}
