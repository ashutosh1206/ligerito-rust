use crate::binary_field::BinaryField;
use std::ops::{Add, Mul};

pub fn fft<F>(v: &mut [F], twiddles: &[F], parallel: Option<bool>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + Send + Sync,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(is_power_of_2(v.len()));

    if parallel.unwrap_or_else(|| false) {
        fft_twiddles(v, twiddles, Some(1));
    } else {
        fft_twiddles_parallel(v, twiddles, Some(1), None);
    }
}

pub fn ifft<F>(v: &mut [F], twiddles: &[F])
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
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

fn fft_twiddles_parallel<F>(
    v: &mut [F],
    twiddles: &[F],
    idx: Option<usize>,
    thread_depth: Option<usize>,
) where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + Send + Sync,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    if v.len() == 1 {
        return;
    }
    let idx = idx.unwrap_or_else(|| 1);

    let thread_depth =
        thread_depth.unwrap_or_else(|| rayon::current_num_threads().ilog2() as usize);

    fft_mul(v, twiddles[idx - 1]);
    let (u, w) = split_half(v);

    if thread_depth > 0 {
        // Parallel case
        rayon::join(
            || fft_twiddles_parallel(u, twiddles, Some(idx * 2), Some(thread_depth - 1)),
            || fft_twiddles_parallel(w, twiddles, Some(idx * 2 + 1), Some(thread_depth - 1)),
        );
    } else {
        // Serial case - fall back to regular fft_twiddles
        fft_twiddles(u, twiddles, Some(idx * 2));
        fft_twiddles(w, twiddles, Some(idx * 2 + 1));
    }
}

fn fft_twiddles<F>(v: &mut [F], twiddles: &[F], idx: Option<usize>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
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
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
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

fn fft_mul<F>(v: &mut [F], lambda: F)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let (u, w) = split_half(v);

    for i in 0..u.len() {
        u[i] = u[i] + lambda * w[i];
        w[i] = w[i] + u[i];
    }
}

fn ifft_mul<F>(v: &mut [F], lambda: F)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let (u, w) = split_half(v);

    for i in 0..u.len() {
        w[i] = w[i] + u[i];
        u[i] = u[i] + lambda * w[i];
    }
}

fn split_half<F>(v: &mut [F]) -> (&mut [F], &mut [F]) {
    let mid = v.len() / 2;
    v.split_at_mut(mid)
}

pub fn is_power_of_2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

pub fn eval_sk_at_vks<F>(n: usize) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(is_power_of_2(n));
    let num_subspaces = n.ilog2() as usize;

    let mut sk_vks: Vec<F> = Vec::with_capacity(num_subspaces + 1);
    sk_vks.push(F::one());

    let mut layer: Vec<F> = (1..=num_subspaces)
        .map(|i| F::new(F::ValueType::try_from(1 << i).unwrap()))
        .collect();
    let mut cur_len = num_subspaces;

    for i in 1..=num_subspaces {
        for j in 1..=cur_len {
            if j == 1 {
                sk_vks.push(layer[0] * layer[0] + sk_vks[i - 1] * layer[0]);
            } else {
                layer[j - 2] = layer[j - 1] * layer[j - 1] + sk_vks[i - 1] * layer[j - 1];
            }
        }
        cur_len -= 1;
    }

    sk_vks
}

pub fn evaluate_basis<F>(basis_len: usize, sks_vks: &[F], x: F) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(is_power_of_2(basis_len));
    let num_subspaces = basis_len.ilog2() as usize;

    let mut sks_at_x: Vec<F> = Vec::with_capacity(num_subspaces);
    sks_at_x.push(x);

    for i in 2..=num_subspaces {
        sks_at_x.push(next_s(sks_at_x[i - 2], sks_vks[i - 2]));
    }

    let mut basis: Vec<F> = vec![F::zero(); basis_len];
    basis[0] = F::one();
    for i in 1..=num_subspaces {
        let current_len: usize = 1 << (i - 1);
        for j in 1..=current_len {
            basis[j + current_len - 1] = sks_at_x[i - 1] * basis[j - 1];
        }
    }

    basis
}

pub fn evaluate_scaled_basis_inplace<F>(
    sks_at_x: &mut [F],
    basis: &mut [F],
    sks_vks: &[F],
    x: F,
    alpha: F,
) where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let num_subspaces = basis.len().ilog2() as usize;

    sks_at_x[0] = x;
    for i in 2..=num_subspaces {
        sks_at_x[i - 1] = next_s(sks_at_x[i - 2], sks_vks[i - 2]);
    }

    basis[0] = alpha;
    for i in 1..=num_subspaces {
        let current_len: usize = 1 << (i - 1);
        for j in 1..=current_len {
            basis[j + current_len - 1] = sks_at_x[i - 1] * basis[j - 1];
        }
    }
}

pub fn compute_pis<F>(pis_len: usize, sks_vks: &[F]) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    assert!(pis_len == 1 << (sks_vks.len() - 1));

    let mut pis: Vec<F> = vec![F::zero(); pis_len];
    pis[0] = F::one();

    for i in 2..=sks_vks.len() {
        let current_len: usize = 1 << (i - 2);
        for j in 1..=current_len {
            pis[j + current_len - 1] = sks_vks[i - 2] * pis[j - 1];
        }
    }

    pis
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

    #[test]
    fn test_fft_ifft_roundtrip() {
        // Test that ifft(fft(v)) == v
        let log_n = 4; // Small test: 2^4 = 16 elements
        let twiddles = compute_twiddles::<BinaryElem16>(log_n, None);

        // Create a test vector with some values
        let mut v = vec![
            BinaryElem16::new(1),
            BinaryElem16::new(2),
            BinaryElem16::new(3),
            BinaryElem16::new(4),
            BinaryElem16::new(5),
            BinaryElem16::new(6),
            BinaryElem16::new(7),
            BinaryElem16::new(8),
            BinaryElem16::new(9),
            BinaryElem16::new(10),
            BinaryElem16::new(11),
            BinaryElem16::new(12),
            BinaryElem16::new(13),
            BinaryElem16::new(14),
            BinaryElem16::new(15),
            BinaryElem16::new(16),
        ];
        let original = v.clone();

        println!(
            "Original: {:?}",
            v.iter().map(|x| x.value).collect::<Vec<_>>()
        );

        // Apply FFT
        fft(&mut v, &twiddles, None);
        println!(
            "After FFT: {:?}",
            v.iter().map(|x| x.value).collect::<Vec<_>>()
        );

        // Apply IFFT
        ifft(&mut v, &twiddles);
        println!(
            "After IFFT: {:?}",
            v.iter().map(|x| x.value).collect::<Vec<_>>()
        );

        // Should get back original
        assert_eq!(v, original, "FFT->IFFT should be identity transform");
    }
}
