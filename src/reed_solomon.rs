use crate::binary_fft::is_power_of_2;
use crate::binary_fft::{compute_pis, compute_twiddles, eval_sk_at_vks};
use crate::binary_field::BinaryField;
use crate::{fft, ifft};
use std::ops::{Add, Mul};

pub struct ReedSolomonEncoding<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    log_message_length: usize,
    log_block_length: usize,
    twiddles: Vec<F>,
    pis: Vec<F>,
}

fn short_from_long_tw<F>(l_tw: &[F], n: usize, k: usize) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut s_tw = vec![F::zero(); (1 << k) - 1];
    let mut jump: usize = 1 << (n - k);
    s_tw[0] = l_tw[jump - 1];

    let mut idx: usize = 2;
    for i in 1..k {
        jump *= 2;
        let take: usize = 1 << i;

        for j in 0..take {
            s_tw[idx + j - 1] = l_tw[jump + j - 1];
        }
        idx += take;
    }

    s_tw
}

impl<F> ReedSolomonEncoding<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + Send + Sync,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(message_length: usize, block_length: usize) -> Self {
        // See if you can combine these two asserts using &&
        assert!(is_power_of_2(message_length));
        assert!(is_power_of_2(block_length));
        assert!(message_length < block_length);

        let log_message_length = message_length.ilog2() as usize;
        let log_block_length = block_length.ilog2() as usize;

        Self {
            log_block_length,
            log_message_length,
            twiddles: compute_twiddles::<F>(log_block_length, None),
            pis: compute_pis(message_length, &eval_sk_at_vks::<F>(message_length)),
        }
    }

    pub fn message_length(&self) -> usize {
        1 << self.log_message_length
    }

    pub fn block_length(&self) -> usize {
        1 << self.log_block_length
    }

    pub fn log_message_length(&self) -> usize {
        self.log_message_length
    }

    pub fn log_block_length(&self) -> usize {
        self.log_block_length
    }

    fn short_from_long_twiddles(&self) -> Vec<F> {
        short_from_long_tw(
            &self.twiddles,
            self.log_block_length,
            self.log_message_length,
        )
    }

    pub fn encode(&self, message: &[F]) -> Vec<F> {
        assert!(message.len() == self.message_length());
        let mut message_coeffs = vec![F::zero(); self.block_length()];
        message_coeffs[0..self.message_length()].copy_from_slice(message);

        let s_tw = self.short_from_long_twiddles();
        ifft(&mut message_coeffs[0..self.message_length()], &s_tw);
        fft(&mut message_coeffs, &self.twiddles);

        message_coeffs
    }

    pub fn encode_non_systematic(&self, message: &mut [F]) {
        assert!(message.len() == self.block_length());

        for i in 0..self.message_length() {
            message[i] = message[i] * self.pis[i];
        }

        fft(message, &self.twiddles);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_field::BinaryElem16;

    #[test]
    fn test_reed_solomon_systematic() {
        // Mirroring the test case from corresponding Julia implementation
        let rs = ReedSolomonEncoding::<BinaryElem16>::new(1024, 4096); // 2^10, 2^12

        // Create a simple test message instead of random
        let message: Vec<BinaryElem16> = (1..=1024).map(|i| BinaryElem16::new(i as u16)).collect();

        let encoded = rs.encode(&message);
        assert_eq!(message, encoded[0..message.len()]);
    }

    #[test]
    fn test_reed_solomon_non_systematic() {
        use crate::binary_fft::{eval_sk_at_vks, evaluate_basis};

        let rs = ReedSolomonEncoding::<BinaryElem16>::new(1024, 4096); // 2^10, 2^12

        // Create test vectors (simulating a_cfs and b_cfs)
        let a_cfs: Vec<BinaryElem16> = (1..=1024).map(|i| BinaryElem16::new(i as u16)).collect();
        let b_cfs: Vec<BinaryElem16> = (1025..=2048).map(|i| BinaryElem16::new(i as u16)).collect();

        // Linear combination parameter
        let l = BinaryElem16::new(42);

        // Compute c_cfs = a_cfs + l * b_cfs
        let c_cfs: Vec<BinaryElem16> = a_cfs
            .iter()
            .zip(b_cfs.iter())
            .map(|(&a, &b)| a + l * b)
            .collect();

        // Pad to block length for non-systematic encoding
        let mut a_padded = vec![BinaryElem16::zero(); 4096];
        let mut b_padded = vec![BinaryElem16::zero(); 4096];
        a_padded[0..1024].copy_from_slice(&a_cfs);
        b_padded[0..1024].copy_from_slice(&b_cfs);

        // Encode non-systematically
        rs.encode_non_systematic(&mut a_padded);
        rs.encode_non_systematic(&mut b_padded);

        // Test random row query - using a fixed index instead of random
        let x = 100usize;
        let xf = BinaryElem16::new(x as u16);

        // Compute expected value using basis evaluation
        let sks_vks = eval_sk_at_vks::<BinaryElem16>(1024);
        let basis = evaluate_basis(1024, &sks_vks, xf);
        let c_at_x: BinaryElem16 = c_cfs
            .iter()
            .zip(basis.iter())
            .map(|(&c, &b)| c * b)
            .fold(BinaryElem16::zero(), |acc, val| acc + val);

        // Compute result from encoded values
        let a_row_opening = a_padded[x];
        let b_row_opening = b_padded[x];
        let result = a_row_opening + l * b_row_opening;

        assert_eq!(result, c_at_x);
    }
}
