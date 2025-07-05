use crate::binary_fft::{compute_pis, compute_twiddles, eval_sk_at_vks};

use crate::binary_fft::is_power_of_2;
use crate::binary_field::BinaryField;
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

impl<F> ReedSolomonEncoding<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
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
}
