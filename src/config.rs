use crate::binary_field::BinaryField;
use crate::data_structures::{ProverConfig, VerifierConfig};
use crate::reed_solomon::ReedSolomonEncoding;
use std::ops::{Add, Mul};

pub fn hardcoded_config_20<F>() -> ProverConfig<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + From<u128> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let recursive_steps = 1;
    let inv_rate = 4;

    let mut dims: Vec<(usize, usize)> = Vec::with_capacity(recursive_steps);
    let initial_dims: (usize, usize) = (1 << 14, 1 << 6);
    dims.push((1 << 10, 1 << 4));

    let mut ks: Vec<usize> = Vec::with_capacity(recursive_steps);
    let initial_k = 6;
    ks.push(4);

    let initial_reed_solomon: ReedSolomonEncoding<F> =
        ReedSolomonEncoding::new(initial_dims.0, initial_dims.0 * inv_rate);
    let mut reed_solomon_codes: Vec<ReedSolomonEncoding<F>> = Vec::with_capacity(recursive_steps);

    for i in 0..recursive_steps {
        reed_solomon_codes.push(ReedSolomonEncoding::new(dims[i].0, dims[i].0 * inv_rate));
    }

    ProverConfig {
        recursive_steps,
        initial_dims,
        dims,
        initial_k,
        ks,
        initial_reed_solomon,
        reed_solomon_codes,
    }
}

pub fn hardcoded_config_20_verifier() -> VerifierConfig {
    let recursive_steps = 1;

    let mut log_rs_dims: Vec<usize> = Vec::with_capacity(recursive_steps);
    let initial_dim = 14;
    log_rs_dims.push(10);

    let mut ks = Vec::with_capacity(recursive_steps);
    let initial_k = 6;
    ks.push(4);

    VerifierConfig {
        recursive_steps,
        initial_dim,
        log_dims: log_rs_dims,
        initial_k,
        ks,
    }
}
