use crate::binary_field::BinaryField;
use crate::data_structures::ProverConfig;
use std::ops::{Add, Mul};

pub fn prover<F>(config: ProverConfig<F>, poly: Vec<F>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
}
