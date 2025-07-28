use crate::binary_field::BinaryField;
use crate::reed_solomon::ReedSolomonEncoding;
use std::ops::{Add, Mul};

pub struct ProverConfig<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    recursive_steps: usize,
    initial_dims: (usize, usize),
    dims: Vec<(usize, usize)>,
    initial_k: usize,
    ks: Vec<usize>,
    initial_reed_solomon: ReedSolomonEncoding<F>,
    reed_solomon_codes: Vec<ReedSolomonEncoding<F>>,
}

pub struct RecursiveLigeroWitness<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
{
    pub flat_mat: Vec<F>,
    pub tree: Vec<Vec<u8>>,
}
