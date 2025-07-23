use crate::binary_field::BinaryField;
use std::ops::{Add, Mul};

pub struct RecursiveLigeroWitness<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
{
    pub flat_mat: Vec<F>,
    pub tree: Vec<Vec<u8>>,
}
