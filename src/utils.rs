use crate::binary_field::BinaryField;
use std::{
    ops::{Add, Mul},
    thread::current,
};

fn evaluate_lagrange_basis<F>(rs: Vec<F>) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
{
    let one_elem = F::one();
    let mut current_layer = vec![one_elem + rs[0], rs[0]];
    let mut layer_len: usize = 2;

    for i in 1..rs.len() {
        let next_layer_size = 2 * layer_len;
        let mut next_layer = vec![F::zero(); next_layer_size];

        let ri_p_one = one_elem + rs[i];
        for j in 0..layer_len {
            next_layer[2 * j] = current_layer[j] * ri_p_one;
            next_layer[2 * j + 1] = current_layer[j] * rs[i];
        }

        current_layer = next_layer;
        layer_len = next_layer_size;
    }

    current_layer
}
