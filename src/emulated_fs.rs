use crate::binary_field::BinaryField;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::ops::{Add, Mul};

pub struct FS {
    hasher: Sha256,
    counter: u32,
}

impl FS {
    pub fn new(seed: u32) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(&seed.to_le_bytes());
        Self { hasher, counter: 0 }
    }

    pub fn absorb(&mut self, elems: &[u8]) {
        self.hasher.update(elems);
    }

    pub fn squeeze(&mut self) -> u128 {
        self.hasher.update(&self.counter.to_le_bytes());
        self.counter += 1;
        let hasher_copy = self.hasher.clone();
        let result = hasher_copy.finalize();
        u128::from_le_bytes([
            result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7],
            result[8], result[9], result[10], result[11], result[12], result[13], result[14],
            result[15],
        ])
    }

    pub fn get_field<F>(&mut self) -> F
    where
        F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
        F::ValueType: TryFrom<u128> + Copy,
        <F::ValueType as TryFrom<u128>>::Error: std::fmt::Debug,
    {
        let value = self.squeeze();
        // For binary-field elements smaller than 64 bits, so we need to mask to fit
        // This should be cryptographically-safe considering we just want randomness
        match std::mem::size_of::<F::ValueType>() {
            1 => F::new(F::ValueType::try_from(value as u8 as u128).unwrap()),
            2 => F::new(F::ValueType::try_from(value as u16 as u128).unwrap()),
            4 => F::new(F::ValueType::try_from(value as u32 as u128).unwrap()),
            8 => F::new(F::ValueType::try_from(value as u64 as u128).unwrap()),
            16 => F::new(F::ValueType::try_from(value).unwrap()),
            _ => panic!(
                "Unsupported field size {}",
                std::mem::size_of::<F::ValueType>()
            ),
        }
    }

    pub fn get_query(&mut self, n: usize) -> usize {
        let value = self.squeeze() as usize;
        value % n
    }

    pub fn get_distinct_queries(&mut self, n: usize, s: usize) -> Vec<usize> {
        let mut queries: Vec<usize> = Vec::with_capacity(s);
        let mut seen: HashSet<usize> = HashSet::new();

        while queries.len() < s {
            let q = self.get_query(n);
            if !seen.contains(&q) {
                queries.push(q);
                seen.insert(q);
            }
        }

        queries
    }
}
