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

    pub fn squeeze(&mut self) -> usize {
        self.hasher.update(&self.counter.to_le_bytes());
        self.counter += 1;
        let hasher_copy = self.hasher.clone();
        let result = hasher_copy.finalize();
        u64::from_le_bytes([
            result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7],
        ]) as usize
    }

    pub fn get_field<F>(&mut self) -> F
    where
        F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
        F::ValueType: TryFrom<usize> + Copy,
        <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
    {
        let value = self.squeeze();
        F::new(F::ValueType::try_from(value).unwrap())
    }

    pub fn get_query(&mut self, n: usize) -> usize {
        let value = self.squeeze();
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
