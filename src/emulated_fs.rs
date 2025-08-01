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
        F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
        F::ValueType: TryFrom<u128>,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_field::BinaryElem16;

    #[test]
    fn test_fs_creation() {
        let fs1 = FS::new(12345);
        let fs2 = FS::new(12345);

        // Same seed should create identical initial state
        assert_eq!(fs1.counter, fs2.counter);
        assert_eq!(fs1.counter, 0);
    }

    #[test]
    fn test_deterministic_behavior() {
        let mut fs1 = FS::new(42);
        let mut fs2 = FS::new(42);

        // Same operations should produce same results
        let val1 = fs1.squeeze();
        let val2 = fs2.squeeze();
        assert_eq!(val1, val2);

        let val3 = fs1.squeeze();
        let val4 = fs2.squeeze();
        assert_eq!(val3, val4);

        // But different calls should produce different values
        assert_ne!(val1, val3);
    }

    #[test]
    fn test_absorb_affects_output() {
        let mut fs1 = FS::new(1234);
        let mut fs2 = FS::new(1234);

        // Absorb different data
        fs1.absorb(b"hello");
        fs2.absorb(b"world");

        // Should produce different outputs
        let val1 = fs1.squeeze();
        let val2 = fs2.squeeze();
        assert_ne!(val1, val2);
    }

    #[test]
    fn test_absorb_order_matters() {
        let mut fs1 = FS::new(1234);
        let mut fs2 = FS::new(1234);

        // Absorb same data in different order
        fs1.absorb(b"hello");
        fs1.absorb(b"world");

        fs2.absorb(b"world");
        fs2.absorb(b"hello");

        // Should produce different outputs
        let val1 = fs1.squeeze();
        let val2 = fs2.squeeze();
        assert_ne!(val1, val2);
    }

    #[test]
    fn test_state_preservation() {
        let mut fs = FS::new(5678);

        // Absorb some data
        fs.absorb(b"initial data");
        let val1 = fs.squeeze();

        // Absorb more data - should build on previous state
        fs.absorb(b"additional data");
        let val2 = fs.squeeze();

        // Create fresh FS and try to get val2 directly
        let mut fs_fresh = FS::new(5678);
        fs_fresh.absorb(b"initial data");
        fs_fresh.squeeze(); // consume first squeeze
        fs_fresh.absorb(b"additional data");
        let val2_fresh = fs_fresh.squeeze();

        assert_eq!(val2, val2_fresh);
    }

    #[test]
    fn test_get_field() {
        let mut fs = FS::new(9999);

        // Should generate field elements
        let field1: BinaryElem16 = fs.get_field();
        let field2: BinaryElem16 = fs.get_field();

        // Should be different values
        assert_ne!(field1.value, field2.value);

        // Should be deterministic
        let mut fs2 = FS::new(9999);
        let field1_repeat: BinaryElem16 = fs2.get_field();
        let field2_repeat: BinaryElem16 = fs2.get_field();

        assert_eq!(field1.value, field1_repeat.value);
        assert_eq!(field2.value, field2_repeat.value);
    }

    #[test]
    fn test_get_query() {
        let mut fs = FS::new(1111);

        // Test with various ranges
        for n in [1, 5, 10, 100] {
            for _ in 0..20 {
                let query = fs.get_query(n);
                assert!(query < n, "Query {} should be less than {}", query, n);
            }
        }
    }

    #[test]
    fn test_get_query_deterministic() {
        let mut fs1 = FS::new(2222);
        let mut fs2 = FS::new(2222);

        // Same inputs should give same outputs
        for n in [3, 7, 15] {
            let q1 = fs1.get_query(n);
            let q2 = fs2.get_query(n);
            assert_eq!(q1, q2);
        }
    }

    #[test]
    fn test_get_distinct_queries() {
        let mut fs = FS::new(3333);

        let n = 10;
        let s = 5;
        let queries = fs.get_distinct_queries(n, s);

        // Should return exactly s queries
        assert_eq!(queries.len(), s);

        // All should be in valid range
        for &q in &queries {
            assert!(q < n, "Query {} should be less than {}", q, n);
        }

        // All should be unique
        let unique_count = queries.iter().collect::<HashSet<_>>().len();
        assert_eq!(unique_count, s);
    }

    #[test]
    fn test_get_distinct_queries_edge_cases() {
        let mut fs = FS::new(4444);

        // Request all possible values
        let n = 5;
        let s = 5;
        let queries = fs.get_distinct_queries(n, s);

        assert_eq!(queries.len(), s);
        let unique_count = queries.iter().collect::<HashSet<_>>().len();
        assert_eq!(unique_count, s);

        // All values from 0 to n-1 should be present
        let mut sorted_queries = queries.clone();
        sorted_queries.sort();
        assert_eq!(sorted_queries, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_get_distinct_queries_deterministic() {
        let mut fs1 = FS::new(5555);
        let mut fs2 = FS::new(5555);

        let queries1 = fs1.get_distinct_queries(20, 8);
        let queries2 = fs2.get_distinct_queries(20, 8);

        // Should produce same results
        assert_eq!(queries1, queries2);
    }

    #[test]
    fn test_counter_increments() {
        let mut fs = FS::new(6666);

        assert_eq!(fs.counter, 0);

        fs.squeeze();
        assert_eq!(fs.counter, 1);

        fs.squeeze();
        assert_eq!(fs.counter, 2);

        // get_field should also increment counter
        let _: BinaryElem16 = fs.get_field();
        assert_eq!(fs.counter, 3);

        // get_query should also increment counter
        fs.get_query(10);
        assert_eq!(fs.counter, 4);
    }

    #[test]
    fn test_different_seeds_different_output() {
        let mut fs1 = FS::new(111);
        let mut fs2 = FS::new(222);

        let val1 = fs1.squeeze();
        let val2 = fs2.squeeze();

        assert_ne!(val1, val2);

        let field1: BinaryElem16 = fs1.get_field();
        let field2: BinaryElem16 = fs2.get_field();

        assert_ne!(field1.value, field2.value);
    }

    #[test]
    fn test_integration_fiat_shamir_pattern() {
        // Simulate a typical Fiat-Shamir transcript
        let mut fs = FS::new(12345);

        // Prover commits to some values
        fs.absorb(b"commitment_1");
        fs.absorb(b"commitment_2");

        // Extract challenge
        let challenge1: BinaryElem16 = fs.get_field();

        // Prover responds with more data
        fs.absorb(&challenge1.to_bytes());
        fs.absorb(b"response_data");

        // Extract another challenge
        let challenge2: BinaryElem16 = fs.get_field();

        // Get some query indices
        let queries = fs.get_distinct_queries(100, 10);

        // Verify deterministic behavior by repeating
        let mut fs2 = FS::new(12345);
        fs2.absorb(b"commitment_1");
        fs2.absorb(b"commitment_2");
        let challenge1_repeat: BinaryElem16 = fs2.get_field();
        fs2.absorb(&challenge1_repeat.to_bytes());
        fs2.absorb(b"response_data");
        let challenge2_repeat: BinaryElem16 = fs2.get_field();
        let queries_repeat = fs2.get_distinct_queries(100, 10);

        assert_eq!(challenge1.value, challenge1_repeat.value);
        assert_eq!(challenge2.value, challenge2_repeat.value);
        assert_eq!(queries, queries_repeat);
    }
}
