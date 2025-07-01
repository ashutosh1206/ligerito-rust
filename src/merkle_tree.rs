use crate::binary_fft::is_power_of_2;
use sha2::{Digest, Sha256};

fn hash_array<T>(x: Vec<T>) -> Vec<u8>
where
    T: AsRef<[u8]>,
{
    let n = x.len();

    // Each element in the vector has a 32 byte hash, each of
    // which is stored as an element in xh
    let mut xh: Vec<u8> = Vec::with_capacity(32 * n);

    for leaf in x {
        let hash = Sha256::digest(leaf);
        xh.extend_from_slice(&hash);
    }

    xh
}

pub fn build_merkle_tree<T>(leaves: Vec<T>) -> Vec<Vec<u8>>
where
    T: AsRef<[u8]>,
{
    assert!(is_power_of_2(leaves.len()));
    if leaves.is_empty() {
        return Vec::new();
    }

    let n = leaves.len();
    let lh = hash_array(leaves);
    let mut layers = vec![lh];

    let mut current_layer = layers[0].clone();
    let mut current_n = n;

    while current_n > 1 {
        let next_n = (current_n + 1) / 2;
        let mut next_layer: Vec<u8> = Vec::with_capacity(32 * next_n);

        for i in 0..next_n {
            let left_hash = &current_layer[i * 2 * 32..(i * 2 * 32) + 32];
            let right_hash = &current_layer[(i * 2 + 1) * 32..((i * 2 + 1) * 32) + 32];

            let mut combined = Vec::with_capacity(64);
            combined.extend_from_slice(left_hash);
            combined.extend_from_slice(right_hash);

            let hash = Sha256::digest(&combined);
            next_layer.extend_from_slice(&hash);
        }

        layers.push(next_layer.clone());
        current_layer = next_layer;
        current_n = next_n;
    }

    return layers;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_merkle_tree() {
        let leaves = vec!["hello", "world", "rust", "test"];
        let layers = build_merkle_tree(leaves);

        assert_eq!(layers.len(), 3);
        assert_eq!(layers[0].len(), 128);
        assert_eq!(layers[1].len(), 64);
        assert_eq!(layers[2].len(), 32);

        // root hash validates entire tree construction
        // Expected value created from output of actual Julia's build_merkle_tree
        let expected_root = vec![
            27, 31, 141, 106, 120, 231, 92, 57, 126, 100, 218, 70, 17, 208, 139, 253, 52, 148, 110,
            160, 229, 27, 10, 236, 193, 146, 31, 236, 141, 245, 76, 152,
        ];
        assert_eq!(layers[2], expected_root);
    }
}
