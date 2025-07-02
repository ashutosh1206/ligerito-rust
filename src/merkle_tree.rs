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

fn ith_layer(
    current_layer: &Vec<u8>,
    queries_len: usize,
    queries: &mut Vec<usize>,
    proof: &mut Vec<Vec<u8>>,
) -> usize {
    let mut next_queries_len: usize = 0;
    let mut i = 0;

    while i < queries_len {
        let query = queries[i];
        let sibling = query ^ 1;

        next_queries_len += 1;
        queries[next_queries_len - 1] = query >> 1;

        let sibling_hash = &current_layer[sibling * 32..sibling * 32 + 32];
        if i == queries_len - 1 {
            proof.push(sibling_hash.to_vec());
            break;
        }

        if query % 2 != 0 {
            proof.push(sibling_hash.to_vec());
            i += 1;
        } else if queries[i + 1] != sibling {
            proof.push(sibling_hash.to_vec());
            i += 1;
        } else {
            i += 2;
        }
    }

    next_queries_len
}

pub fn prove(tree: Vec<Vec<u8>>, queries: Vec<usize>) -> Vec<Vec<u8>> {
    let mut proof: Vec<Vec<u8>> = Vec::new();
    let depth = tree.len() - 1;

    let mut queries_buff = queries.clone();
    let mut queries_cnt = queries.len();

    for i in 0..depth {
        let current_layer = &tree[i];
        queries_cnt = ith_layer(current_layer, queries_cnt, &mut queries_buff, &mut proof);
    }

    proof
}

fn hash_siblings(left: &[u8], right: &[u8]) -> Vec<u8> {
    let mut combined = Vec::with_capacity(64);
    combined.extend_from_slice(left);
    combined.extend_from_slice(right);
    Sha256::digest(&combined).to_vec()
}

fn verify_ith_layer(
    layer: &mut Vec<Vec<u8>>,
    queries: &mut Vec<usize>,
    curr_cnt: usize,
    proof: &Vec<Vec<u8>>,
    mut proof_cnt: usize,
) -> (usize, usize) {
    let mut next_cnt = 0;
    let mut i = 0;

    while i < curr_cnt {
        let query = queries[i];
        let sibling = query ^ 1;

        next_cnt += 1;
        queries[next_cnt - 1] = query >> 1;

        if i == curr_cnt - 1 {
            proof_cnt += 1;
            if query % 2 != 0 {
                layer[next_cnt - 1] = hash_siblings(&proof[proof_cnt - 1], &layer[i]);
            } else {
                layer[next_cnt - 1] = hash_siblings(&layer[i], &proof[proof_cnt - 1]);
            }
            break;
        }

        if query % 2 != 0 {
            proof_cnt += 1;
            layer[next_cnt - 1] = hash_siblings(&proof[proof_cnt - 1], &layer[i]);
            i += 1;
        } else if queries[i + 1] != sibling {
            proof_cnt += 1;
            layer[next_cnt - 1] = hash_siblings(&layer[i], &proof[proof_cnt - 1]);
            i += 1;
        } else {
            layer[next_cnt - 1] = hash_siblings(&layer[i], &layer[i + 1]);
            i += 2;
        }
    }

    (next_cnt, proof_cnt)
}

pub fn verify<T>(
    root: Vec<u8>,
    batched_proof: Vec<Vec<u8>>,
    depth: usize,
    leaves: Vec<T>,
    leaf_indices: Vec<usize>,
) -> bool
where
    T: AsRef<[u8]>,
{
    let mut proof = batched_proof.clone();
    let mut layer: Vec<Vec<u8>> = leaves
        .iter()
        .map(|leaf| Sha256::digest(leaf).to_vec())
        .collect();
    let mut queries = leaf_indices.clone();
    let mut curr_cnt = queries.len();
    let mut proof_cnt = 0usize;

    for _ in 0..depth {
        (curr_cnt, proof_cnt) =
            verify_ith_layer(&mut layer, &mut queries, curr_cnt, &proof, proof_cnt);
    }

    curr_cnt == 1 && layer[0] == root
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

    #[test]
    fn test_prove_sparse_queries() {
        let leaves = vec!["hello", "world", "rust", "test"];
        let tree = build_merkle_tree(leaves);

        // Query leaves at indices 0 and 2 (0-indexed, corresponds to Julia's [1, 3])
        let queries = vec![0, 2];
        let proof = prove(tree, queries);

        // Expected proof from actual Julia MerkleTree implementation
        let expected_proof = vec![
            vec![
                72, 110, 164, 98, 36, 209, 187, 79, 182, 128, 243, 79, 124, 154, 217, 106, 143, 36,
                236, 136, 190, 115, 234, 142, 90, 108, 101, 38, 14, 156, 184, 167,
            ],
            vec![
                159, 134, 208, 129, 136, 76, 125, 101, 154, 47, 234, 160, 197, 90, 208, 21, 163,
                191, 79, 27, 43, 11, 130, 44, 209, 93, 108, 21, 176, 240, 10, 8,
            ],
        ];

        assert_eq!(proof.len(), 2);
        assert_eq!(proof, expected_proof);
    }

    #[test]
    fn test_prove_consecutive_queries() {
        let leaves = vec!["hello", "world", "rust", "test"];
        let tree = build_merkle_tree(leaves);

        // Query leaves at indices 0, 1, and 2 (0-indexed, corresponds to Julia's [1, 2, 3])
        let queries = vec![0, 1, 2];
        let proof = prove(tree, queries);

        // Expected proof from actual Julia MerkleTree implementation
        let expected_proof = vec![vec![
            159, 134, 208, 129, 136, 76, 125, 101, 154, 47, 234, 160, 197, 90, 208, 21, 163, 191,
            79, 27, 43, 11, 130, 44, 209, 93, 108, 21, 176, 240, 10, 8,
        ]];

        assert_eq!(proof.len(), 1);
        assert_eq!(proof, expected_proof);
    }

    #[test]
    fn test_prove_and_verify_roundtrip() {
        let leaves = vec!["hello", "world", "rust", "test"];
        let tree = build_merkle_tree(leaves.clone());

        // Get the root hash
        let root = tree.last().unwrap().clone();
        let depth = tree.len() - 1;

        // Test case 1: Sparse queries [0, 2]
        let queries = vec![0, 2];
        let proof = prove(tree.clone(), queries.clone());
        let queried_leaves = vec![leaves[0], leaves[2]];

        let is_valid = verify(root.clone(), proof, depth, queried_leaves, queries);
        assert!(is_valid, "Verification should succeed for sparse queries");

        // Test case 2: Consecutive queries [0, 1, 2]
        let queries2 = vec![0, 1, 2];
        let proof2 = prove(tree.clone(), queries2.clone());
        let queried_leaves2 = vec![leaves[0], leaves[1], leaves[2]];

        let is_valid2 = verify(root, proof2, depth, queried_leaves2, queries2);
        assert!(
            is_valid2,
            "Verification should succeed for consecutive queries"
        );
    }
}
