use crate::binary_field::BinaryField;
use crate::data_structures::{ProverConfig, VerifierConfig};
use crate::reed_solomon::ReedSolomonEncoding;
use std::ops::{Add, Mul};

pub fn hardcoded_config_20<F>() -> ProverConfig<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize> + TryFrom<u128> + Copy,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_field::BinaryElem32;
    use crate::prover::prover;
    use crate::verifier::verifier;

    #[test]
    fn test_ligerito_full() {
        let config = hardcoded_config_20::<BinaryElem32>();
    }

    #[test]
    fn test_merkle_with_ligero_data() {
        use crate::merkle_tree::*;
        use crate::ligero::*;
        use crate::binary_field::BinaryElem32;
        
        // Create test polynomial - need m*n elements
        let m = 16384;
        let n = 64;
        let poly: Vec<BinaryElem32> = (0..(m*n)).map(|i| BinaryElem32::new(i as u32)).collect();
        let prover_config = hardcoded_config_20::<BinaryElem32>();
        
        // Test ligero_commit
        let wtns = ligero_commit(&poly, prover_config.initial_dims.0, prover_config.initial_dims.1, &prover_config.initial_reed_solomon);
        
        println!("Ligero commit test:");
        println!("  Tree depth: {}", wtns.tree.len() - 1);
        println!("  Num rows: {}", wtns.num_rows);
        println!("  Num cols: {}", wtns.num_cols);
        
        // Test with a few specific queries AND some random ones
        let queries = vec![0, 1, 100, 50459, 50294, 1959];
        let opened_rows: Vec<Vec<BinaryElem32>> = queries
            .iter()
            .map(|&q| extract_row(&wtns.flat_mat, q, wtns.num_rows, wtns.num_cols))
            .collect();
            
        // Convert to bytes the same way as in the verifier
        let opened_rows_bytes: Vec<Vec<u8>> = opened_rows
            .iter()
            .map(|row| row.iter().flat_map(|elem| elem.to_bytes()).collect())
            .collect();
            
        let proof = prove(&wtns.tree, &queries);
        let root = wtns.tree.last().unwrap();
        let depth = wtns.tree.len() - 1;
        
        println!("  Proof length: {}", proof.len());
        println!("  Root: {:?}", &root[..16]);
        
        let is_valid = verify(root, &proof, depth, &opened_rows_bytes, &queries);
        println!("  Verification result: {}", is_valid);
        
        if !is_valid {
            println!("  First query: {}", queries[0]);
            println!("  First opened row length: {}", opened_rows_bytes[0].len());
            println!("  First opened row first 16 bytes: {:?}", &opened_rows_bytes[0][..16]);
            println!("  Proof first element: {:?}", &proof[0][..16]);
        }
        
        assert!(is_valid, "Ligero Merkle verification should pass");
    }

    #[test]
    fn test_prove_verify_config_20() {
        // Julia: poly = rand(BinaryElem32, 2^20)
        let poly: Vec<BinaryElem32> = (0..(1 << 20))
            .map(|_| BinaryElem32::new(rand::random()))
            .collect();

        let prover_config: ProverConfig<BinaryElem32> = hardcoded_config_20::<BinaryElem32>();
        let verifier_config = hardcoded_config_20_verifier();

        // Run prover
        let proof = prover(prover_config, poly);
        
        // Run verifier
        let verification_result = verifier(verifier_config, proof);

        assert!(verification_result, "Verification should pass");
    }
}
