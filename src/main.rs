use cryptoutils::data_structures::ProofSize;
use cryptoutils::*;
use std::time::Instant;

fn format_bytes(bytes: usize) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut unit_index = 0;

    while size >= 1024.0 && unit_index < UNITS.len() - 1 {
        size /= 1024.0;
        unit_index += 1;
    }

    if unit_index == 0 {
        format!("{} {}", bytes, UNITS[unit_index])
    } else {
        format!("{:.1} {}", size, UNITS[unit_index])
    }
}

fn main() {
    let total_start = Instant::now();
    // Same test as test_prove_verify_config_20 but as a binary
    let poly: Vec<BinaryElem32> = (0..(1 << 20))
        .map(|_| BinaryElem32::new(rand::random()))
        .collect();

    let prover_config: ProverConfig<BinaryElem32> = hardcoded_config_20::<BinaryElem32>();
    let verifier_config = hardcoded_config_20_verifier();

    // Use original prover directly
    let prover_start = Instant::now();
    let proof = prover(prover_config, poly);
    let prover_duration = prover_start.elapsed();

    let proof_size = proof.proof_size();
    println!("Proof size: {}", format_bytes(proof_size));

    let verifier_start = Instant::now();
    let verification_result = verifier(verifier_config, proof);
    let verifier_duration = verifier_start.elapsed();

    let total_duration = total_start.elapsed();
    println!("=== Performance Results ===");
    println!("Prover time:    {:>8.2?}", prover_duration);
    println!("Verifier time:  {:>8.2?}", verifier_duration);
    println!("Total time:     {:>8.2?}", total_duration);
    println!("===========================");

    println!("Verification result: {}", verification_result);
    assert!(verification_result, "Verification should pass");
}
