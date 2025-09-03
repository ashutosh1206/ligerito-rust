use cryptoutils::*;

fn main() {
    // Same test as test_prove_verify_config_20 but as a binary
    let poly: Vec<BinaryElem32> = (0..(1 << 20))
        .map(|_| BinaryElem32::new(rand::random()))
        .collect();

    let prover_config: ProverConfig<BinaryElem32> = hardcoded_config_20::<BinaryElem32>();
    let verifier_config = hardcoded_config_20_verifier();

    // Use original prover directly
    let proof = prover(prover_config, poly);
    let verification_result = verifier(verifier_config, proof);

    println!("Verification result: {}", verification_result);
    assert!(verification_result, "Verification should pass");
}
