use crate::binary_field::BinaryField;
use crate::data_structures::{FinalizedLigeritoProof, VerifierConfig};
use crate::emulated_fs::FS;
use crate::ligerito::SumcheckVerifierInstance;
use crate::merkle_tree::verify;
use crate::sumcheck_polys::induce_sumcheck_poly;
use crate::{MultiLinearPoly, eval_sk_at_vks, ligero_verify};
use std::ops::{Add, Mul};

pub fn verifier<F>(config: VerifierConfig, proof: FinalizedLigeritoProof<F>) -> bool
where
    F: Copy
        + BinaryField
        + Add<Output = F>
        + Mul<Output = F>
        + PartialEq
        + std::fmt::Debug
        + Send
        + Sync,
    F::ValueType: TryFrom<usize> + TryFrom<u128> + Copy,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
    <F::ValueType as TryFrom<u128>>::Error: std::fmt::Debug,
{
    let mut fs = FS::new(1234);
    let s = 148;
    let log_inv_rate = 2;

    fs.absorb(&proof.initial_ligero_cm.root);

    let partial_evals_0: Vec<F> = (0..config.initial_k).map(|_| fs.get_field()).collect();
    fs.absorb(&proof.recursive_commitments[0].root);

    let mut depth = config.initial_dim + log_inv_rate;
    let mut queries = fs.get_distinct_queries(1 << depth, s);

    let mut opened_rows_bytes: Vec<Vec<u8>> = proof
        .initial_ligero_proof
        .opened_rows
        .iter()
        .map(|row| row.iter().flat_map(|elem| elem.to_bytes()).collect())
        .collect();

    let is_valid = verify(
        &proof.initial_ligero_cm.root,
        &proof.initial_ligero_proof.merkle_proof,
        depth,
        &opened_rows_bytes,
        &queries,
    );
    assert_eq!(is_valid, true);

    let alpha: F = fs.get_field();
    let mut sks_vks = eval_sk_at_vks::<F>(1 << config.initial_dim);
    let (mut basis_poly, mut enforced_sum) = induce_sumcheck_poly(
        config.initial_dim,
        &sks_vks,
        &proof.initial_ligero_proof.opened_rows,
        &partial_evals_0,
        &queries,
        alpha,
    );

    let (mut sumcheck_verifier, g1) = SumcheckVerifierInstance::new(
        MultiLinearPoly::new(basis_poly),
        enforced_sum,
        proof.sumcheck_transcript.tr,
    );
    fs.absorb(&g1.to_bytes());

    for i in 0..config.recursive_steps {
        let mut rs: Vec<F> = Vec::with_capacity(config.ks[i]);
        for _ in 0..config.ks[i] {
            let ri: F = fs.get_field();
            fs.absorb(&sumcheck_verifier.fold(ri).to_bytes());
            rs.push(ri)
        }

        let root = &proof.recursive_commitments[i].root;
        if i == config.recursive_steps - 1 {
            let yr_u8: Vec<u8> = proof
                .final_ligero_proof
                .yr
                .iter()
                .flat_map(|elem| elem.to_bytes())
                .collect();
            fs.absorb(&yr_u8);

            depth = config.log_dims[i] + log_inv_rate;
            queries = fs.get_distinct_queries(1 << depth, s);
            opened_rows_bytes = proof
                .final_ligero_proof
                .opened_rows
                .iter()
                .map(|row| row.iter().flat_map(|elem| elem.to_bytes()).collect())
                .collect();

            assert_eq!(
                verify(
                    root,
                    &proof.final_ligero_proof.merkle_proof,
                    depth,
                    &opened_rows_bytes,
                    &queries,
                ),
                true
            );

            ligero_verify(
                &queries,
                &proof.final_ligero_proof.opened_rows,
                &proof.final_ligero_proof.yr,
                &rs,
            );

            let final_r: F = fs.get_field();
            let f = MultiLinearPoly::new(proof.final_ligero_proof.yr.clone());
            let partial_eval = f.partial_eval(&vec![final_r]);
            let f_eval = partial_eval.evals();

            return sumcheck_verifier.verify_partial(final_r, f_eval);
        }

        fs.absorb(&proof.recursive_commitments[i + 1].root);

        depth = config.log_dims[i] + log_inv_rate;
        let liger_proof = &proof.recursive_proofs[i];
        let queries = fs.get_distinct_queries(1 << depth, s);

        let liger_opened_rows_bytes: Vec<Vec<u8>> = liger_proof
            .opened_rows
            .iter()
            .map(|row| row.iter().flat_map(|elem| elem.to_bytes()).collect())
            .collect();

        assert_eq!(
            verify(
                root,
                &liger_proof.merkle_proof,
                depth,
                &liger_opened_rows_bytes,
                &queries,
            ),
            true
        );

        let alpha: F = fs.get_field();
        sks_vks = eval_sk_at_vks(1 << config.log_dims[i]);
        (basis_poly, enforced_sum) = induce_sumcheck_poly(
            config.log_dims[i],
            &sks_vks,
            &liger_proof.opened_rows,
            &rs,
            &queries,
            alpha,
        );
        let gl_i = sumcheck_verifier.introduce_new(MultiLinearPoly::new(basis_poly), enforced_sum);
        fs.absorb(&gl_i.to_bytes());

        let beta: F = fs.get_field();
        sumcheck_verifier.glue(beta);
    }

    panic!("Verifier completed without entering final step - check config.recursive_steps");
}
