use crate::binary_fft::eval_sk_at_vks;
use crate::binary_field::BinaryField;
use crate::data_structures::{
    FinalLigeroProof, FinalizedLigeritoProof, LigeritoProof, ProverConfig,
    RecursiveLigeroCommitment, RecursiveLigeroProof, SumcheckTranscript,
};
use crate::emulated_fs::FS;
use crate::ligerito::SumcheckProverInstance;
use crate::ligero::{extract_row, ligero_commit};
use crate::merkle_tree::prove;
use crate::multilinear_poly::MultiLinearPoly;
use crate::sumcheck_polys::induce_sumcheck_poly;
use std::ops::{Add, Mul};

pub fn prover<F>(config: ProverConfig<F>, poly: Vec<F>) -> FinalizedLigeritoProof<F>
where
    F: Copy
        + BinaryField
        + Add<Output = F>
        + Mul<Output = F>
        + PartialEq
        + std::fmt::Debug
        + Send
        + Sync,
    F::ValueType: TryFrom<usize> + TryFrom<u128> + Copy + std::fmt::Debug,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
    <F::ValueType as TryFrom<u128>>::Error: std::fmt::Debug,
{
    let mut fs = FS::new(1234);
    let s = 148;

    let wtns_0 = ligero_commit(
        &poly,
        config.initial_dims.0,
        config.initial_dims.1,
        &config.initial_reed_solomon,
    );
    let mut proof: LigeritoProof<F> = LigeritoProof {
        initial_ligero_cm: RecursiveLigeroCommitment {
            root: wtns_0.tree.last().unwrap().clone(),
        },
        initial_ligero_proof: None,
        recursive_commitments: Vec::with_capacity(config.recursive_steps),
        recursive_proofs: Vec::with_capacity(config.recursive_steps),
        final_ligero_proof: None,
        sumcheck_transcript: None,
    };
    fs.absorb(&proof.initial_ligero_cm.root);

    let partial_evals_0: Vec<F> = (0..config.initial_k).map(|_| fs.get_field()).collect();

    let mut f = MultiLinearPoly::new(poly);
    f = f.partial_eval(&partial_evals_0);

    let wtns_1 = ligero_commit(
        f.evals(),
        config.dims[0].0,
        config.dims[0].1,
        &config.reed_solomon_codes[0],
    );
    let cm_1 = RecursiveLigeroCommitment {
        root: wtns_1.tree.last().unwrap().clone(),
    };
    fs.absorb(&cm_1.root);
    proof.recursive_commitments.push(cm_1);

    let mut queries = fs.get_distinct_queries(wtns_0.num_rows, s);
    let alpha: F = fs.get_field();

    let mut sks_vks = eval_sk_at_vks::<F>(1 << f.num_vars());

    let mut opened_rows: Vec<Vec<F>> = queries
        .iter()
        .map(|&q| extract_row(&wtns_0.flat_mat, q, wtns_0.num_rows, wtns_0.num_cols))
        .collect();
    let mut mtree_proof = prove(&wtns_0.tree, &queries);

    let (mut basis_poly, mut enforced_sum) = induce_sumcheck_poly(
        f.num_vars(),
        &sks_vks,
        &opened_rows,
        &partial_evals_0,
        &queries,
        alpha,
    );
    proof.initial_ligero_proof = Some(RecursiveLigeroProof {
        opened_rows,
        merkle_proof: mtree_proof,
    });

    let (mut sumcheck_prover, s1) =
        SumcheckProverInstance::new(f, MultiLinearPoly::new(basis_poly), enforced_sum);
    fs.absorb(&s1.to_bytes());

    let mut wtns_prev = wtns_1;
    for i in 0..config.recursive_steps {
        let mut rs: Vec<F> = Vec::with_capacity(config.ks[i]);
        for _ in 0..config.ks[i] {
            let ri: F = fs.get_field();
            fs.absorb(&sumcheck_prover.fold(ri).to_bytes());
            rs.push(ri);
        }

        if i == config.recursive_steps - 1 {
            let evals_u8: Vec<u8> = sumcheck_prover
                .f
                .evals()
                .iter()
                .flat_map(|elem| elem.to_bytes())
                .collect();
            fs.absorb(&evals_u8);

            queries = fs.get_distinct_queries(wtns_prev.num_rows, s);
            opened_rows = queries
                .iter()
                .map(|&q| {
                    extract_row(
                        &wtns_prev.flat_mat,
                        q,
                        wtns_prev.num_rows,
                        wtns_prev.num_cols,
                    )
                })
                .collect();
            mtree_proof = prove(&wtns_prev.tree, &queries);
            proof.final_ligero_proof = Some(FinalLigeroProof {
                yr: sumcheck_prover.f.evals().clone(),
                opened_rows,
                merkle_proof: mtree_proof,
            });
            proof.sumcheck_transcript = Some(SumcheckTranscript {
                tr: sumcheck_prover.transcript.clone(), // FIXME: see if passing reference is better here
            });

            return FinalizedLigeritoProof {
                initial_ligero_cm: proof.initial_ligero_cm,
                recursive_commitments: proof.recursive_commitments,
                initial_ligero_proof: proof.initial_ligero_proof.unwrap(), // Handle case where this is None
                recursive_proofs: proof.recursive_proofs,
                final_ligero_proof: proof.final_ligero_proof.unwrap(), // Handle case where this is None
                sumcheck_transcript: proof.sumcheck_transcript.unwrap(), // Handle case where this is None
            };
        }

        let wtns_next = ligero_commit(
            sumcheck_prover.f.evals(),
            config.dims[i + 1].0,
            config.dims[i + 1].1,
            &config.reed_solomon_codes[i + 1],
        );
        let cm_next = RecursiveLigeroCommitment {
            root: wtns_next.tree.last().unwrap().clone(),
        };
        fs.absorb(&cm_next.root);
        proof.recursive_commitments.push(cm_next);

        queries = fs.get_distinct_queries(wtns_prev.num_rows, s);
        let alpha: F = fs.get_field();

        sks_vks = eval_sk_at_vks::<F>(1 << sumcheck_prover.f.num_vars());
        opened_rows = queries
            .iter()
            .map(|&q| {
                extract_row(
                    &wtns_prev.flat_mat,
                    q,
                    wtns_prev.num_rows,
                    wtns_prev.num_cols,
                )
            })
            .collect();
        mtree_proof = prove(&wtns_prev.tree, &queries);

        (basis_poly, enforced_sum) = induce_sumcheck_poly(
            sumcheck_prover.f.num_vars(),
            &sks_vks,
            &opened_rows,
            &rs,
            &queries,
            alpha,
        );
        proof.recursive_proofs.push(RecursiveLigeroProof {
            opened_rows,
            merkle_proof: mtree_proof,
        });

        fs.absorb(
            &sumcheck_prover
                .introduce_new(MultiLinearPoly::new(basis_poly), enforced_sum)
                .to_bytes(),
        );

        sumcheck_prover.glue(fs.get_field());

        wtns_prev = wtns_next;
    }

    panic!("Prover completed without entering final step - check config.recursive_steps");
}
