use crate::binary_fft::eval_sk_at_vks;
use crate::binary_field::BinaryField;
use crate::data_structures::{LigeritoProof, ProverConfig, RecursiveLigeroCommitment};
use crate::emulated_fs::FS;
use crate::ligero::ligero_commit;
use crate::multilinear_poly::MultiLinearPoly;
use std::ops::{Add, Mul};

pub fn prover<F>(config: ProverConfig<F>, poly: Vec<F>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize> + From<u128>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut fs = FS::new(1234);
    let S = 148;

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
    f = f.partial_eval(partial_evals_0);

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

    let queries = fs.get_distinct_queries(wtns_0.num_rows, S);
    let alpha: F = fs.get_field();

    let sks_vks = eval_sk_at_vks::<F>(f.num_vars());

    // let opened_rows =
}
