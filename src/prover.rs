use crate::data_structures::{ProverConfig, RecursiveLigeroCommitment};
use crate::emulated_fs::FS;
use crate::ligero::ligero_commit;
use crate::{LigeritoProof, binary_field::BinaryField};
use std::ops::{Add, Mul};

pub fn prover<F>(config: ProverConfig<F>, poly: Vec<F>)
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut fs = FS::new(1234);
    let S = 148;

    let wtns_0 = ligero_commit(
        poly,
        config.initial_dims.0,
        config.initial_dims.1,
        config.initial_reed_solomon,
    );
    let proof: LigeritoProof<F> = LigeritoProof {
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

    // let partial_evals_0 =
}
