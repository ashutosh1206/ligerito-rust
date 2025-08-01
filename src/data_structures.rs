use crate::binary_field::BinaryField;
use crate::reed_solomon::ReedSolomonEncoding;
// use crate::
use std::ops::{Add, Mul};

pub struct ProverConfig<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub recursive_steps: usize,
    pub initial_dims: (usize, usize),
    pub dims: Vec<(usize, usize)>,
    pub initial_k: usize,
    pub ks: Vec<usize>,
    pub initial_reed_solomon: ReedSolomonEncoding<F>,
    pub reed_solomon_codes: Vec<ReedSolomonEncoding<F>>,
}

pub struct VerifierConfig {
    pub recursive_steps: usize,
    pub initial_dim: usize,
    pub log_dims: Vec<usize>,
    pub initial_k: usize,
    pub ks: Vec<usize>,
}

pub struct RecursiveLigeroWitness<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
{
    pub flat_mat: Vec<F>,
    pub tree: Vec<Vec<u8>>,
    pub num_rows: usize,
    pub num_cols: usize,
}

pub struct RecursiveLigeroCommitment {
    pub root: Vec<u8>,
}

pub struct RecursiveLigeroProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    // should this be a flat array that we constructured in src/ligero.rs? Probably not
    pub opened_rows: Vec<Vec<F>>,
    pub merkle_proof: Vec<Vec<u8>>,
}

pub struct FinalLigeroProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub yr: Vec<F>,
    pub opened_rows: Vec<Vec<F>>,
    pub merkle_proof: Vec<Vec<u8>>,
}

pub struct SumcheckTranscript<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub tr: Vec<(F, F, F)>,
}

pub struct LigeritoProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub initial_ligero_cm: RecursiveLigeroCommitment,
    pub initial_ligero_proof: Option<RecursiveLigeroProof<F>>,
    pub recursive_commitments: Vec<RecursiveLigeroCommitment>,
    pub recursive_proofs: Vec<RecursiveLigeroProof<F>>,
    pub final_ligero_proof: Option<FinalLigeroProof<F>>,
    pub sumcheck_transcript: Option<SumcheckTranscript<F>>,
}

pub struct FinalizedLigeritoProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub initial_ligero_cm: RecursiveLigeroCommitment,
    pub initial_ligero_proof: RecursiveLigeroProof<F>,
    pub recursive_commitments: Vec<RecursiveLigeroCommitment>,
    pub recursive_proofs: Vec<RecursiveLigeroProof<F>>,
    pub final_ligero_proof: FinalLigeroProof<F>,
    pub sumcheck_transcript: SumcheckTranscript<F>,
}
