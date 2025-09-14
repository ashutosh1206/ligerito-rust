use crate::binary_field::BinaryField;
use crate::reed_solomon::ReedSolomonEncoding;
// use crate::
use std::ops::{Add, Mul};

pub trait ProofSize {
    fn proof_size(&self) -> usize;
}

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

impl<F> ProofSize for RecursiveLigeroProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    fn proof_size(&self) -> usize {
        let opened_rows_size = self
            .opened_rows
            .iter()
            .map(|row| row.len() * std::mem::size_of::<F>())
            .sum::<usize>();
        let merkle_proof_size = self
            .merkle_proof
            .iter()
            .map(|proof| proof.len())
            .sum::<usize>();
        opened_rows_size + merkle_proof_size
    }
}

impl<F> ProofSize for FinalLigeroProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    fn proof_size(&self) -> usize {
        let yr_size = self.yr.len() * std::mem::size_of::<F>();
        let opened_rows_size = self
            .opened_rows
            .iter()
            .map(|row| row.len() * std::mem::size_of::<F>())
            .sum::<usize>();
        let merkle_proof_size = self
            .merkle_proof
            .iter()
            .map(|proof| proof.len())
            .sum::<usize>();
        yr_size + opened_rows_size + merkle_proof_size
    }
}

impl<F> ProofSize for SumcheckTranscript<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    fn proof_size(&self) -> usize {
        self.tr.len() * std::mem::size_of::<(F, F, F)>()
    }
}

impl ProofSize for RecursiveLigeroCommitment {
    fn proof_size(&self) -> usize {
        self.root.len()
    }
}

impl<F> ProofSize for FinalizedLigeritoProof<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    fn proof_size(&self) -> usize {
        let initial_cm_size = self.initial_ligero_cm.proof_size();
        let initial_proof_size = self.initial_ligero_proof.proof_size();
        let recursive_cms_size = self
            .recursive_commitments
            .iter()
            .map(|cm| cm.proof_size())
            .sum::<usize>();
        let recursive_proofs_size = self
            .recursive_proofs
            .iter()
            .map(|proof| proof.proof_size())
            .sum::<usize>();
        let final_proof_size = self.final_ligero_proof.proof_size();
        let sumcheck_size = self.sumcheck_transcript.proof_size();

        initial_cm_size
            + initial_proof_size
            + recursive_cms_size
            + recursive_proofs_size
            + final_proof_size
            + sumcheck_size
    }
}
