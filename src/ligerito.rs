use crate::binary_field::BinaryField;
use crate::multilinear_poly::{MultiLinearPoly, eval_013_product};
use crate::sumcheck::QuadraticEvals;
use crate::{QuadraticPoly, fold_quadratic, quadratic_from_evals};
use std::ops::{Add, Mul};

pub struct SumcheckProverInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    f: MultiLinearPoly<F>,
    basis_polys: Vec<MultiLinearPoly<F>>,
    sum: F,
    transcript: Vec<(F, F, F)>,
    to_be_glued: Option<MultiLinearPoly<F>>,
}

impl<F> SumcheckProverInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(f: MultiLinearPoly<F>, b1: MultiLinearPoly<F>, h1: F) -> (Self, QuadraticEvals<F>) {
        let mut transcript: Vec<(F, F, F)> = vec![];
        let (s0, s1, s2) = eval_013_product(&f, &b1);
        assert_eq!(s0 + s1, h1);
        transcript.push((s0, s1, s2));
        let quad_evals = QuadraticEvals::new(s0, s1, s2);

        (
            Self {
                f,
                basis_polys: vec![b1],
                sum: h1,
                transcript,
                to_be_glued: None,
            },
            quad_evals,
        )
    }

    pub fn introduce_new(&mut self, bi: MultiLinearPoly<F>, h: F) -> QuadraticEvals<F> {
        let (s0, s1, s2) = eval_013_product(&self.f, &bi);
        assert_eq!(s0 + s1, h);
        self.transcript.push((s0, s1, s2));
        self.to_be_glued = Some(bi);

        QuadraticEvals::new(s0, s1, s2)
    }

    pub fn glue(&mut self, alpha: F) {
        assert!(self.to_be_glued.is_some());
        let mut poly = self.to_be_glued.take().unwrap();
        poly.scale_evals(alpha);
        self.basis_polys.push(poly);
        // self.to_be_glued.take() already leaves None in its place
        // So this is redundant and a relic of Julia port
        self.to_be_glued = None;
    }

    pub fn eval_01x_product(&self, x: Option<F>) -> (F, F, F) {
        let x = x.unwrap_or_else(|| F::new(F::ValueType::try_from(3).unwrap()));

        let val_at_0 = self.f.partial_eval_at_0();
        let val_at_1 = self.f.partial_eval_at_1();
        let f0 = val_at_0.evals().clone();
        let f1 = val_at_1.evals().clone();
        let f2 = evals_at_x(&f0, &f1, x);

        let bs_val_at_0 = self.basis_polys[0].partial_eval_at_0();
        let bs_val_at_1 = self.basis_polys[0].partial_eval_at_1();
        let mut b0s_sum = bs_val_at_0.evals().clone();
        let mut b1s_sum = bs_val_at_1.evals().clone();
        let mut b2s_sum = evals_at_x(&b0s_sum, &b1s_sum, x);

        for i in 1..self.basis_polys.len() {
            let bi_val_at_0 = self.basis_polys[i].partial_eval_at_0();
            let bi_val_at_1 = self.basis_polys[i].partial_eval_at_1();
            let b0i = bi_val_at_0.evals().clone();
            let b1i = bi_val_at_1.evals().clone();
            let b2i = evals_at_x(&b0i, &b1i, x);

            for i in 0..b0s_sum.len() {
                b0s_sum[i] = b0s_sum[i] + b0i[i];
            }
            for i in 0..b1s_sum.len() {
                b1s_sum[i] = b1s_sum[i] + b1i[i];
            }
            for i in 0..b2s_sum.len() {
                b2s_sum[i] = b2s_sum[i] + b2i[i];
            }
        }
        let s0 = f0
            .iter()
            .zip(b0s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        let s1 = f1
            .iter()
            .zip(b1s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        let s2 = f2
            .iter()
            .zip(b2s_sum.iter())
            .map(|(&f_val, &b_val)| f_val * b_val)
            .fold(F::zero(), |acc, x| acc + x);
        (s0, s1, s2)
    }

    pub fn fold(&mut self, r: F) -> QuadraticEvals<F> {
        self.f = self.f.partial_eval(&vec![r]);
        for i in 0..self.basis_polys.len() {
            self.basis_polys[i] = self.basis_polys[i].partial_eval(&vec![r]);
        }
        let (s0, s1, s2) = self.eval_01x_product(None);
        self.transcript.push((s0, s1, s2));
        QuadraticEvals::new(s0, s1, s2)
    }
}

fn evals_at_x<F>(evals_at_0: &Vec<F>, evals_at_1: &Vec<F>, alpha: F) -> Vec<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    evals_at_0
        .iter()
        .zip(evals_at_1.iter())
        .map(|(&val0, &val1)| alpha * val1 + (F::one() + alpha) * val0)
        .collect()
}

pub struct SumcheckVerifierInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    basis_polys: Vec<MultiLinearPoly<F>>,
    separation_challenges: Vec<F>,
    sum: F,
    transcript: Vec<(F, F, F)>,
    ris: Vec<F>,
    tr_reader: usize,
    running_poly: Option<QuadraticPoly<F>>,
    to_glue: Option<QuadraticPoly<F>>,
}

impl<F> SumcheckVerifierInstance<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F> + PartialEq + std::fmt::Debug,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new(
        b1: MultiLinearPoly<F>,
        h1: F,
        transcript: Vec<(F, F, F)>,
    ) -> (Self, QuadraticEvals<F>) {
        let mut verifier = Self {
            basis_polys: vec![b1],
            separation_challenges: vec![F::one()],
            sum: h1,
            transcript,
            ris: vec![],
            tr_reader: 0,
            running_poly: None,
            to_glue: None,
        };
        let (g0, g1, g2) = verifier.read_tr();
        assert_eq!(g0 + g1, verifier.sum);

        verifier.running_poly = Some(quadratic_from_evals(g0, g1, g2, None));
        (verifier, QuadraticEvals::new(g0, g1, g2))
    }

    pub fn read_tr(&mut self) -> (F, F, F) {
        assert!(self.tr_reader < self.transcript.len());
        let (g0, g1, g2) = self.transcript[self.tr_reader];
        self.tr_reader += 1;
        return (g0, g1, g2);
    }

    pub fn fold(&mut self, r: F) -> QuadraticEvals<F> {
        self.ris.push(r);
        assert!(self.running_poly.is_some());
        self.sum = self.running_poly.take().unwrap().eval_quadratic(r);

        let (g0, g1, g2) = self.read_tr();
        assert_eq!(g0 + g1, self.sum);
        self.running_poly = Some(quadratic_from_evals(g0, g1, g2, None));
        QuadraticEvals::new(g0, g1, g2)
    }

    pub fn introduce_new(&mut self, bi: MultiLinearPoly<F>, h: F) -> QuadraticEvals<F> {
        let (g0, g1, g2) = self.read_tr();
        assert_eq!(g0 + g1, h);

        self.basis_polys.push(bi);
        self.to_glue = Some(quadratic_from_evals(g0, g1, g2, None));

        QuadraticEvals::new(g0, g1, g2)
    }

    pub fn glue(&mut self, alpha: F) {
        assert!(self.running_poly.is_some());
        assert!(self.to_glue.is_some());

        self.separation_challenges.push(alpha);
        self.running_poly = Some(fold_quadratic(
            self.running_poly.take().unwrap(),
            self.to_glue.as_ref().unwrap().clone(),
            alpha,
        ))
    }

    pub fn evaluate_basis_polys(&mut self, r: F) -> F {
        self.ris.push(r);
        let mut b_eval = self.basis_polys[0].partial_eval(&self.ris).evals()[0];

        for i in 1..self.basis_polys.len() {
            let n = self.basis_polys[i].num_vars();
            let eval_pts = self.ris[self.ris.len() - n..].to_vec();
            let bi_eval = self.basis_polys[i].partial_eval(&eval_pts).evals()[0];

            b_eval = b_eval + self.separation_challenges[i] * bi_eval;
        }

        b_eval
    }

    pub fn verify(&mut self, r: F, f_eval: F) -> bool {
        assert!(self.running_poly.is_some());
        let running_poly = self.running_poly.as_ref().unwrap().clone();
        self.sum = running_poly.eval_quadratic(r);
        let basis_evals = self.evaluate_basis_polys(r);
        f_eval * basis_evals == self.sum
    }

    pub fn evaluate_basis_polys_partially(&mut self, r: F, k: usize) -> Vec<F> {
        self.ris.push(r);
        let partial_eval_res = self.basis_polys[0].partial_eval(&self.ris);
        let mut acc = partial_eval_res.evals().to_vec();

        for i in 1..self.basis_polys.len() {
            let n = self.basis_polys[i].num_vars();
            // self.ris.len() >= n - k so `self.ris.len() - (n - k)` will never underflow
            // but `self.ris.len() - n + k` can underflow, causing panic
            let eval_pts = &self.ris[self.ris.len() - (n - k)..];
            let partial_eval_res = self.basis_polys[i].partial_eval(&eval_pts);
            let bi_evals = partial_eval_res.evals();
            let alpha = self.separation_challenges[i];

            assert_eq!(acc.len(), bi_evals.len());
            acc = acc
                .iter()
                .zip(bi_evals.iter())
                .map(|(&acci, &bievali)| acci + alpha * bievali)
                .collect();
        }

        acc
    }

    pub fn verify_partial(&mut self, r: F, f_partial_eval: Vec<F>) -> bool {
        let k = f_partial_eval.len().ilog2() as usize;
        assert!(self.running_poly.is_some());
        self.sum = self.running_poly.as_ref().unwrap().eval_quadratic(r);
        let basis_evals = self.evaluate_basis_polys_partially(r, k);

        assert_eq!(f_partial_eval.len(), basis_evals.len());
        f_partial_eval
            .iter()
            .zip(basis_evals.iter())
            .map(|(&fevali, &bevali)| fevali * bevali)
            .fold(F::zero(), |acc, x| acc + x)
            == self.sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_field::BinaryElem16;
    use crate::binary_field::random;
    use sha2::{Digest, Sha256};
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    // Simple Fiat-Shamir transcript simulator
    struct FSTranscript {
        hasher: Sha256,
        counter: u32,
    }

    impl FSTranscript {
        fn new(seed: u32) -> Self {
            let mut hasher = Sha256::new();
            hasher.update(&seed.to_le_bytes());
            Self { hasher, counter: 0 }
        }

        fn absorb(&mut self, evals: &QuadraticEvals<BinaryElem16>) {
            // In real implementation, you'd properly serialize the field elements
            // For now, we'll use a simple approach
            let mut hasher = DefaultHasher::new();
            evals.hash(&mut hasher);
            let hash = hasher.finish();
            self.hasher.update(&hash.to_le_bytes());
        }

        fn squeeze(&mut self) -> BinaryElem16 {
            self.hasher.update(&self.counter.to_le_bytes());
            self.counter += 1;
            let result = self.hasher.finalize_reset();
            let value = u16::from_le_bytes([result[0], result[1]]);
            self.hasher = Sha256::new();
            BinaryElem16::new(value)
        }
    }

    fn random_poly(k: usize) -> MultiLinearPoly<BinaryElem16> {
        let evals: Vec<BinaryElem16> = (0..(1 << k)).map(|_| random::<BinaryElem16>()).collect();
        MultiLinearPoly::new(evals)
    }

    fn inner_product(
        f: &MultiLinearPoly<BinaryElem16>,
        b: &MultiLinearPoly<BinaryElem16>,
    ) -> BinaryElem16 {
        f.evals()
            .iter()
            .zip(b.evals().iter())
            .map(|(&fi, &bi)| fi * bi)
            .fold(BinaryElem16::zero(), |acc, x| acc + x)
    }

    #[test]
    fn test_ligerito_emulator() {
        let k = 12;
        let glues = vec![2, 5, 9];
        let bs: Vec<MultiLinearPoly<BinaryElem16>> =
            glues.iter().map(|&gi| random_poly(k - gi)).collect();

        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        // Store challenges and hs for verifier
        let mut prover_challenges = Vec::new();
        let mut verifier_challenges = Vec::new();
        let mut hs = Vec::new();

        // === PROVER - Complete execution ===
        let mut fs_prover = FSTranscript::new(1234);
        let (mut prover, s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        fs_prover.absorb(&s1);

        let mut folds = 0;
        let mut gl_idx = 0;

        for _i in 0..(k - 1) {
            let ri = fs_prover.squeeze();
            prover_challenges.push(ri);
            let si = prover.fold(ri);
            fs_prover.absorb(&si);
            folds += 1;

            if gl_idx < glues.len() && folds == glues[gl_idx] {
                let bi = bs[gl_idx].clone();
                let hi = inner_product(&prover.f, &bi);
                hs.push(hi);

                let gl_i = prover.introduce_new(bi, hi);
                fs_prover.absorb(&gl_i);

                let alpha = fs_prover.squeeze();
                prover.glue(alpha);

                gl_idx += 1;
            }
        }

        // === VERIFIER - Use complete transcript ===
        let mut fs_verifier = FSTranscript::new(1234);
        let mut folds = 0;
        let mut gl_idx = 0;
        let (mut verifier, g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());
        fs_verifier.absorb(&g1);

        for i in 0..(k - 1) {
            let ri = fs_verifier.squeeze();
            verifier_challenges.push(ri);
            let gi = verifier.fold(ri);
            fs_verifier.absorb(&gi);
            folds += 1;

            if gl_idx < glues.len() && folds == glues[gl_idx] {
                let bi = bs[gl_idx].clone();
                let hi = hs[gl_idx];

                let gl_i = verifier.introduce_new(bi, hi);
                fs_verifier.absorb(&gl_i);

                let alpha = fs_verifier.squeeze();
                verifier.glue(alpha);

                gl_idx += 1;
            }
        }

        // Final check - emulate oracle access to f
        let final_ri = fs_verifier.squeeze();
        verifier_challenges.push(final_ri);
        let f_eval = f.partial_eval(&verifier_challenges).evals()[0];

        // Perform final verification
        let ok = verifier.verify(final_ri, f_eval);
        assert!(ok, "Ligerito protocol verification should pass");

        // Verify that prover and verifier used same challenges
        assert_eq!(
            prover_challenges,
            verifier_challenges[..prover_challenges.len()]
        );
    }

    #[test]
    fn test_ligerito_emulator_no_glues() {
        // Test without any glues first to isolate the issue
        let k = 5;
        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        // Store challenges for verifier
        let mut prover_challenges = Vec::new();
        let mut verifier_challenges = Vec::new();

        // === PROVER - Complete execution ===
        let mut fs_prover = FSTranscript::new(1234);
        let (mut prover, s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        fs_prover.absorb(&s1);

        for _i in 0..(k - 1) {
            let ri = fs_prover.squeeze();
            prover_challenges.push(ri);
            let si = prover.fold(ri);
            fs_prover.absorb(&si);
        }

        // === VERIFIER - Use complete transcript ===
        let mut fs_verifier = FSTranscript::new(1234);
        let (mut verifier, g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());
        fs_verifier.absorb(&g1);

        for _i in 0..(k - 1) {
            let ri = fs_verifier.squeeze();
            verifier_challenges.push(ri);
            let gi = verifier.fold(ri);
            fs_verifier.absorb(&gi);
        }

        // Final check - emulate oracle access to f
        let final_ri = fs_verifier.squeeze();
        verifier_challenges.push(final_ri);
        let f_eval = f.partial_eval(&verifier_challenges).evals()[0];

        // Perform final verification
        let ok = verifier.verify(final_ri, f_eval);
        assert!(ok, "Ligerito protocol without glues should pass");

        // Verify that prover and verifier used same challenges
        assert_eq!(
            prover_challenges,
            verifier_challenges[..prover_challenges.len()]
        );
    }

    #[test]
    fn test_ligerito_simple() {
        // Start with a smaller test first
        let k = 5;
        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        // === PROVER - Generate complete transcript ===
        let (mut prover, _s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);

        // Generate some deterministic challenges
        let mut challenges = Vec::new();
        for i in 0..(k - 1) {
            let ri = BinaryElem16::new((i + 1) as u16);
            challenges.push(ri);
            prover.fold(ri);
        }

        // === VERIFIER - Use complete transcript ===
        let (mut verifier, _g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());

        for &ri in &challenges {
            verifier.fold(ri);
        }

        // Final check - emulate oracle access to f
        let final_r = BinaryElem16::new(42);
        let eval_points: Vec<BinaryElem16> = challenges
            .iter()
            .chain(std::iter::once(&final_r))
            .cloned()
            .collect();
        let f_eval = f.partial_eval(&eval_points).evals()[0];

        // Perform final verification
        let ok = verifier.verify(final_r, f_eval);
        assert!(ok, "Simple ligerito protocol verification should pass");
    }

    #[test]
    fn test_basic_sumcheck_debug() {
        // Minimal test to debug the assertion issue
        let k = 3;
        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        // === PROVER ===
        let (mut prover, s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        println!("Initial s1: {:?}", s1);

        // Just do one fold
        let r1 = BinaryElem16::new(1);
        let s2 = prover.fold(r1);
        println!("After first fold s2: {:?}", s2);
        println!("Prover transcript: {:?}", prover.transcript);

        // === VERIFIER ===
        let (mut verifier, g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());
        println!("Initial g1: {:?}", g1);
        println!("Verifier initial sum: {:?}", verifier.sum);
        println!("Verifier initial tr_reader: {:?}", verifier.tr_reader);
        println!(
            "Verifier initial transcript length: {:?}",
            verifier.transcript.len()
        );

        let g2 = verifier.fold(r1);
        println!("After first fold g2: {:?}", g2);

        // Try a second fold
        let r2 = BinaryElem16::new(2);
        let s3 = prover.fold(r2);
        println!("After second fold s3: {:?}", s3);

        println!("About to do second verifier fold...");
        println!("Verifier sum before fold: {:?}", verifier.sum);
        println!("Verifier tr_reader: {:?}", verifier.tr_reader);
        println!("Transcript length: {:?}", verifier.transcript.len());

        // The issue is here - let me check if the transcript has enough entries
        if verifier.tr_reader >= verifier.transcript.len() {
            println!(
                "ERROR: tr_reader {} >= transcript.len() {}",
                verifier.tr_reader,
                verifier.transcript.len()
            );
            return;
        }

        let g3 = verifier.fold(r2);
        println!("After second fold g3: {:?}", g3);
    }

    #[test]
    fn test_ligerito_deterministic() {
        // Test with deterministic challenges instead of Fiat-Shamir
        let k = 3;
        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        // Use deterministic challenges
        let challenges: Vec<BinaryElem16> = (1..(k as u16)).map(BinaryElem16::new).collect();

        // === PROVER - Complete execution ===
        let (mut prover, _s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);

        for &ri in &challenges {
            prover.fold(ri);
        }

        // === VERIFIER - Use complete transcript ===
        let (mut verifier, _g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());

        for &ri in &challenges {
            verifier.fold(ri);
        }

        // Final check - emulate oracle access to f
        let final_r = BinaryElem16::new(42);
        let mut all_challenges = challenges.clone();
        all_challenges.push(final_r);
        let f_eval = f.partial_eval(&all_challenges).evals()[0];

        // Perform final verification
        let ok = verifier.verify(final_r, f_eval);
        assert!(
            ok,
            "Ligerito protocol with deterministic challenges should pass"
        );
    }

    #[test]
    fn test_ligerito_minimal() {
        // Super minimal test with k=2 to trace through exactly
        let k = 2;
        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_product(&f, &b1);

        println!("=== MINIMAL TEST k=2 ===");
        println!("f.evals = {:?}", f.evals());
        println!("b1.evals = {:?}", b1.evals());
        println!("h = {:?}", h);

        // === PROVER - Complete execution ===
        let (mut prover, s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        println!("Initial s1: {:?}", s1);

        // Only one fold for k=2
        let r1 = BinaryElem16::new(1);
        let s2 = prover.fold(r1);
        println!("After fold s2: {:?}", s2);

        // === VERIFIER - Use complete transcript ===
        let (mut verifier, g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());
        println!("Initial g1: {:?}", g1);

        let g2 = verifier.fold(r1);
        println!("After fold g2: {:?}", g2);

        // Final check
        let final_r = BinaryElem16::new(42);
        let f_eval = f.partial_eval(&vec![r1, final_r]).evals()[0];
        println!("f_eval = {:?}", f_eval);

        let ok = verifier.verify(final_r, f_eval);
        assert!(ok, "Minimal test should pass");
    }

    #[test]
    fn test_quadratic_construction() {
        // Test the quadratic_from_evals function independently
        let at0 = BinaryElem16::new(100);
        let at1 = BinaryElem16::new(200);
        let atx = BinaryElem16::new(150);

        let quad = quadratic_from_evals(at0, at1, atx, None);

        println!("Testing quadratic construction:");
        println!("at0={:?}, at1={:?}, atx={:?}", at0, at1, atx);

        // Verify the quadratic evaluates correctly at the three points
        let eval_0 = quad.eval_quadratic(BinaryElem16::zero());
        let eval_1 = quad.eval_quadratic(BinaryElem16::one());
        let eval_3 = quad.eval_quadratic(BinaryElem16::new(3));

        println!(
            "quad(0)={:?}, quad(1)={:?}, quad(3)={:?}",
            eval_0, eval_1, eval_3
        );

        assert_eq!(eval_0, at0);
        assert_eq!(eval_1, at1);
        assert_eq!(eval_3, atx);

        // Also test at another point
        let eval_2 = quad.eval_quadratic(BinaryElem16::new(2));
        println!("quad(2)={:?}", eval_2);
    }

    // Helper function to match Julia's inner function
    fn inner_with_rs(
        f: &MultiLinearPoly<BinaryElem16>,
        b: &MultiLinearPoly<BinaryElem16>,
        rs: &[BinaryElem16],
    ) -> BinaryElem16 {
        let n1 = f.num_vars();
        let n2 = b.num_vars();

        // Julia: eval_pts = rs[1:n1 - n2] (1-based indexing)
        // Rust: eval_pts = rs[0..(n1 - n2)] (0-based indexing)
        // Handle case where n1 - n2 might be 0 (empty slice) or larger than rs.len()
        let eval_len = if n1 >= n2 { n1 - n2 } else { 0 };
        let eval_len = eval_len.min(rs.len()); // Don't exceed rs length
        let eval_pts = &rs[0..eval_len];

        let fp = if eval_pts.is_empty() {
            // If no evaluation points, return f unchanged
            f.clone()
        } else {
            f.partial_eval(&eval_pts)
        };
        assert_eq!(fp.num_vars(), b.num_vars());

        fp.evals()
            .iter()
            .zip(b.evals().iter())
            .map(|(&fi, &bi)| fi * bi)
            .fold(BinaryElem16::zero(), |acc, x| acc + x)
    }

    #[test]
    fn test_ligerito_partial_emulator() {
        // Port of the Julia ligerito_partial_emulator.jl test
        let k = 12;
        let rs: Vec<BinaryElem16> = (0..(k - 4)).map(|_| random::<BinaryElem16>()).collect(); // rs = rand(BinaryElem16, k - 4)
        let glues = vec![2, 5];
        let bs: Vec<MultiLinearPoly<BinaryElem16>> =
            glues.iter().map(|&gi| random_poly(k - gi)).collect(); // bs = [random_poly(BinaryElem16, k - gi) for gi in glues]

        let separation_challenges: Vec<BinaryElem16> =
            (0..glues.len()).map(|_| random::<BinaryElem16>()).collect(); // separation_challenges = rand(BinaryElem16, length(glues))

        let f = random_poly(k);
        let b1 = random_poly(k);
        let h = inner_with_rs(&f, &b1, &rs); // h = inner(f, b1, rs)

        // Store hs for verifier
        let mut hs = Vec::new();

        // === PROVER ===
        let (mut prover, _s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        let mut folds = 0;
        let mut gl_idx = 0;

        // prover folds rs.len() - 1 times in total (Julia: for i in 1:(length(rs) - 1))
        for i in 0..(rs.len() - 1) {
            prover.fold(rs[i]);
            folds += 1;

            if gl_idx < glues.len() && folds == glues[gl_idx] {
                let bi = bs[gl_idx].clone();
                let alpha = separation_challenges[gl_idx];

                // Julia: hi = inner_from_running(prover, bi)
                let hi = inner_product(&prover.f, &bi);
                hs.push(hi);

                prover.introduce_new(bi, hi);
                prover.glue(alpha);

                gl_idx += 1;
            }
        }

        // === VERIFIER ===
        let mut folds = 0;
        let mut gl_idx = 0;
        let (mut verifier, _g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());

        // verifier folds rs.len() - 1 times before the final check
        for i in 0..(rs.len() - 1) {
            verifier.fold(rs[i]);
            folds += 1;

            if gl_idx < glues.len() && folds == glues[gl_idx] {
                let bi = bs[gl_idx].clone();
                let hi = hs[gl_idx];
                let alpha = separation_challenges[gl_idx];

                verifier.introduce_new(bi, hi);
                verifier.glue(alpha);

                gl_idx += 1;
            }
        }

        // Final check - emulate oracle access to f
        // Julia: f_partial_eval = partial_eval(f, rs).evals
        let f_partial_eval = f.partial_eval(&rs).evals().clone();

        // Julia: ok = verify_partial(verifier, rs[end], f_partial_eval)
        let ok = verifier.verify_partial(rs[rs.len() - 1], f_partial_eval);
        assert!(ok, "Ligerito partial emulator verification should pass");
    }

    #[test]
    fn test_specific_failing_case() {
        // Use the exact values from the failing test to debug
        let k = 3;

        // Create specific polynomials with known values to trace the issue
        let f_evals = vec![
            BinaryElem16::new(1),
            BinaryElem16::new(2),
            BinaryElem16::new(3),
            BinaryElem16::new(4),
            BinaryElem16::new(5),
            BinaryElem16::new(6),
            BinaryElem16::new(7),
            BinaryElem16::new(8),
        ];
        let f = MultiLinearPoly::new(f_evals);

        let b1_evals = vec![
            BinaryElem16::new(10),
            BinaryElem16::new(20),
            BinaryElem16::new(30),
            BinaryElem16::new(40),
            BinaryElem16::new(50),
            BinaryElem16::new(60),
            BinaryElem16::new(70),
            BinaryElem16::new(80),
        ];
        let b1 = MultiLinearPoly::new(b1_evals);
        let h = inner_product(&f, &b1);

        println!("=== SPECIFIC FAILING CASE ===");
        println!("f.evals = {:?}", f.evals());
        println!("b1.evals = {:?}", b1.evals());
        println!("h = {:?}", h);

        // === PROVER ===
        let (mut prover, s1) = SumcheckProverInstance::new(f.clone(), b1.clone(), h);
        println!("Initial s1: {:?}", s1);

        let r1 = BinaryElem16::new(1);
        let s2 = prover.fold(r1);
        println!("After fold 1 s2: {:?}", s2);

        let r2 = BinaryElem16::new(2);
        let s3 = prover.fold(r2);
        println!("After fold 2 s3: {:?}", s3);

        // === VERIFIER ===
        let (mut verifier, g1) = SumcheckVerifierInstance::new(b1, h, prover.transcript.clone());
        println!("Initial g1: {:?}", g1);

        let g2 = verifier.fold(r1);
        println!("After fold 1 g2: {:?}", g2);

        // This should pass but might not
        let g3 = verifier.fold(r2);
        println!("After fold 2 g3: {:?}", g3);
    }
}
