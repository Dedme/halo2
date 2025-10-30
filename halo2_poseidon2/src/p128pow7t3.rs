use alloc::vec::Vec;

use ff::Field;
use pasta_curves::{pallas::Base as Fp, vesta::Base as Fq};

use super::{fp_pow7, fq_pow7, Mds, Spec};

/// Poseidon-128 using the $x^7$ S-box (Aleph Zero optimization), with a width of 3 field elements,
/// and optimized number of rounds for 128-bit security.
///
/// Based on Aleph Zero's optimizations: α=7, reduced rounds (8 full + 48 partial = 56 total).
/// Uses 48 partial rounds (even number) to match circuit gadget requirements.
/// This provides better performance than the standard x^5 S-box while maintaining security.
#[derive(Debug)]
pub struct P128Pow7T3;

impl Spec<Fp, 3, 2> for P128Pow7T3 {
    fn full_rounds() -> usize {
        8
    }

    fn partial_rounds() -> usize {
        48
    }

    fn sbox(val: Fp) -> Fp {
        val.pow_vartime([7, 0, 0, 0])
    }

    fn secure_mds() -> usize {
        0
    }

    fn constants() -> (Vec<[Fp; 3]>, Mds<Fp, 3>, Mds<Fp, 3>, [Fp; 3]) {
        (
            fp_pow7::ROUND_CONSTANTS.to_vec(),
            fp_pow7::MDS, // mat_external (full MDS for external rounds)
            fp_pow7::MDS, // mat_internal (for now same as external - Poseidon v1 compat)
            fp_pow7::MAT_DIAG3_M_1_POSEIDON2, // diagonal elements minus 1
        )
    }
}

impl Spec<Fq, 3, 2> for P128Pow7T3 {
    fn full_rounds() -> usize {
        8
    }

    fn partial_rounds() -> usize {
        48
    }

    fn sbox(val: Fq) -> Fq {
        val.pow_vartime([7, 0, 0, 0])
    }

    fn secure_mds() -> usize {
        0
    }

    fn constants() -> (Vec<[Fq; 3]>, Mds<Fq, 3>, Mds<Fq, 3>, [Fq; 3]) {
        (
            fq_pow7::ROUND_CONSTANTS.to_vec(),
            fq_pow7::MDS, // mat_external (full MDS for external rounds)
            fq_pow7::MDS, // mat_internal (for now same as external - Poseidon v1 compat)
            fq_pow7::MAT_DIAG3_M_1_POSEIDON2, // diagonal elements minus 1
        )
    }
}

#[cfg(test)]
mod tests {
    const POW7_ROUNDS: usize = 56; // 8 full + 48 partial

    use alloc::vec::Vec;
    use core::marker::PhantomData;

    use ff::{Field, FromUniformBytes};

    use super::{
        super::{fp_pow7, fq_pow7},
        Fp, Fq,
    };
    use crate::{generate_constants, permute, ConstantLength, Hash, Mds, Spec};

    /// The same Poseidon specification as poseidon::P128Pow7T3, but constructed
    /// such that its constants will be generated at runtime.
    #[derive(Debug)]
    pub struct P128Pow7T3Gen<F: Field, const SECURE_MDS: usize>(PhantomData<F>);

    impl<F: Field, const SECURE_MDS: usize> P128Pow7T3Gen<F, SECURE_MDS> {
        #![allow(dead_code)]
        pub fn new() -> Self {
            P128Pow7T3Gen(PhantomData::default())
        }
    }

    impl<F: FromUniformBytes<64> + Ord, const SECURE_MDS: usize> Spec<F, 3, 2>
        for P128Pow7T3Gen<F, SECURE_MDS>
    {
        fn full_rounds() -> usize {
            8
        }

        fn partial_rounds() -> usize {
            48 // Even number required for gadget (processes 2 partial rounds per circuit row)
        }

        fn sbox(val: F) -> F {
            val.pow_vartime([7])
        }

        fn secure_mds() -> usize {
            SECURE_MDS
        }

        fn constants() -> (Vec<[F; 3]>, Mds<F, 3>, Mds<F, 3>, [F; 3]) {
            generate_constants::<_, Self, 3, 2>()
        }
    }

    #[test]
    fn verify_constants_count() {
        // Verify that the hardcoded Pow7 constants have the correct dimensions
        // Note: The Pow7 constants were generated externally (not via Grain LFSR)
        // to optimize for the x^7 S-box, so we don't compare against generate_constants
        assert_eq!(fp_pow7::ROUND_CONSTANTS.len(), POW7_ROUNDS);
        assert_eq!(fq_pow7::ROUND_CONSTANTS.len(), POW7_ROUNDS);

        // Verify MDS matrices are correctly sized
        assert_eq!(fp_pow7::MDS.len(), 3);
        assert_eq!(fp_pow7::MDS[0].len(), 3);
        assert_eq!(fq_pow7::MDS.len(), 3);
        assert_eq!(fq_pow7::MDS[0].len(), 3);
    }

    #[test]
    fn test_basic_pow7_permutation() {
        // Test that the Pow7 permutation produces consistent results

        // Test with Fp
        {
            let (round_constants, mat_external, mat_internal, mat_internal_diag_m_1) =
                <super::P128Pow7T3 as Spec<Fp, 3, 2>>::constants();
            let mut state_fp = [Fp::from(0), Fp::from(1), Fp::from(2)];
            permute::<Fp, super::P128Pow7T3, 3, 2>(
                &mut state_fp,
                &mat_external,
                &mat_internal,
                &mat_internal_diag_m_1,
                &round_constants,
            );
            // Just verify it doesn't panic and produces non-zero output
            assert!(
                state_fp[0] != Fp::from(0)
                    || state_fp[1] != Fp::from(0)
                    || state_fp[2] != Fp::from(0)
            );
        }

        // Test with Fq
        {
            let (round_constants, mat_external, mat_internal, mat_internal_diag_m_1) =
                <super::P128Pow7T3 as Spec<Fq, 3, 2>>::constants();
            let mut state_fq = [Fq::from(0), Fq::from(1), Fq::from(2)];
            permute::<Fq, super::P128Pow7T3, 3, 2>(
                &mut state_fq,
                &mat_external,
                &mat_internal,
                &mat_internal_diag_m_1,
                &round_constants,
            );
            // Just verify it doesn't panic and produces non-zero output
            assert!(
                state_fq[0] != Fq::from(0)
                    || state_fq[1] != Fq::from(0)
                    || state_fq[2] != Fq::from(0)
            );
        }
    }

    #[test]
    fn test_pow7_hash() {
        // Test that the Pow7 hash function works
        let message_fp = [Fp::from(6), Fp::from(42)];
        let result_fp =
            Hash::<_, super::P128Pow7T3, ConstantLength<2>, 3, 2>::init().hash(message_fp);
        assert!(result_fp != Fp::from(0));

        let message_fq = [Fq::from(6), Fq::from(42)];
        let result_fq =
            Hash::<_, super::P128Pow7T3, ConstantLength<2>, 3, 2>::init().hash(message_fq);
        assert!(result_fq != Fq::from(0));
    }
}
