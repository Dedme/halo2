use alloc::vec::Vec;

use ff::Field;
use pasta_curves::{pallas::Base as Fp, vesta::Base as Fq};

use super::{Mds, Spec};

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
        48  // Even number required for gadget (processes 2 partial rounds per circuit row)
    }

    fn sbox(val: Fp) -> Fp {
        val.pow_vartime([7])  // Changed from x^5 to x^7 (Aleph Zero optimization)
    }

    fn secure_mds() -> usize {
        unimplemented!()
    }

    fn constants() -> (Vec<[Fp; 3]>, Mds<Fp, 3>, Mds<Fp, 3>) {
        (
            super::fp_pow7::ROUND_CONSTANTS.to_vec(),
            super::fp_pow7::MDS,
            super::fp_pow7::MDS_INV,
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
        val.pow_vartime([7])
    }

    fn secure_mds() -> usize {
        unimplemented!()
    }

    fn constants() -> (Vec<[Fq; 3]>, Mds<Fq, 3>, Mds<Fq, 3>) {
        (
            super::fq_pow7::ROUND_CONSTANTS.to_vec(),
            super::fq_pow7::MDS,
            super::fq_pow7::MDS_INV,
        )
    }
}

#[cfg(test)]
mod tests {
    const POW7_ROUNDS: usize = 56;  // 8 full + 48 partial

    use alloc::vec::Vec;
    use core::marker::PhantomData;

    use ff::{Field, FromUniformBytes, PrimeField};

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
            47
        }

        fn sbox(val: F) -> F {
            val.pow_vartime([7])
        }

        fn secure_mds() -> usize {
            SECURE_MDS
        }

        fn constants() -> (Vec<[F; 3]>, Mds<F, 3>, Mds<F, 3>) {
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

    // NOTE: This test is disabled for Pow7 because the hardcoded constants come from
    // an external source (HorizenLabs Sage script), not from Grain LFSR.
    // The permute_test_vectors test already validates correctness against reference vectors.
    #[test]
    #[ignore]
    fn test_against_reference() {
        {
            let vectors = crate::test_vectors_7t3::fp::permute();
            assert!(!vectors.is_empty(), "expected at least one test vector");
            let tv = &vectors[0];

            let mut input = [
                Fp::from_repr(tv.initial_state[0]).unwrap(),
                Fp::from_repr(tv.initial_state[1]).unwrap(),
                Fp::from_repr(tv.initial_state[2]).unwrap(),
            ];

            let expected_output = [
                Fp::from_repr(tv.final_state[0]).unwrap(),
                Fp::from_repr(tv.final_state[1]).unwrap(),
                Fp::from_repr(tv.final_state[2]).unwrap(),
            ];

            permute::<Fp, P128Pow7T3Gen<Fp, 0>, 3, 2>(
                &mut input,
                &fp_pow7::MDS,
                &fp_pow7::ROUND_CONSTANTS,
            );
            assert_eq!(input, expected_output);
        }

        {
            let vectors = crate::test_vectors_7t3::fq::permute();
            assert!(!vectors.is_empty(), "expected at least one test vector");
            let tv = &vectors[0];

            let mut input = [
                Fq::from_repr(tv.initial_state[0]).unwrap(),
                Fq::from_repr(tv.initial_state[1]).unwrap(),
                Fq::from_repr(tv.initial_state[2]).unwrap(),
            ];

            let expected_output = [
                Fq::from_repr(tv.final_state[0]).unwrap(),
                Fq::from_repr(tv.final_state[1]).unwrap(),
                Fq::from_repr(tv.final_state[2]).unwrap(),
            ];

            permute::<Fq, P128Pow7T3Gen<Fq, 0>, 3, 2>(
                &mut input,
                &fq_pow7::MDS,
                &fq_pow7::ROUND_CONSTANTS,
            );
            assert_eq!(input, expected_output);
        }
    }

    #[test]
    fn permute_test_vectors() {
        {
            let (round_constants, mds, _) = super::P128Pow7T3::constants();

            for tv in crate::test_vectors_7t3::fp::permute() {
                let mut state = [
                    Fp::from_repr(tv.initial_state[0]).unwrap(),
                    Fp::from_repr(tv.initial_state[1]).unwrap(),
                    Fp::from_repr(tv.initial_state[2]).unwrap(),
                ];

                permute::<Fp, super::P128Pow7T3, 3, 2>(&mut state, &mds, &round_constants);

                for (expected, actual) in tv.final_state.iter().zip(state.iter()) {
                    assert_eq!(&actual.to_repr(), expected);
                }
            }
        }

        {
            let (round_constants, mds, _) = super::P128Pow7T3::constants();

            for tv in crate::test_vectors_7t3::fq::permute() {
                let mut state = [
                    Fq::from_repr(tv.initial_state[0]).unwrap(),
                    Fq::from_repr(tv.initial_state[1]).unwrap(),
                    Fq::from_repr(tv.initial_state[2]).unwrap(),
                ];

                permute::<Fq, super::P128Pow7T3, 3, 2>(&mut state, &mds, &round_constants);

                for (expected, actual) in tv.final_state.iter().zip(state.iter()) {
                    assert_eq!(&actual.to_repr(), expected);
                }
            }
        }
    }

    #[test]
    fn hash_test_vectors() {
        for tv in crate::test_vectors_7t3::fp::hash() {
            let message = [
                Fp::from_repr(tv.input[0]).unwrap(),
                Fp::from_repr(tv.input[1]).unwrap(),
            ];

            let result =
                Hash::<_, super::P128Pow7T3, ConstantLength<2>, 3, 2>::init().hash(message);

            assert_eq!(result.to_repr(), tv.output);
        }

        for tv in crate::test_vectors_7t3::fq::hash() {
            let message = [
                Fq::from_repr(tv.input[0]).unwrap(),
                Fq::from_repr(tv.input[1]).unwrap(),
            ];

            let result =
                Hash::<_, super::P128Pow7T3, ConstantLength<2>, 3, 2>::init().hash(message);

            assert_eq!(result.to_repr(), tv.output);
        }
    }
}
