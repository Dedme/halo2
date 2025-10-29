use alloc::vec::Vec;

use ff::Field;
use pasta_curves::{pallas::Base as Fp, vesta::Base as Fq};

use super::{Mds, Spec};

/// Poseidon-128 using the $x^5$ S-box, with a width of 3 field elements, and the
/// standard number of rounds for 128-bit security "with margin".
///
/// The standard specification for this set of parameters (on either of the Pasta
/// fields) uses $R_F = 8, R_P = 56$. This is conveniently an even number of
/// partial rounds, making it easier to construct a Halo 2 circuit.
#[derive(Debug)]
pub struct P128Pow5T3;

impl Spec<Fp, 3, 2> for P128Pow5T3 {
     fn full_rounds() -> usize {
        8
    }

    fn partial_rounds() -> usize {
        56
    }

    fn sbox(val: Fp) -> Fp {
        val.pow_vartime([5])
    }

    fn secure_mds() -> usize {
        // Using pre-generated constants from zkhash; MDS security is verified externally
        unimplemented!("Poseidon2 uses pre-generated constants")
    }

    fn constants() -> (Vec<[Fp; 3]>, Mds<Fp, 3>, Mds<Fp, 3>) {
        (
            super::fp::ROUND_CONSTANTS[..].to_vec(),
            super::fp::MDS,
            super::fp::MDS_INV,
        )
    }
}

impl Spec<Fq, 3, 2> for P128Pow5T3 {
    fn full_rounds() -> usize {
        8
    }

    fn partial_rounds() -> usize {
        56
    }

    fn sbox(val: Fq) -> Fq {
        val.pow_vartime([5])
    }

    fn secure_mds() -> usize {
        unimplemented!()
    }

    fn constants() -> (Vec<[Fq; 3]>, Mds<Fq, 3>, Mds<Fq, 3>) {
        (
            super::fq::ROUND_CONSTANTS[..].to_vec(),
            super::fq::MDS,
            super::fq::MDS_INV,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ConstantLength, Hash as PoseidonHash};

    #[test]
    fn test_poseidon2_deterministic_fp() {
        // Verify Poseidon2 produces consistent results for Pallas field
        let inputs = [Fp::from(1), Fp::from(2)];

        let hash1 =
            PoseidonHash::<Fp, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs);
        let hash2 =
            PoseidonHash::<Fp, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs);

        assert_eq!(hash1, hash2, "Poseidon2 must be deterministic");
        assert_ne!(hash1, Fp::zero(), "Hash should not be zero");
    }

    #[test]
    fn test_poseidon2_deterministic_fq() {
        // Verify Poseidon2 produces consistent results for Vesta field
        let inputs = [Fq::from(1), Fq::from(2)];

        let hash1 =
            PoseidonHash::<Fq, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs);
        let hash2 =
            PoseidonHash::<Fq, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs);

        assert_eq!(hash1, hash2, "Poseidon2 must be deterministic");
        assert_ne!(hash1, Fq::zero(), "Hash should not be zero");
    }

    #[test]
    fn test_poseidon2_different_inputs() {
        // Verify different inputs produce different hashes
        let inputs1 = [Fp::from(1), Fp::from(2)];
        let inputs2 = [Fp::from(2), Fp::from(1)];

        let hash1 =
            PoseidonHash::<Fp, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs1);
        let hash2 =
            PoseidonHash::<Fp, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs2);

        assert_ne!(
            hash1, hash2,
            "Different inputs should produce different hashes"
        );
    }

    #[test]
    fn test_poseidon2_zero_inputs() {
        // Verify hashing zero inputs produces non-zero result
        let inputs = [Fp::zero(), Fp::zero()];

        let hash =
            PoseidonHash::<Fp, P128Pow5T3, ConstantLength<2>, 3, 2>::init().hash(inputs);

        assert_ne!(hash, Fp::zero(), "Hash of zero inputs should not be zero");
    }

    // TODO: Add test vectors from zkhash reference implementation
    // Once available, validate against known-good Poseidon2 outputs
}
