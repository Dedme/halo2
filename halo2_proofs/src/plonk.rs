//! This module provides an implementation of a variant of (Turbo)[PLONK][plonk]
//! that is designed specifically for the polynomial commitment scheme described
//! in the [Halo][halo] paper.
//!
//! [halo]: https://eprint.iacr.org/2019/1021
//! [plonk]: https://eprint.iacr.org/2019/953

use blake2b_simd::Params as Blake2bParams;
use group::ff::{Field, FromUniformBytes, PrimeField};

use crate::arithmetic::CurveAffine;
use crate::helpers::{CurveRead, CurveWrite};
use crate::io_utils;
use crate::poly::{
    Coeff, EvaluationDomain, ExtendedLagrangeCoeff, LagrangeCoeff, PinnedEvaluationDomain,
    Polynomial,
};
use crate::transcript::{ChallengeScalar, EncodedChallenge, Transcript};

mod assigned;
mod circuit;
mod error;
mod keygen;
mod lookup;
pub(crate) mod permutation;
mod vanishing;

mod prover;
mod verifier;

pub use assigned::*;
pub use circuit::*;
pub use error::*;
pub use keygen::*;
pub use prover::*;
pub use verifier::*;

use std::io::{self, Read, Write};

/// This is a verifying key which allows for the verification of proofs for a
/// particular circuit.
#[derive(Clone, Debug)]
pub struct VerifyingKey<C: CurveAffine> {
    domain: EvaluationDomain<C::Scalar>,
    fixed_commitments: Vec<C>,
    permutation: permutation::VerifyingKey<C>,
    cs: ConstraintSystem<C::Scalar>,
    /// Cached maximum degree of `cs` (which doesn't change after construction).
    cs_degree: usize,
    /// The representative of this `VerifyingKey` in transcripts.
    transcript_repr: C::Scalar,
}

const VERIFICATION_KEY_SERIALIZATION_VERSION: u32 = 1;
const PROVING_KEY_SERIALIZATION_VERSION: u32 = 1;

impl<C: CurveAffine> VerifyingKey<C>
where
    C::Scalar: PrimeField + FromUniformBytes<64>,
{
    fn from_parts(
        domain: EvaluationDomain<C::Scalar>,
        fixed_commitments: Vec<C>,
        permutation: permutation::VerifyingKey<C>,
        cs: ConstraintSystem<C::Scalar>,
    ) -> Self {
        // Compute cached values.
        let cs_degree = cs.degree();

        let mut vk = Self {
            domain,
            fixed_commitments,
            permutation,
            cs,
            cs_degree,
            // Temporary, this is not pinned.
            transcript_repr: C::Scalar::ZERO,
        };

        let mut hasher = Blake2bParams::new()
            .hash_length(64)
            .personal(b"Halo2-Verify-Key")
            .to_state();

        let s = format!("{:?}", vk.pinned());

        hasher.update(&(s.len() as u64).to_le_bytes());
        hasher.update(s.as_bytes());

        // Hash in final Blake2bState
        vk.transcript_repr = C::Scalar::from_uniform_bytes(hasher.finalize().as_array());

        vk
    }

    /// Serializes the verifying key to the provided writer.
    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        io_utils::write_u32(writer, VERIFICATION_KEY_SERIALIZATION_VERSION)?;
        io_utils::write_u32(writer, self.domain.get_k())?;
        let j = (self.domain.get_quotient_poly_degree() as u32) + 1;
        io_utils::write_u32(writer, j)?;
        io_utils::write_vec(writer, &self.fixed_commitments, |w, commitment| {
            commitment.write(w)
        })?;
        self.permutation.write(writer)?;
        self.cs.write(writer)
    }

    /// Deserializes a verifying key from the provided reader.
    pub fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let version = io_utils::read_u32(reader)?;
        if version != VERIFICATION_KEY_SERIALIZATION_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported verifying key version {version}"),
            ));
        }

        let k = io_utils::read_u32(reader)?;
        let j = io_utils::read_u32(reader)?;
        let domain = EvaluationDomain::new(j, k);
        let fixed_commitments = io_utils::read_vec(reader, |r| C::read(r))?;
        let permutation = permutation::VerifyingKey::read(reader)?;
        let cs = ConstraintSystem::read(reader)?;

        Ok(VerifyingKey::from_parts(
            domain,
            fixed_commitments,
            permutation,
            cs,
        ))
    }
}

impl<C: CurveAffine> VerifyingKey<C> {
    /// Hashes a verification key into a transcript.
    pub fn hash_into<E: EncodedChallenge<C>, T: Transcript<C, E>>(
        &self,
        transcript: &mut T,
    ) -> io::Result<()> {
        transcript.common_scalar(self.transcript_repr)?;

        Ok(())
    }

    /// Obtains a pinned representation of this verification key that contains
    /// the minimal information necessary to reconstruct the verification key.
    pub fn pinned(&self) -> PinnedVerificationKey<'_, C> {
        PinnedVerificationKey {
            base_modulus: C::Base::MODULUS,
            scalar_modulus: C::Scalar::MODULUS,
            domain: self.domain.pinned(),
            fixed_commitments: &self.fixed_commitments,
            permutation: &self.permutation,
            cs: self.cs.pinned(),
        }
    }
}

/// Minimal representation of a verification key that can be used to identify
/// its active contents.
#[allow(dead_code)]
#[derive(Debug)]
pub struct PinnedVerificationKey<'a, C: CurveAffine> {
    base_modulus: &'static str,
    scalar_modulus: &'static str,
    domain: PinnedEvaluationDomain<'a, C::Scalar>,
    cs: PinnedConstraintSystem<'a, C::Scalar>,
    fixed_commitments: &'a Vec<C>,
    permutation: &'a permutation::VerifyingKey<C>,
}
/// This is a proving key which allows for the creation of proofs for a
/// particular circuit.
#[derive(Clone, Debug)]
pub struct ProvingKey<C: CurveAffine> {
    vk: VerifyingKey<C>,
    l0: Polynomial<C::Scalar, ExtendedLagrangeCoeff>,
    l_blind: Polynomial<C::Scalar, ExtendedLagrangeCoeff>,
    l_last: Polynomial<C::Scalar, ExtendedLagrangeCoeff>,
    fixed_values: Vec<Polynomial<C::Scalar, LagrangeCoeff>>,
    fixed_polys: Vec<Polynomial<C::Scalar, Coeff>>,
    fixed_cosets: Vec<Polynomial<C::Scalar, ExtendedLagrangeCoeff>>,
    permutation: permutation::ProvingKey<C>,
}

impl<C: CurveAffine> ProvingKey<C> {
    /// Get the underlying [`VerifyingKey`].
    pub fn get_vk(&self) -> &VerifyingKey<C> {
        &self.vk
    }
}

impl<C: CurveAffine> ProvingKey<C>
where
    C::Scalar: PrimeField + FromUniformBytes<64>,
{
    /// Serializes the proving key, including its verifying key, to the writer.
    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        io_utils::write_u32(writer, PROVING_KEY_SERIALIZATION_VERSION)?;
        self.vk.write(writer)?;
        io_utils::write_polynomial(writer, &self.l0)?;
        io_utils::write_polynomial(writer, &self.l_blind)?;
        io_utils::write_polynomial(writer, &self.l_last)?;
        io_utils::write_vec(writer, &self.fixed_values, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })?;
        io_utils::write_vec(writer, &self.fixed_polys, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })?;
        io_utils::write_vec(writer, &self.fixed_cosets, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })?;
        self.permutation.write(writer)
    }

    /// Deserializes a proving key, reconstructing the embedded verifying key as well.
    pub fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let version = io_utils::read_u32(reader)?;
        if version != PROVING_KEY_SERIALIZATION_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported proving key version {version}"),
            ));
        }

        let vk = VerifyingKey::read(reader)?;
        let l0 = io_utils::read_polynomial(reader)?;
        let l_blind = io_utils::read_polynomial(reader)?;
        let l_last = io_utils::read_polynomial(reader)?;
        let fixed_values = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;
        let fixed_polys = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;
        let fixed_cosets = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;
        let permutation = permutation::ProvingKey::read(reader)?;

        Ok(ProvingKey {
            vk,
            l0,
            l_blind,
            l_last,
            fixed_values,
            fixed_polys,
            fixed_cosets,
            permutation,
        })
    }
}

impl<C: CurveAffine> VerifyingKey<C> {
    /// Get the underlying [`EvaluationDomain`].
    pub fn get_domain(&self) -> &EvaluationDomain<C::Scalar> {
        &self.domain
    }
}

#[derive(Clone, Copy, Debug)]
struct Theta;
type ChallengeTheta<F> = ChallengeScalar<F, Theta>;

#[derive(Clone, Copy, Debug)]
struct Beta;
type ChallengeBeta<F> = ChallengeScalar<F, Beta>;

#[derive(Clone, Copy, Debug)]
struct Gamma;
type ChallengeGamma<F> = ChallengeScalar<F, Gamma>;

#[derive(Clone, Copy, Debug)]
struct Y;
type ChallengeY<F> = ChallengeScalar<F, Y>;

#[derive(Clone, Copy, Debug)]
struct X;
type ChallengeX<F> = ChallengeScalar<F, X>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::SimpleFloorPlanner;
    use crate::plonk::{keygen_pk, keygen_vk};
    use crate::poly::{commitment::Params, Basis};
    use group::ff::PrimeField;
    use pasta_curves::EqAffine;

    #[derive(Clone, Copy)]
    struct MyCircuit;

    impl<F: Field> Circuit<F> for MyCircuit {
        type Config = ();
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            *self
        }

        fn configure(_meta: &mut ConstraintSystem<F>) -> Self::Config {}

        fn synthesize(
            &self,
            _config: Self::Config,
            _layouter: impl crate::circuit::Layouter<F>,
        ) -> Result<(), Error> {
            Ok(())
        }
    }

    fn assert_polynomials_equal<F, B>(lhs: &Polynomial<F, B>, rhs: &Polynomial<F, B>)
    where
        F: PrimeField,
        B: Basis,
    {
        let lhs_values: Vec<F> = lhs.iter().copied().collect();
        let rhs_values: Vec<F> = rhs.iter().copied().collect();
        assert_eq!(lhs_values, rhs_values);
    }

    #[test]
    fn verifying_and_proving_key_roundtrip() {
        let params: Params<EqAffine> = Params::new(3);
        let vk = keygen_vk(&params, &MyCircuit).expect("keygen_vk should not fail");
        let pk = keygen_pk(&params, vk.clone(), &MyCircuit).expect("keygen_pk should not fail");

        let mut vk_bytes = Vec::new();
        vk.write(&mut vk_bytes)
            .expect("verifying key serialization");
        let mut vk_reader = &vk_bytes[..];
        let vk_roundtrip =
            VerifyingKey::read(&mut vk_reader).expect("verifying key deserialization");

        assert_eq!(
            format!("{:?}", vk.pinned()),
            format!("{:?}", vk_roundtrip.pinned())
        );
        assert_eq!(vk.fixed_commitments, vk_roundtrip.fixed_commitments);
        assert_eq!(vk.cs_degree, vk_roundtrip.cs_degree);
        assert_eq!(vk.transcript_repr, vk_roundtrip.transcript_repr);
        assert_eq!(
            format!("{:?}", vk.cs.pinned()),
            format!("{:?}", vk_roundtrip.cs.pinned())
        );

        let mut pk_bytes = Vec::new();
        pk.write(&mut pk_bytes).expect("proving key serialization");
        let mut pk_reader = &pk_bytes[..];
        let pk_roundtrip: ProvingKey<EqAffine> =
            ProvingKey::read(&mut pk_reader).expect("proving key deserialization");

        assert_eq!(
            format!("{:?}", pk.vk.pinned()),
            format!("{:?}", pk_roundtrip.vk.pinned())
        );
        assert_polynomials_equal(&pk.l0, &pk_roundtrip.l0);
        assert_polynomials_equal(&pk.l_blind, &pk_roundtrip.l_blind);
        assert_polynomials_equal(&pk.l_last, &pk_roundtrip.l_last);

        assert_eq!(pk.fixed_values.len(), pk_roundtrip.fixed_values.len());
        for (expected, actual) in pk.fixed_values.iter().zip(pk_roundtrip.fixed_values.iter()) {
            assert_polynomials_equal(expected, actual);
        }

        assert_eq!(pk.fixed_polys.len(), pk_roundtrip.fixed_polys.len());
        for (expected, actual) in pk.fixed_polys.iter().zip(pk_roundtrip.fixed_polys.iter()) {
            assert_polynomials_equal(expected, actual);
        }

        assert_eq!(pk.fixed_cosets.len(), pk_roundtrip.fixed_cosets.len());
        for (expected, actual) in pk.fixed_cosets.iter().zip(pk_roundtrip.fixed_cosets.iter()) {
            assert_polynomials_equal(expected, actual);
        }

        assert_eq!(
            format!("{:?}", pk.permutation),
            format!("{:?}", pk_roundtrip.permutation)
        );
    }
}
