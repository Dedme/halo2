use super::circuit::{self, Any, Column};
use crate::{
    arithmetic::CurveAffine,
    helpers::{CurveRead, CurveWrite},
    io_utils,
    poly::{Coeff, ExtendedLagrangeCoeff, LagrangeCoeff, Polynomial},
};
use group::ff::PrimeField;
use std::io::{self, Read, Write};

pub(crate) mod keygen;
pub(crate) mod prover;
pub(crate) mod verifier;

/// A permutation argument.
#[derive(Debug, Clone)]
pub(crate) struct Argument {
    /// A sequence of columns involved in the argument.
    columns: Vec<Column<Any>>,
}

impl Argument {
    pub(crate) fn new() -> Self {
        Argument { columns: vec![] }
    }

    /// Returns the minimum circuit degree required by the permutation argument.
    /// The argument may use larger degree gates depending on the actual
    /// circuit's degree and how many columns are involved in the permutation.
    pub(crate) fn required_degree(&self) -> usize {
        // degree 2:
        // l_0(X) * (1 - z(X)) = 0
        //
        // We will fit as many polynomials p_i(X) as possible
        // into the required degree of the circuit, so the
        // following will not affect the required degree of
        // this middleware.
        //
        // (1 - (l_last(X) + l_blind(X))) * (
        //   z(\omega X) \prod (p(X) + \beta s_i(X) + \gamma)
        // - z(X) \prod (p(X) + \delta^i \beta X + \gamma)
        // )
        //
        // On the first sets of columns, except the first
        // set, we will do
        //
        // l_0(X) * (z(X) - z'(\omega^(last) X)) = 0
        //
        // where z'(X) is the permutation for the previous set
        // of columns.
        //
        // On the final set of columns, we will do
        //
        // degree 3:
        // l_last(X) * (z'(X)^2 - z'(X)) = 0
        //
        // which will allow the last value to be zero to
        // ensure the argument is perfectly complete.

        // There are constraints of degree 3 regardless of the
        // number of columns involved.
        3
    }

    pub(crate) fn add_column(&mut self, column: Column<Any>) {
        if !self.columns.contains(&column) {
            self.columns.push(column);
        }
    }

    pub(crate) fn get_columns(&self) -> Vec<Column<Any>> {
        self.columns.clone()
    }

    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        let columns = self.get_columns();
        io_utils::write_vec(writer, &columns, |w, column| {
            io_utils::write_usize(w, column.index())?;
            let column_type = match column.column_type() {
                Any::Advice => 0u8,
                Any::Fixed => 1u8,
                Any::Instance => 2u8,
            };
            io_utils::write_u8(w, column_type)
        })
    }

    pub(crate) fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let columns = io_utils::read_vec(reader, |r| {
            let index = io_utils::read_usize(r)?;
            let column_type = match io_utils::read_u8(r)? {
                0 => Any::Advice,
                1 => Any::Fixed,
                2 => Any::Instance,
                other => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid column type tag {other}"),
                    ))
                }
            };
            Ok(circuit::make_any_column(index, column_type))
        })?;

        let mut argument = Argument::new();
        for column in columns {
            argument.add_column(column);
        }

        Ok(argument)
    }
}

/// The verifying key for a single permutation argument.
#[derive(Clone, Debug)]
pub(crate) struct VerifyingKey<C: CurveAffine> {
    commitments: Vec<C>,
}

/// The proving key for a single permutation argument.
#[derive(Clone, Debug)]
pub(crate) struct ProvingKey<C: CurveAffine> {
    permutations: Vec<Polynomial<C::Scalar, LagrangeCoeff>>,
    polys: Vec<Polynomial<C::Scalar, Coeff>>,
    pub(super) cosets: Vec<Polynomial<C::Scalar, ExtendedLagrangeCoeff>>,
}

impl<C: CurveAffine> VerifyingKey<C> {
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        io_utils::write_vec(writer, &self.commitments, |w, commitment| {
            commitment.write(w)
        })
    }

    pub(crate) fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let commitments = io_utils::read_vec(reader, |r| C::read(r))?;
        Ok(VerifyingKey { commitments })
    }
}

impl<C: CurveAffine> ProvingKey<C>
where
    C::Scalar: PrimeField,
{
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        io_utils::write_vec(writer, &self.permutations, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })?;
        io_utils::write_vec(writer, &self.polys, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })?;
        io_utils::write_vec(writer, &self.cosets, |w, poly| {
            io_utils::write_polynomial(w, poly)
        })
    }

    pub(crate) fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let permutations = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;
        let polys = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;
        let cosets = io_utils::read_vec(reader, |r| io_utils::read_polynomial(r))?;

        Ok(ProvingKey {
            permutations,
            polys,
            cosets,
        })
    }
}
