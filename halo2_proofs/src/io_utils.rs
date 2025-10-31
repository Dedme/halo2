use std::convert::TryFrom;
use std::io::{self, Read, Write};

use group::ff::PrimeField;

use crate::poly::{Basis, Polynomial};

pub(crate) fn write_u32<W: Write>(writer: &mut W, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub(crate) fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn write_u64<W: Write>(writer: &mut W, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub(crate) fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

pub(crate) fn write_usize<W: Write>(writer: &mut W, value: usize) -> io::Result<()> {
    let value = u64::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "value does not fit into u64 during serialization",
        )
    })?;
    write_u64(writer, value)
}

pub(crate) fn read_usize<R: Read>(reader: &mut R) -> io::Result<usize> {
    let value = read_u64(reader)?;
    usize::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "value does not fit into usize during deserialization",
        )
    })
}

pub(crate) fn write_i32<W: Write>(writer: &mut W, value: i32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

pub(crate) fn read_i32<R: Read>(reader: &mut R) -> io::Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

pub(crate) fn write_bool<W: Write>(writer: &mut W, value: bool) -> io::Result<()> {
    writer.write_all(&[value as u8])
}

pub(crate) fn read_bool<R: Read>(reader: &mut R) -> io::Result<bool> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    match buf[0] {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "boolean value must be 0 or 1",
        )),
    }
}

pub(crate) fn write_u8<W: Write>(writer: &mut W, value: u8) -> io::Result<()> {
    writer.write_all(&[value])
}

pub(crate) fn read_u8<R: Read>(reader: &mut R) -> io::Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

pub(crate) fn write_len<W: Write>(writer: &mut W, len: usize) -> io::Result<()> {
    let len = u32::try_from(len).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "length exceeds u32::MAX during serialization",
        )
    })?;
    write_u32(writer, len)
}

pub(crate) fn read_len<R: Read>(reader: &mut R) -> io::Result<usize> {
    let len = read_u32(reader)?;
    usize::try_from(len).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "length does not fit into usize during deserialization",
        )
    })
}

pub(crate) fn write_field<F: PrimeField, W: Write>(writer: &mut W, value: &F) -> io::Result<()> {
    writer.write_all(value.to_repr().as_ref())
}

pub(crate) fn read_field<F: PrimeField, R: Read>(reader: &mut R) -> io::Result<F> {
    let mut repr = F::Repr::default();
    reader.read_exact(repr.as_mut())?;
    Option::from(F::from_repr(repr)).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid field element encoding during deserialization",
        )
    })
}

pub(crate) fn write_vec<T, W>(
    writer: &mut W,
    values: &[T],
    mut write_item: impl FnMut(&mut W, &T) -> io::Result<()>,
) -> io::Result<()>
where
    W: Write,
{
    write_len(writer, values.len())?;
    for value in values {
        write_item(writer, value)?;
    }
    Ok(())
}

pub(crate) fn read_vec<T, R>(
    reader: &mut R,
    mut read_item: impl FnMut(&mut R) -> io::Result<T>,
) -> io::Result<Vec<T>>
where
    R: Read,
{
    let len = read_len(reader)?;
    let mut values = Vec::with_capacity(len);
    for _ in 0..len {
        values.push(read_item(reader)?);
    }
    Ok(values)
}

pub(crate) fn write_str<W: Write>(writer: &mut W, value: &str) -> io::Result<()> {
    write_len(writer, value.len())?;
    writer.write_all(value.as_bytes())
}

pub(crate) fn read_string<R: Read>(reader: &mut R) -> io::Result<String> {
    let len = read_len(reader)?;
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    String::from_utf8(buf).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid UTF-8 string during deserialization: {err}"),
        )
    })
}

pub(crate) fn write_polynomial<F, B, W>(writer: &mut W, poly: &Polynomial<F, B>) -> io::Result<()>
where
    F: PrimeField,
    B: Basis,
    W: Write,
{
    write_len(writer, poly.num_coeffs())?;
    for value in poly.iter() {
        write_field(writer, value)?;
    }
    Ok(())
}

pub(crate) fn read_polynomial<F, B, R>(reader: &mut R) -> io::Result<Polynomial<F, B>>
where
    F: PrimeField,
    B: Basis,
    R: Read,
{
    let len = read_len(reader)?;
    let mut values = Vec::with_capacity(len);
    for _ in 0..len {
        values.push(read_field(reader)?);
    }
    Ok(Polynomial::from_raw(values))
}
