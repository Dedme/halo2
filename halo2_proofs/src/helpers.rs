use std::io;

use pasta_curves::arithmetic::CurveAffine;

pub(crate) trait CurveRead: CurveAffine {
    /// Reads a compressed element from the buffer and attempts to parse it
    /// using `from_bytes`.
    fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        let mut compressed = Self::Repr::default();
        reader.read_exact(compressed.as_mut())?;
        Option::from(Self::from_bytes(&compressed))
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "invalid point encoding in proof"))
    }
}

impl<C: CurveAffine> CurveRead for C {}

pub(crate) trait CurveWrite: CurveAffine {
    /// Writes the compressed element to the provided buffer using `to_bytes`.
    fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(self.to_bytes().as_ref())
    }
}

impl<C: CurveAffine> CurveWrite for C {}
