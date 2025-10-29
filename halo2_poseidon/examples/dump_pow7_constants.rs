use core::marker::PhantomData;

use ff::{Field, PrimeField};
use halo2_poseidon::{generate_constants, Mds, Spec};
use pasta_curves::{pallas::Base as Fp, vesta::Base as Fq};

#[derive(Debug)]
struct Pow7Fp<const SECURE_MDS: usize>(PhantomData<Fp>);

impl<const SECURE_MDS: usize> Spec<Fp, 3, 2> for Pow7Fp<SECURE_MDS> {
    fn full_rounds() -> usize {
        8
    }

    fn partial_rounds() -> usize {
        48
    }

    fn sbox(val: Fp) -> Fp {
        val.pow_vartime([7])
    }

    fn secure_mds() -> usize {
        SECURE_MDS
    }

    fn constants() -> (Vec<[Fp; 3]>, Mds<Fp, 3>, Mds<Fp, 3>) {
        generate_constants::<Fp, Self, 3, 2>()
    }
}

#[derive(Debug)]
struct Pow7Fq<const SECURE_MDS: usize>(PhantomData<Fq>);

impl<const SECURE_MDS: usize> Spec<Fq, 3, 2> for Pow7Fq<SECURE_MDS> {
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
        SECURE_MDS
    }

    fn constants() -> (Vec<[Fq; 3]>, Mds<Fq, 3>, Mds<Fq, 3>) {
        generate_constants::<Fq, Self, 3, 2>()
    }
}

fn to_raw_u64s<F: PrimeField>(value: &F) -> [u64; 4] {
    let repr = value.to_repr();
    let bytes = repr.as_ref();
    let mut limbs = [0u64; 4];

    for (idx, limb) in limbs.iter_mut().enumerate() {
        let start = idx * 8;
        let mut chunk = [0u8; 8];
        chunk.copy_from_slice(&bytes[start..start + 8]);
        *limb = u64::from_le_bytes(chunk);
    }

    limbs
}

fn format_hex(limb: u64) -> String {
    let raw = format!("{:016x}", limb);
    let mut formatted = String::from("0x");

    for (idx, chunk) in raw.as_bytes().chunks(4).enumerate() {
        if idx != 0 {
            formatted.push('_');
        }
        formatted.push_str(core::str::from_utf8(chunk).expect("valid hex"));
    }

    formatted
}

fn print_constants_fp() {
    let (round_constants, _mds, _mds_inv) = Pow7Fp::<0>::constants();

    println!(
        "pub(crate) const ROUND_CONSTANTS_POW7: [[pallas::Base; 3]; {}] = [",
        round_constants.len()
    );

    for row in &round_constants {
        println!("    [");
        for element in row {
            let limbs = to_raw_u64s(element);
            println!("        pallas::Base::from_raw([");
            for limb in limbs {
                println!("            {},", format_hex(limb));
            }
            println!("        ]),");
        }
        println!("    ],");
    }
    println!("];\n");
}

fn print_constants_fq() {
    let (round_constants, _mds, _mds_inv) = Pow7Fq::<0>::constants();

    println!(
        "pub(crate) const ROUND_CONSTANTS_POW7: [[vesta::Base; 3]; {}] = [",
        round_constants.len()
    );

    for row in &round_constants {
        println!("    [");
        for element in row {
            let limbs = to_raw_u64s(element);
            println!("        vesta::Base::from_raw([");
            for limb in limbs {
                println!("            {},", format_hex(limb));
            }
            println!("        ]),");
        }
        println!("    ],");
    }
    println!("];");
}

fn main() {
    print_constants_fp();
    print_constants_fq();
}
