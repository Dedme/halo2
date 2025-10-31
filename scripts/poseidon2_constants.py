"""Generate Poseidon2 constants for halo2 from the Sage reference script.

This script drives the upstream `generate_parameters_grain.sage` generator, parses
its output, and updates the corresponding Rust constant tables inside the
`halo2_poseidon2` crate. Usage example:

```bash
python scripts/poseidon2_constants.py \
    --sage-script /path/to/pasta-hadeshash/code/generate_parameters_grain.sage \
    --field vesta
```

Run once per field (``pallas`` or ``vesta``) or omit the ``--field`` flag to
regenerate both.
"""

import argparse
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List


REPO_ROOT = Path(__file__).resolve().parent.parent


@dataclass(frozen=True)
class FieldConfig:
    name: str
    base_type: str
    rust_file: Path
    sage_parameters: List[str]


FIELD_CONFIGS = {
    "vesta": FieldConfig(
        name="vesta",
        base_type="vesta::Base",
        rust_file=REPO_ROOT / "halo2_poseidon2" / "src" / "fq.rs",
        sage_parameters=[
            "1",
            "0",
            "255",
            "3",
            "8",
            "56",
            "0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001",
        ],
    ),
    "pallas": FieldConfig(
        name="pallas",
        base_type="pallas::Base",
        rust_file=REPO_ROOT / "halo2_poseidon2" / "src" / "fp.rs",
        sage_parameters=[
            "1",
            "0",
            "255",
            "3",
            "8",
            "56",
            "0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001",
        ],
    ),
}


def run_sage(sage_script: Path, parameters: List[str]) -> str:
    command = [
        "sage",
        "-c",
        f'load("{sage_script.as_posix()}"); main({json.dumps(parameters)})',
    ]
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    if result.stderr.strip():
        print(result.stderr)
    return result.stdout


def parse_sage_output(output: str) -> tuple[int, List[int], List[int], int]:
    rc_section_match = re.search(
        r"Round constants for GF\(p\):(?P<section>.+?)Prime number:",
        output,
        flags=re.DOTALL,
    )
    if not rc_section_match:
        raise ValueError("Unable to locate round constants in Sage output")
    rc_section = rc_section_match.group("section")
    rc_hex = re.findall(r"0x[0-9a-fA-F]+", rc_section)
    if len(rc_hex) % 3 != 0:
        raise ValueError("Round constant count is not divisible by rate=3")
    rc_ints = [int(value, 16) for value in rc_hex]

    prime_line = output.split("Prime number:", 1)[1]
    prime_hex_matches = re.findall(r"0x[0-9a-fA-F]+", prime_line)
    if not prime_hex_matches:
        raise ValueError("Unable to determine prime from Sage output")
    prime = int(prime_hex_matches[-1], 16)

    mds_match = re.search(
        r"MDS matrix:(?P<section>.+?)Inverse MDS matrix:",
        output,
        flags=re.DOTALL,
    )
    if not mds_match:
        raise ValueError("Unable to locate MDS matrix in Sage output")
    mds_hex = re.findall(r"0x[0-9a-fA-F]+", mds_match.group("section"))
    if len(mds_hex) != 9:
        raise ValueError("Expected 9 MDS entries for width 3")
    mds_ints = [int(value, 16) for value in mds_hex]

    round_constants_count_match = re.search(
        r"Number of round constants:\s*(\d+)", output
    )
    if not round_constants_count_match:
        raise ValueError("Unable to read round constant count header")
    round_constants_count = int(round_constants_count_match.group(1))

    return round_constants_count, rc_ints, mds_ints, prime


def to_limbs(n: int) -> List[int]:
    bytes_le = n.to_bytes(32, "little")
    return [int.from_bytes(bytes_le[i * 8 : (i + 1) * 8], "little") for i in range(4)]


def fmt_limbs(limbs: List[int]) -> List[str]:
    return [format(value, "#018x") for value in limbs]


def format_diagonal_minus_one(base_type: str, mds: List[int], prime: int) -> str:
    diag_indices = [0, 4, 8]
    diag_values = [mds[index] for index in diag_indices]
    diag_minus_one = [((value - 1) % prime) for value in diag_values]

    lines = [
        f"pub(crate) const MAT_DIAG3_M_1_POSEIDON2: [{base_type}; 3] = [",
    ]
    for value in diag_minus_one:
        lines.append(f"    {base_type}::from_raw([")
        for limb in fmt_limbs(to_limbs(value)):
            lines.append(f"        {limb},")
        lines.append("    ]),")
    lines.append("];\n")
    return "\n".join(lines)


def format_round_constants(base_type: str, rc_ints: List[int]) -> str:
    groups = [rc_ints[i : i + 3] for i in range(0, len(rc_ints), 3)]
    lines = [
        f"pub(crate) const ROUND_CONSTANTS: [[{base_type}; 3]; {len(groups)}] = [",
    ]
    for triple in groups:
        lines.append("    [")
        for value in triple:
            lines.append(f"        {base_type}::from_raw([")
            for limb in fmt_limbs(to_limbs(value)):
                lines.append(f"            {limb},")
            lines.append("        ]),")
        lines.append("    ],")
    lines.append("];\n")
    return "\n".join(lines)


def format_mds(base_type: str, mds_ints: List[int]) -> str:
    lines = [
        f"pub(crate) const MDS: [[{base_type}; 3]; 3] = [",
    ]
    for row_start in range(0, len(mds_ints), 3):
        lines.append("    [")
        for value in mds_ints[row_start : row_start + 3]:
            lines.append(f"        {base_type}::from_raw([")
            for limb in fmt_limbs(to_limbs(value)):
                lines.append(f"            {limb},")
            lines.append("        ]),")
        lines.append("    ],")
    lines.append("];\n")
    return "\n".join(lines)


def replace_block(content: str, block_name: str, replacement: str) -> str:
    start_token = f"pub(crate) const {block_name}"
    start = content.find(start_token)
    if start == -1:
        raise ValueError(f"Unable to locate existing block for {block_name}")

    eq_index = content.find("=", start)
    if eq_index == -1:
        raise ValueError(f"Malformed constant definition for {block_name}")

    i = content.find("[", eq_index)
    if i == -1:
        raise ValueError(f"Unable to locate array start for {block_name}")

    depth = 0
    end = None
    while i < len(content):
        char = content[i]
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                end = content.find(";", i)
                if end == -1:
                    raise ValueError(f"Unable to find terminator for {block_name}")
                end += 1
                while end < len(content) and content[end] in "\r\n":
                    end += 1
                break
        i += 1

    if end is None:
        raise ValueError(f"Failed to compute block end for {block_name}")

    replacement_block = replacement.rstrip("\n") + "\n"
    return content[:start] + replacement_block + content[end:]


def update_rust_file(
    rust_file: Path,
    base_type: str,
    rc_ints: List[int],
    mds_ints: List[int],
    prime: int,
) -> None:
    content = rust_file.read_text()

    diag_block = format_diagonal_minus_one(base_type, mds_ints, prime)
    rc_block = format_round_constants(base_type, rc_ints)
    mds_block = format_mds(base_type, mds_ints)

    content = replace_block(content, "MAT_DIAG3_M_1_POSEIDON2", diag_block)
    content = replace_block(content, "ROUND_CONSTANTS", rc_block)
    content = replace_block(content, "MDS", mds_block)

    rust_file.write_text(content)


def generate_for_field(sage_script: Path, config: FieldConfig) -> None:
    print(f"Generating Poseidon2 constants for {config.name}…")
    output = run_sage(sage_script, config.sage_parameters)
    round_count, rc_ints, mds_ints, prime = parse_sage_output(output)

    if round_count != len(rc_ints):
        raise ValueError(
            f"Round constant count mismatch: header={round_count} actual={len(rc_ints)}"
        )

    update_rust_file(config.rust_file, config.base_type, rc_ints, mds_ints, prime)
    print(f"Updated {config.rust_file} with {round_count} round constants.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sage-script",
        type=Path,
        required=True,
        help="Path to generate_parameters_grain.sage",
    )
    parser.add_argument(
        "--field",
        choices=sorted(FIELD_CONFIGS.keys()),
        action="append",
        help="Field(s) to regenerate. Defaults to all",
    )
    args = parser.parse_args()

    sage_script = args.sage_script.resolve()
    if not sage_script.exists():
        raise FileNotFoundError(f"Sage script not found: {sage_script}")

    fields = args.field or sorted(FIELD_CONFIGS.keys())
    for field in fields:
        generate_for_field(sage_script, FIELD_CONFIGS[field])


if __name__ == "__main__":
    main()