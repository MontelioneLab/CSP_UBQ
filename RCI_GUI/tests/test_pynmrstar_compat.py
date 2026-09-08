import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPTS_DIR))

from format_conversion import convert_file, detect_supported_format, list_input_entities, parse_supported_input


EXAMPLE_DIR = ROOT / "RCI-Wishart-calculation" / "BMRB-to-RCI-examples" / "2JVD_bmr15476-nmrstar2-example"
OUTPUT_SUFFIXES = ["RCI.txt", "S2.txt", "MD_RMSD.txt", "NMR_RMSD.txt"]
SHIFTY_CONVERTER = ROOT / "RCI-Wishart-calculation" / "Janet-conversion-scripts" / "nmrstar3toSHIFTY.py"
SCRIPT_PYTHON = ROOT / ".venv" / "bin" / "python"
ONE_TO_THREE = {
	"A": "ALA",
	"C": "CYS",
	"D": "ASP",
	"E": "GLU",
	"F": "PHE",
	"G": "GLY",
	"H": "HIS",
	"I": "ILE",
	"K": "LYS",
	"L": "LEU",
	"M": "MET",
	"N": "ASN",
	"P": "PRO",
	"Q": "GLN",
	"R": "ARG",
	"S": "SER",
	"T": "THR",
	"V": "VAL",
	"W": "TRP",
	"Y": "TYR",
}


def run_script(script_name, sample_name, sample_path=None, extra_args=None):
	script_path = ROOT / script_name
	if sample_path is None:
		sample_path = EXAMPLE_DIR / sample_name
	if extra_args is None:
		extra_args = []
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		shutil.copy2(sample_path, temp_path / sample_name)
		process = subprocess.run(
			[str(SCRIPT_PYTHON), str(script_path), "-b", sample_name, *extra_args],
			cwd=temp_path,
			text=True,
			capture_output=True,
		)
		outputs = {}
		for suffix in OUTPUT_SUFFIXES:
			output_path = temp_path / f"{sample_name}.{suffix}"
			if output_path.exists():
				outputs[suffix] = output_path.read_text()
		return process, outputs


def _assert_output_text_close(left_text: str, right_text: str, *, abs_tol: float = 1e-12) -> None:
	left_lines = [line.split() for line in left_text.splitlines() if line.strip()]
	right_lines = [line.split() for line in right_text.splitlines() if line.strip()]
	assert len(left_lines) == len(right_lines)
	for left_parts, right_parts in zip(left_lines, right_lines):
		assert len(left_parts) == len(right_parts)
		for left_token, right_token in zip(left_parts, right_parts):
			try:
				left_value = float(left_token)
				right_value = float(right_token)
			except ValueError:
				assert left_token == right_token
				continue
			assert math.isclose(left_value, right_value, rel_tol=0.0, abs_tol=abs_tol)


def _assert_output_sets_close(left_outputs: dict[str, str], right_outputs: dict[str, str]) -> None:
	assert set(left_outputs) == set(OUTPUT_SUFFIXES)
	assert set(right_outputs) == set(OUTPUT_SUFFIXES)
	for suffix in OUTPUT_SUFFIXES:
		_assert_output_text_close(left_outputs[suffix], right_outputs[suffix])


def _make_headerless_three_letter_shifty(source_path: Path, output_path: Path) -> None:
	lines = []
	for raw_line in source_path.read_text().splitlines():
		stripped = raw_line.strip()
		if (not stripped) or stripped.startswith("#NUM") or stripped.startswith("#SEQ"):
			continue
		parts = stripped.split()
		if parts[0].startswith("#"):
			continue
		parts[1] = ONE_TO_THREE[parts[1]]
		lines.append("\t".join(parts))
	output_path.write_text("\n".join(lines) + "\n")


def test_active_script_matches_between_original_nmrstar_versions():
	v21_process, v21_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", "bmr15476_21.str.txt")
	assert v21_process.returncode == 0, v21_process.stderr or v21_process.stdout

	v3_process, v3_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", "bmr15476_3.str.txt")
	assert v3_process.returncode == 0, v3_process.stderr or v3_process.stdout

	_assert_output_sets_close(v21_outputs, v3_outputs)


def test_shifty_matches_converter_output():
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		generated_shifty = temp_path / "generated-shifty.txt"
		subprocess.run(
			[
				str(SCRIPT_PYTHON),
				str(SHIFTY_CONVERTER),
				str(EXAMPLE_DIR / "bmr15476_3.str.txt"),
			],
			text=True,
			check=True,
			stdout=generated_shifty.open("w"),
		)

		generated_process, generated_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", generated_shifty.name, generated_shifty)
		bundled_process, bundled_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", "15476-shifty.txt")

	assert generated_process.returncode == 0, generated_process.stderr or generated_process.stdout
	assert bundled_process.returncode == 0, bundled_process.stderr or bundled_process.stdout
	_assert_output_sets_close(generated_outputs, bundled_outputs)


def test_script_outputs_match_across_converted_supported_formats():
	source_v3 = ROOT / "inputs" / "bmr15476_3.str.txt"
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		generated_v21 = temp_path / "generated_v21.str.txt"
		generated_shifty = temp_path / "generated_shifty.txt"

		convert_file(source_v3, "NMRSTAR21", generated_v21)
		convert_file(source_v3, "SHIFTY", generated_shifty)

		assert detect_supported_format(source_v3) == "NMRSTAR3"
		assert detect_supported_format(generated_v21) == "NMRSTAR21"
		assert detect_supported_format(generated_shifty) == "SHIFTY"
		first_line = generated_shifty.read_text().splitlines()[0]
		assert first_line == "#NUM\tAA\tHA\tCA\tCB\tCO\tN\tHN"

		v3_process, v3_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", source_v3.name, source_v3)
		v21_process, v21_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", generated_v21.name, generated_v21)
		shifty_process, shifty_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", generated_shifty.name, generated_shifty)

	assert v3_process.returncode == 0, v3_process.stderr or v3_process.stdout
	assert v21_process.returncode == 0, v21_process.stderr or v21_process.stdout
	assert shifty_process.returncode == 0, shifty_process.stderr or shifty_process.stdout
	assert set(v3_outputs) == set(OUTPUT_SUFFIXES)
	_assert_output_sets_close(v21_outputs, v3_outputs)
	_assert_output_sets_close(shifty_outputs, v3_outputs)


def test_headerless_three_letter_shifty_is_detected_and_runs():
	source_shifty = ROOT / "inputs" / "bmr15476_shifty.txt"
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		headerless_shifty = temp_path / "test_shifty.txt"
		_make_headerless_three_letter_shifty(source_shifty, headerless_shifty)

		assert detect_supported_format(headerless_shifty) == "SHIFTY"

		process, outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", headerless_shifty.name, headerless_shifty)

	assert process.returncode == 0, process.stderr or process.stdout
	assert set(outputs) == set(OUTPUT_SUFFIXES)


def test_generated_shifty_writes_dense_numeric_rows_without_seq_comments():
	source_v3 = ROOT / "inputs" / "bmr15476_3.str.txt"
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		generated_shifty = temp_path / "generated_shifty.txt"

		convert_file(source_v3, "SHIFTY", generated_shifty)
		lines = generated_shifty.read_text().splitlines()

	assert lines[0] == "#NUM\tAA\tHA\tCA\tCB\tCO\tN\tHN"
	assert not any(line.startswith("#SEQ") for line in lines)
	assert len(lines) == 55
	assert lines[-1].split()[0] == "54"


def test_shifty_glycine_ha_roundtrip_duplicates_to_ha2_and_ha3():
	source_v3 = ROOT / "inputs" / "bmr15476_3.str.txt"
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		generated_shifty = temp_path / "generated_shifty.txt"
		roundtrip_v21 = temp_path / "roundtrip_v21.str.txt"

		convert_file(source_v3, "SHIFTY", generated_shifty)
		convert_file(generated_shifty, "NMRSTAR21", roundtrip_v21)

		lines = generated_shifty.read_text().splitlines()
		assert lines[0] == "#NUM\tAA\tHA\tCA\tCB\tCO\tN\tHN"

		gly20 = next(line for line in lines if line.split()[0] == "20")
		gly40 = next(line for line in lines if line.split()[0] == "40")
		for gly_line in [gly20, gly40]:
			parts = gly_line.split()
			ha = float(parts[2])
			assert ha != 0.0

		source_dataset = parse_supported_input(source_v3)
		roundtrip_dataset = parse_supported_input(roundtrip_v21)

	def gly_rows(dataset):
		return sorted(
			[
				(row.residue_number, row.residue_name, row.atom_name, row.shift_value)
				for row in dataset.shifts
				if row.residue_number in {20, 40} and row.atom_name in {"HA2", "HA3"}
			]
		)

	source_gly_rows = gly_rows(source_dataset)
	roundtrip_gly_rows = gly_rows(roundtrip_dataset)
	assert [row[:3] for row in source_gly_rows] == [row[:3] for row in roundtrip_gly_rows]
	expected_ha_by_residue = {}
	for gly_line in [gly20, gly40]:
		parts = gly_line.split()
		expected_ha_by_residue[int(parts[0])] = float(parts[2])
	for residue_number, _, atom_name, roundtrip_shift in roundtrip_gly_rows:
		assert atom_name in {"HA2", "HA3"}
		assert math.isclose(roundtrip_shift, expected_ha_by_residue[residue_number], rel_tol=0.0, abs_tol=1e-12)


def test_30786_lists_multiple_entities():
	entities = list_input_entities(ROOT / "inputs" / "30786_3.str.txt")
	assert [entity.entity_id for entity in entities] == ["1", "2"]
	assert entities[0].residue_count == 96
	assert entities[1].residue_count == 23


def test_30786_entity_specific_outputs_match_across_converted_formats():
	source_v3 = ROOT / "inputs" / "30786_3.str.txt"
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_path = Path(temp_dir)
		for entity_id in ["1", "2"]:
			generated_v21 = temp_path / f"30786_entity_{entity_id}.nmrstar21.str.txt"
			generated_shifty = temp_path / f"30786_entity_{entity_id}.shifty.txt"

			convert_file(source_v3, "NMRSTAR21", generated_v21, entity_id=entity_id)
			convert_file(source_v3, "SHIFTY", generated_shifty, entity_id=entity_id)

			v3_process, v3_outputs = run_script(
				"scripts/rci_v_1c_PyNMR-STAR.py",
				source_v3.name,
				source_v3,
				extra_args=["-entity", entity_id],
			)
			v21_process, v21_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", generated_v21.name, generated_v21)
			shifty_process, shifty_outputs = run_script("scripts/rci_v_1c_PyNMR-STAR.py", generated_shifty.name, generated_shifty)

			assert v3_process.returncode == 0, v3_process.stderr or v3_process.stdout
			assert v21_process.returncode == 0, v21_process.stderr or v21_process.stdout
			assert shifty_process.returncode == 0, shifty_process.stderr or shifty_process.stdout
			_assert_output_sets_close(v21_outputs, v3_outputs)
			_assert_output_sets_close(shifty_outputs, v3_outputs)
