from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

try:
    import pynmrstar
except ImportError:  # pragma: no cover - optional at import time
    pynmrstar = None


SHIFTY_HEADER_PREFIX = "#NUM"
SHIFTY_HEADER_SECOND_COLUMN = "AA"
SHIFTY_HEADER_COLUMNS = ["HA", "CA", "CB", "CO", "N", "HN"]
SHIFTY_LEGACY_COLUMN_INDEX = {"HA": 2, "CA": 3, "CB": 4, "CO": 5, "N": 6, "HN": 7}
SUPPORTED_ATOM_ORDER = {"H": 0, "HA": 1, "HA2": 2, "HA3": 3, "C": 4, "CA": 5, "CB": 6, "N": 7}

AA_ONE_TO_THREE = {
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
AA_THREE_TO_ONE = {value: key for key, value in AA_ONE_TO_THREE.items()}


@dataclass(frozen=True)
class ShiftRow:
    residue_number: int
    residue_name: str
    atom_name: str
    shift_value: float


@dataclass
class ShiftDataset:
    source_format: str
    sequence: list[tuple[int, str]]
    shifts: list[ShiftRow]
    entity_id: str | None = None


@dataclass(frozen=True)
class EntityInfo:
    entity_id: str
    residue_count: int
    residue_start: int
    residue_end: int


def _read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def _normalize_nmrstar_value(value):
    if value in (None, ".", "?", ""):
        return None
    return value


def _loop_tag_map(loop_obj) -> dict[str, int]:
    return {tag.split(".")[-1]: index for index, tag in enumerate(loop_obj.get_tag_names())}


def _normalize_entity_id(value) -> str | None:
    normalized = _normalize_nmrstar_value(value)
    if normalized is None:
        return None
    return str(normalized)


def _conversion_atom_name(residue_name: str, atom_name: str) -> str | None:
    upper_atom = atom_name.upper()
    if upper_atom in {"H", "HN", "NH"}:
        return "H"
    if upper_atom == "HA":
        return "HA"
    if residue_name == "GLY" and upper_atom in {"HA2", "HA3"}:
        return upper_atom
    if upper_atom in {"C", "CO"}:
        return "C"
    if upper_atom == "CA":
        return "CA"
    if upper_atom == "CB":
        return "CB"
    if upper_atom == "N":
        return "N"
    return None


def _normalize_conversion_shifts(source_format: str, sequence: dict[int, str], rows: list[ShiftRow]) -> ShiftDataset:
    grouped: dict[tuple[int, str, str], list[float]] = {}
    for row in rows:
        normalized_atom = _conversion_atom_name(row.residue_name, row.atom_name)
        if normalized_atom is None:
            continue
        key = (row.residue_number, row.residue_name, normalized_atom)
        grouped.setdefault(key, []).append(row.shift_value)

    normalized_rows: list[ShiftRow] = []
    for (residue_number, residue_name, atom_name), values in grouped.items():
        normalized_rows.append(
            ShiftRow(
                residue_number=residue_number,
                residue_name=residue_name,
                atom_name=atom_name,
                shift_value=sum(values) / len(values),
            )
        )

    normalized_rows.sort(
        key=lambda row: (
            row.residue_number,
            SUPPORTED_ATOM_ORDER.get(row.atom_name, 99),
            row.atom_name,
        )
    )
    return ShiftDataset(source_format, sorted(sequence.items()), normalized_rows)


def _parse_shifty_header(stripped_line: str) -> dict[str, int]:
    parts = stripped_line.split()
    if len(parts) < 3 or parts[0] != SHIFTY_HEADER_PREFIX or parts[1] != SHIFTY_HEADER_SECOND_COLUMN:
        raise ValueError("Not a SHIFTY header.")
    return {column_name: index for index, column_name in enumerate(parts)}


def _is_shifty_header_line(stripped_line: str) -> bool:
    try:
        _parse_shifty_header(stripped_line)
        return True
    except ValueError:
        return False


def _is_shifty_sequence_line(stripped_line: str) -> bool:
    parts = stripped_line.split()
    return len(parts) >= 3 and parts[0] == "#SEQ"


def _looks_like_shifty_row(stripped_line: str) -> bool:
    parts = stripped_line.split()
    if len(parts) < 8:
        return False
    try:
        int(parts[0])
    except ValueError:
        return False
    return (parts[1] in AA_ONE_TO_THREE) or (parts[1] in AA_THREE_TO_ONE)


def _shifty_residue_to_three_letter(token: str) -> str | None:
    token = token.upper()
    if token in AA_ONE_TO_THREE:
        return AA_ONE_TO_THREE[token]
    if token in AA_THREE_TO_ONE:
        return token
    return None


def _build_entity_infos(sequences_by_entity: dict[str, dict[int, str]]) -> list[EntityInfo]:
    entities: list[EntityInfo] = []
    for entity_id in sorted(sequences_by_entity, key=_entity_sort_key):
        residue_numbers = sorted(sequences_by_entity[entity_id])
        if not residue_numbers:
            continue
        entities.append(
            EntityInfo(
                entity_id=entity_id,
                residue_count=len(residue_numbers),
                residue_start=residue_numbers[0],
                residue_end=residue_numbers[-1],
            )
        )
    return entities


def _entity_sort_key(entity_id: str) -> tuple[int, object]:
    try:
        return (0, int(float(entity_id)))
    except ValueError:
        return (1, entity_id)


def _select_entity_id(sequences_by_entity: dict[str, dict[int, str]], requested_entity_id: str | None) -> str:
    if requested_entity_id is not None:
        if requested_entity_id not in sequences_by_entity:
            raise ValueError(f"Entity {requested_entity_id} was not found in the input file.")
        return requested_entity_id
    if not sequences_by_entity:
        raise ValueError("No usable entity sequences were found in the input file.")
    return sorted(sequences_by_entity, key=_entity_sort_key)[0]


def detect_supported_format(path: Path) -> str:
    for line in _read_lines(path):
        stripped = line.strip()
        if not stripped:
            continue
        if _is_shifty_header_line(stripped):
            return "SHIFTY"
        if _looks_like_shifty_row(stripped):
            return "SHIFTY"
        break
    if pynmrstar is not None:
        try:
            entry = pynmrstar.Entry.from_file(str(path))
            if entry.get_loops_by_category("_Atom_chem_shift"):
                return "NMRSTAR3"
        except Exception:
            pass
    return "NMRSTAR21"


def list_input_entities(path: Path) -> list[EntityInfo]:
    input_format = detect_supported_format(path)
    if input_format != "NMRSTAR3":
        dataset = parse_supported_input(path)
        if not dataset.sequence:
            return []
        residue_numbers = [residue_number for residue_number, _ in dataset.sequence]
        return [
            EntityInfo(
                entity_id=dataset.entity_id or "1",
                residue_count=len(residue_numbers),
                residue_start=min(residue_numbers),
                residue_end=max(residue_numbers),
            )
        ]

    if pynmrstar is None:
        raise RuntimeError("pynmrstar is required to inspect NMR-STAR 3 entities.")
    entry = pynmrstar.Entry.from_file(str(path))
    sequences_by_entity: dict[str, dict[int, str]] = {}
    for loop_obj in entry.get_loops_by_category("_Entity_poly_seq"):
        tag_map = _loop_tag_map(loop_obj)
        residue_tag = None
        for candidate_tag in ["Comp_index_ID", "Num", "Seq_ID"]:
            if candidate_tag in tag_map:
                residue_tag = candidate_tag
                break
        if residue_tag is None or "Mon_ID" not in tag_map:
            continue
        entity_id_tag = "Entity_ID" if "Entity_ID" in tag_map else None
        for row in loop_obj.data:
            residue_number = _normalize_nmrstar_value(row[tag_map[residue_tag]])
            residue_name = _normalize_nmrstar_value(row[tag_map["Mon_ID"]])
            if residue_number is None or residue_name is None:
                continue
            entity_id = _normalize_entity_id(row[tag_map[entity_id_tag]]) if entity_id_tag is not None else "1"
            if entity_id is None:
                entity_id = "1"
            sequences_by_entity.setdefault(entity_id, {})[int(float(residue_number))] = residue_name
    return _build_entity_infos(sequences_by_entity)


def parse_shifty(path: Path, entity_id: str | None = None) -> ShiftDataset:
    sequence: dict[int, str] = {}
    shifts: list[ShiftRow] = []
    header_found = False
    column_index: dict[str, int] | None = None
    for line in _read_lines(path):
        stripped = line.strip()
        if not stripped:
            continue
        if _is_shifty_header_line(stripped):
            header_found = True
            column_index = _parse_shifty_header(stripped)
            continue
        if _is_shifty_sequence_line(stripped):
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    residue_number = int(parts[1])
                except ValueError:
                    continue
                residue_name = parts[2].upper()
                if residue_name in AA_ONE_TO_THREE:
                    sequence[residue_number] = AA_ONE_TO_THREE[residue_name]
                elif residue_name in AA_THREE_TO_ONE:
                    sequence[residue_number] = residue_name
            continue
        if stripped.startswith("#"):
            continue
        if column_index is None and _looks_like_shifty_row(stripped):
            header_found = True
            column_index = dict(SHIFTY_LEGACY_COLUMN_INDEX)
        parts = stripped.split()
        if column_index is None:
            continue
        residue_number = int(parts[0])
        residue_name = _shifty_residue_to_three_letter(parts[1])
        if residue_name is None:
            continue
        sequence[residue_number] = residue_name
        ha_value = _parse_optional_shifty_float(parts, column_index, "HA")
        if ha_value not in (None, 0.0):
            if residue_name == "GLY":
                shifts.append(ShiftRow(residue_number, residue_name, "HA2", ha_value))
                shifts.append(ShiftRow(residue_number, residue_name, "HA3", ha_value))
            else:
                shifts.append(ShiftRow(residue_number, residue_name, "HA", ha_value))

        for column_name, atom_name in [("HN", "H"), ("CO", "C"), ("CA", "CA"), ("CB", "CB"), ("N", "N")]:
            shift_value = _parse_optional_shifty_float(parts, column_index, column_name)
            if shift_value in (None, 0.0):
                continue
            shifts.append(ShiftRow(residue_number, residue_name, atom_name, shift_value))
    if not header_found or not sequence:
        raise ValueError("No usable SHIFTY data were found.")
    dataset = _normalize_conversion_shifts("SHIFTY", sequence, shifts)
    dataset.entity_id = entity_id or "1"
    return dataset


def parse_nmrstar3(path: Path, entity_id: str | None = None) -> ShiftDataset:
    if pynmrstar is None:
        raise RuntimeError("pynmrstar is required to parse NMR-STAR 3 files.")
    entry = pynmrstar.Entry.from_file(str(path))

    sequences_by_entity: dict[str, dict[int, str]] = {}
    for loop_obj in entry.get_loops_by_category("_Entity_poly_seq"):
        tag_map = _loop_tag_map(loop_obj)
        residue_tag = None
        for candidate_tag in ["Comp_index_ID", "Num", "Seq_ID"]:
            if candidate_tag in tag_map:
                residue_tag = candidate_tag
                break
        if residue_tag is None or "Mon_ID" not in tag_map:
            continue
        entity_id_tag = "Entity_ID" if "Entity_ID" in tag_map else None
        for row in loop_obj.data:
            residue_number = _normalize_nmrstar_value(row[tag_map[residue_tag]])
            residue_name = _normalize_nmrstar_value(row[tag_map["Mon_ID"]])
            if residue_number is None or residue_name is None:
                continue
            row_entity_id = _normalize_entity_id(row[tag_map[entity_id_tag]]) if entity_id_tag is not None else "1"
            if row_entity_id is None:
                row_entity_id = "1"
            sequences_by_entity.setdefault(row_entity_id, {})[int(float(residue_number))] = residue_name

    selected_entity_id = _select_entity_id(sequences_by_entity, entity_id)
    sequence = sequences_by_entity[selected_entity_id]

    shifts: list[ShiftRow] = []
    for loop_obj in entry.get_loops_by_category("_Atom_chem_shift"):
        tag_map = _loop_tag_map(loop_obj)
        required_tags = ["Comp_index_ID", "Comp_ID", "Atom_ID", "Val"]
        if any(tag not in tag_map for tag in required_tags):
            continue
        entity_id_tag = "Entity_ID" if "Entity_ID" in tag_map else None
        for row in loop_obj.data:
            row_entity_id = _normalize_entity_id(row[tag_map[entity_id_tag]]) if entity_id_tag is not None else "1"
            if row_entity_id is None:
                row_entity_id = "1"
            if row_entity_id != selected_entity_id:
                continue
            residue_number = _normalize_nmrstar_value(row[tag_map["Comp_index_ID"]])
            residue_name = _normalize_nmrstar_value(row[tag_map["Comp_ID"]])
            atom_name = _normalize_nmrstar_value(row[tag_map["Atom_ID"]])
            shift_value = _normalize_nmrstar_value(row[tag_map["Val"]])
            if None in (residue_number, residue_name, atom_name, shift_value):
                continue
            shifts.append(
                ShiftRow(int(float(residue_number)), residue_name, atom_name, float(shift_value))
            )
    if not sequence or not shifts:
        raise ValueError("No usable NMR-STAR 3 data were found.")
    dataset = _normalize_conversion_shifts("NMRSTAR3", sequence, shifts)
    dataset.entity_id = selected_entity_id
    return dataset


def _parse_nmrstar21_modern_loops(path: Path) -> tuple[dict[int, str], list[ShiftRow]] | None:
    """Parse loop-style NMR-STAR 2.1 files that use _Chem_shift_value columns."""
    sequence: dict[int, str] = {}
    shifts: list[ShiftRow] = []
    lines = list(_read_lines(path))
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if not parts or parts[0] != "loop_":
            index += 1
            continue

        tags: list[str] = []
        index += 1
        while index < len(lines):
            row_parts = lines[index].split()
            if not row_parts:
                index += 1
                continue
            if row_parts[0].startswith("_"):
                tags.append(row_parts[0])
                index += 1
                continue
            break

        tag_index = {tag: position for position, tag in enumerate(tags)}
        if "_Residue_seq_code" in tag_index and "_Residue_label" in tag_index and "_Chem_shift_value" not in tag_index:
            while index < len(lines):
                data_parts = lines[index].split()
                if not data_parts or data_parts[0] == "stop_":
                    break
                offset = 0
                while offset + 1 < len(data_parts):
                    try:
                        residue_number = int(data_parts[offset])
                        residue_name = data_parts[offset + 1]
                    except ValueError:
                        break
                    if residue_name in AA_THREE_TO_ONE:
                        sequence[residue_number] = residue_name
                    offset += 2
                index += 1
            continue

        if "_Chem_shift_value" in tag_index:
            residue_number_idx = tag_index.get("_Residue_seq_code")
            if residue_number_idx is None:
                residue_number_idx = tag_index.get("_Residue_author_seq_code")
            label_idx = tag_index.get("_Residue_label")
            atom_idx = tag_index.get("_Atom_name")
            value_idx = tag_index["_Chem_shift_value"]
            while index < len(lines):
                data_parts = lines[index].split()
                if not data_parts or data_parts[0] == "stop_":
                    break
                if len(data_parts) <= value_idx:
                    index += 1
                    continue
                try:
                    if residue_number_idx is not None:
                        residue_number = int(float(data_parts[residue_number_idx]))
                    else:
                        residue_number = int(float(data_parts[0]))
                    residue_name = data_parts[label_idx] if label_idx is not None else ""
                    atom_name = data_parts[atom_idx] if atom_idx is not None else ""
                    shift_value = float(data_parts[value_idx])
                except (ValueError, IndexError):
                    index += 1
                    continue
                if residue_name in AA_THREE_TO_ONE:
                    shifts.append(ShiftRow(residue_number, residue_name, atom_name, shift_value))
                index += 1
            continue

        index += 1

    if not shifts:
        return None
    return sequence, shifts


def parse_nmrstar21(path: Path, entity_id: str | None = None) -> ShiftDataset:
    modern = _parse_nmrstar21_modern_loops(path)
    if modern is not None:
        sequence, shifts = modern
        if not sequence:
            for row in shifts:
                sequence.setdefault(row.residue_number, row.residue_name)
        dataset = _normalize_conversion_shifts("NMRSTAR21", sequence, shifts)
        dataset.entity_id = entity_id or "1"
        return dataset

    sequence: dict[int, str] = {}
    shifts: list[ShiftRow] = []
    in_sequence = False
    in_shifts = False
    for line in _read_lines(path):
        parts = line.split()
        if not parts:
            continue
        token = parts[0]
        if token == "_Chem_shift_ambiguity_code":
            in_shifts = True
            in_sequence = False
            continue
        if token == "_Residue_label" and not in_shifts:
            in_sequence = True
            in_shifts = False
            continue
        if token == "stop_":
            in_sequence = False
            in_shifts = False
            continue

        if in_sequence and len(parts) >= 2:
            try:
                residue_number = int(parts[0])
            except ValueError:
                continue
            residue_name = next((entry for entry in parts[1:] if len(entry) == 3 and entry in AA_THREE_TO_ONE), None)
            if residue_name is not None:
                sequence[residue_number] = residue_name

        if in_shifts and len(parts) >= 6:
            try:
                residue_number = int(parts[1])
                residue_name = parts[2]
                atom_name = parts[3]
                shift_value = float(parts[5])
            except ValueError:
                continue
            if residue_name in AA_THREE_TO_ONE:
                shifts.append(ShiftRow(residue_number, residue_name, atom_name, shift_value))
    if not sequence:
        for row in shifts:
            sequence.setdefault(row.residue_number, row.residue_name)
    if not sequence or not shifts:
        raise ValueError("No usable NMR-STAR 2.1 data were found.")
    dataset = _normalize_conversion_shifts("NMRSTAR21", sequence, shifts)
    dataset.entity_id = entity_id or "1"
    return dataset


def parse_supported_input(path: Path, entity_id: str | None = None) -> ShiftDataset:
    input_format = detect_supported_format(path)
    if input_format == "SHIFTY":
        return parse_shifty(path, entity_id)
    if input_format == "NMRSTAR3":
        return parse_nmrstar3(path, entity_id)
    return parse_nmrstar21(path, entity_id)


def write_shifty(dataset: ShiftDataset, path: Path) -> None:
    shift_lookup = {(row.residue_number, row.atom_name): row.shift_value for row in dataset.shifts}
    lines = ["\t".join([SHIFTY_HEADER_PREFIX, SHIFTY_HEADER_SECOND_COLUMN, *SHIFTY_HEADER_COLUMNS])]
    for residue_number, residue_name in dataset.sequence:
        residue_code = AA_THREE_TO_ONE[residue_name]
        ha_value = shift_lookup.get((residue_number, "HA"))
        if residue_name == "GLY":
            ha_components = [
                value
                for value in [shift_lookup.get((residue_number, "HA2")), shift_lookup.get((residue_number, "HA3"))]
                if value is not None
            ]
            if ha_components:
                ha_value = sum(ha_components) / len(ha_components)
        row = [
            str(residue_number),
            residue_code,
            _format_shifty_value(ha_value),
            _format_shifty_value(shift_lookup.get((residue_number, "CA"))),
            _format_shifty_value(shift_lookup.get((residue_number, "CB"))),
            _format_shifty_value(shift_lookup.get((residue_number, "C"))),
            _format_shifty_value(shift_lookup.get((residue_number, "N"))),
            _format_shifty_value(shift_lookup.get((residue_number, "H"))),
        ]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_nmrstar3(dataset: ShiftDataset, path: Path) -> None:
    lines = [
        "data_converted",
        "",
        "save_entry_information",
        "   _Entry.Sf_category entry_information",
        "   _Entry.Sf_framecode entry_information",
        "   _Entry.ID 1",
        "   _Entry.Title",
        ";",
        "Converted chemical shift file",
        ";",
        "save_",
        "",
        "save_sequence",
        "   _Saveframe_category entity",
        "   _Saveframe_framecode sequence",
        "   loop_",
        "      _Entity_poly_seq.Comp_index_ID",
        "      _Entity_poly_seq.Mon_ID",
    ]
    for residue_number, residue_name in dataset.sequence:
        lines.append(f"      {residue_number} {residue_name}")
    lines.extend(
        [
            "   stop_",
            "save_",
            "",
            "save_chemical_shifts",
            "   _Saveframe_category assigned_chemical_shifts",
            "   _Saveframe_framecode chemical_shifts",
            "   loop_",
            "      _Atom_chem_shift.Comp_index_ID",
            "      _Atom_chem_shift.Comp_ID",
            "      _Atom_chem_shift.Atom_ID",
            "      _Atom_chem_shift.Val",
        ]
    )
    for row in dataset.shifts:
        lines.append(
            f"      {row.residue_number} {row.residue_name} {row.atom_name} {_format_float(row.shift_value)}"
        )
    lines.extend(["   stop_", "save_", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_nmrstar21(dataset: ShiftDataset, path: Path) -> None:
    lines = [
        "data_converted",
        "",
        "save_entry_information",
        "   _NMR_STAR_version 2.1.1",
        "save_",
        "",
        "loop_",
        "   _Residue_label",
    ]
    for residue_number, residue_name in dataset.sequence:
        lines.append(f"   {residue_number} {residue_name}")
    lines.extend(
        [
            "stop_",
            "",
            "loop_",
            "   _Chem_shift_ambiguity_code",
            "   _Residue_seq_code",
            "   _Residue_label",
            "   _Atom_name",
            "   _Atom_type",
            "   _Chem_shift_value",
        ]
    )
    for index, row in enumerate(dataset.shifts, start=1):
        atom_type = row.atom_name[0]
        lines.append(
            f"   {index} {row.residue_number} {row.residue_name} {row.atom_name} {atom_type} {_format_float(row.shift_value)}"
        )
    lines.extend(["stop_", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def convert_dataset(dataset: ShiftDataset, target_format: str, output_path: Path) -> None:
    if target_format == "SHIFTY":
        write_shifty(dataset, output_path)
    elif target_format == "NMRSTAR3":
        write_nmrstar3(dataset, output_path)
    elif target_format == "NMRSTAR21":
        write_nmrstar21(dataset, output_path)
    else:
        raise ValueError(f"Unsupported target format: {target_format}")


def convert_file(input_path: Path, target_format: str, output_path: Path, entity_id: str | None = None) -> ShiftDataset:
    dataset = parse_supported_input(input_path, entity_id)
    convert_dataset(dataset, target_format, output_path)
    return dataset


def suggested_output_name(input_path: Path, target_format: str) -> str:
    base_name = input_path.stem
    if target_format == "SHIFTY":
        return f"{base_name}.shifty.txt"
    if target_format == "NMRSTAR3":
        return f"{base_name}.nmrstar3.str.txt"
    if target_format == "NMRSTAR21":
        return f"{base_name}.nmrstar21.str.txt"
    raise ValueError(f"Unsupported target format: {target_format}")


def _format_shifty_value(value: float | None) -> str:
    if value is None:
        return "0.00"
    return _format_float(value)


def _format_float(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _parse_optional_shifty_float(parts: list[str], column_index: dict[str, int], column_name: str) -> float | None:
    index = column_index.get(column_name)
    if index is None or index >= len(parts):
        return None
    try:
        return float(parts[index])
    except ValueError:
        return None
