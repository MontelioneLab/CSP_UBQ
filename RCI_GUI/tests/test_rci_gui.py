import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from rci_gui import (
    AdvancedOptions,
    build_cli_args,
    can_delete_run,
    delete_run_directory,
    discover_artifacts,
    read_text_file,
    safe_filename,
)


def test_safe_filename_drops_path_components():
    assert safe_filename("../../example.str") == "example.str"
    assert safe_filename("") == "input.str"


def test_discover_artifacts_collects_expected_outputs(tmp_path: Path):
    run_dir = tmp_path / "run"
    workspace = run_dir / "workspace"
    nested = workspace / "sample.str_txt__run" / "rcExample"
    nested.mkdir(parents=True)
    workspace.mkdir(parents=True, exist_ok=True)

    (workspace / "sample.str").write_text("data", encoding="utf-8")
    (workspace / "sample.str_RCI.gif").write_bytes(b"gif")
    (workspace / "sample.str_S2.gif").write_bytes(b"gif")
    (nested / "log").write_text("warning", encoding="utf-8")
    (workspace / "sample.str.RCI.txt").write_text("rci", encoding="utf-8")
    (run_dir / "stdout.log").write_text("stdout", encoding="utf-8")

    artifacts = discover_artifacts(run_dir, "sample.str")

    assert artifacts["images"] == [
        "workspace/sample.str_RCI.gif",
        "workspace/sample.str_S2.gif",
    ]
    assert "stdout.log" in artifacts["logs"]
    assert "workspace/sample.str_txt__run/rcExample/log" in artifacts["logs"]
    assert artifacts["text_outputs"] == ["workspace/sample.str.RCI.txt"]


def test_read_text_file_missing_path_returns_message(tmp_path: Path):
    message = read_text_file(tmp_path, "missing.log")
    assert message == "(file not found: missing.log)"


def test_build_cli_args_maps_advanced_options():
    options = AdvancedOptions(
        random_coil="Wang",
        neighbor_correction="Wang",
        exclude_unassigned=True,
        terminal_correction="Exclude first 3 and last 3 aa",
        fill_small_gaps=False,
        correct_referencing=True,
        predict_secondary_structure=False,
    )

    args, warnings = build_cli_args("sample.str", options)

    assert args == [
        "-b",
        "sample.str",
        "-mpl",
        "-r",
        "2",
        "-d",
        "0",
        "-no_i",
        "-no_fl6",
        "-nogapfill",
    ]
    assert warnings


def test_build_cli_args_default_options_use_new_rci_defaults():
    args, warnings = build_cli_args("sample.str", AdvancedOptions())

    assert args == [
        "-b",
        "sample.str",
        "-mpl",
        "-r",
        "4",
        "-d",
        "1",
        "-no_i",
        "-end_corr2",
        "-gapfill2",
        "-dynamr",
    ]
    assert warnings == []


def test_discover_artifacts_prefers_svg_outputs(tmp_path: Path):
    run_dir = tmp_path / "run"
    workspace = run_dir / "workspace"
    workspace.mkdir(parents=True)

    (workspace / "sample.str_RCI.gif").write_bytes(b"gif")
    (workspace / "sample.str_RCI.svg").write_text("<svg></svg>", encoding="utf-8")
    (workspace / "sample.str_MD_RMSD.svg").write_text("<svg></svg>", encoding="utf-8")

    artifacts = discover_artifacts(run_dir, "sample.str")

    assert artifacts["images"] == [
        "workspace/sample.str_MD_RMSD.svg",
        "workspace/sample.str_RCI.svg",
    ]


def test_can_delete_run_for_completed_and_failed():
    assert can_delete_run({"status": "completed"})
    assert can_delete_run({"status": "failed"})
    assert not can_delete_run({"status": "running"})


def test_delete_run_directory_removes_run(tmp_path: Path):
    run_dir = tmp_path / "run"
    (run_dir / "workspace").mkdir(parents=True)
    (run_dir / "run.json").write_text("{}", encoding="utf-8")

    delete_run_directory(run_dir)

    assert not run_dir.exists()
