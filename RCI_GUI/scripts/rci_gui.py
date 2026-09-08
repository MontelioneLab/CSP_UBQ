import json
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from format_conversion import convert_file, detect_supported_format, list_input_entities, suggested_output_name


SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
SCRIPT_PATH = SCRIPTS_DIR / "rci_v_1c_PyNMR-STAR.py"
RUNS_DIR = PROJECT_ROOT / "outputs" / "gui_runs"
METADATA_NAME = "run.json"
IMAGE_EXTENSIONS = {".gif", ".jpg", ".jpeg", ".png", ".svg"}
TEXT_OUTPUT_SUFFIXES = (".RCI.txt", ".S2.txt", ".MD_RMSD.txt", ".NMR_RMSD.txt")


@dataclass
class UploadedFile:
    filename: str
    content: bytes


@dataclass
class AdvancedOptions:
    random_coil: str = "Schwarzinger"
    neighbor_correction: str = "Schwarzinger"
    exclude_unassigned: bool = True
    terminal_correction: str = "end_corr2"
    fill_small_gaps: bool = True
    correct_referencing: bool = False
    predict_secondary_structure: bool = True


RANDOM_COIL_FLAGS = {
    "Wishart": "1",
    "Wang": "2",
    "Lukin": "3",
    "Schwarzinger": "4",
}

NEIGHBOR_FLAGS = {
    "Wang": "0",
    "Schwarzinger": "1",
}


def safe_filename(filename: str) -> str:
    candidate = Path(filename).name.strip()
    if not candidate:
        return "input.str"
    return candidate.replace("\x00", "")


def make_run_id() -> str:
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return f"{timestamp}-{uuid.uuid4().hex[:8]}"


def load_metadata(run_dir: Path) -> dict:
    metadata_path = run_dir / METADATA_NAME
    if not metadata_path.exists():
        return {}
    return json.loads(metadata_path.read_text(encoding="utf-8"))


def save_metadata(run_dir: Path, payload: dict) -> None:
    (run_dir / METADATA_NAME).write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def discover_artifacts(run_dir: Path, input_name: str) -> dict:
    images = []
    logs = []
    text_outputs = []
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        relative = path.relative_to(run_dir).as_posix()
        if suffix in IMAGE_EXTENSIONS:
            images.append(relative)
        elif path.name in {"stdout.log", "stderr.log", "command.log", "log"}:
            logs.append(relative)
        elif path.name.endswith(TEXT_OUTPUT_SUFFIXES) or path.name == "coil_neightbor_corr":
            text_outputs.append(relative)

    input_relative = f"workspace/{input_name}"
    text_outputs = [path for path in text_outputs if path != input_relative]
    svg_images = [path for path in images if path.lower().endswith(".svg")]
    if svg_images:
        images = svg_images

    return {
        "images": images,
        "logs": logs,
        "text_outputs": text_outputs,
    }


def advanced_options_to_dict(options: AdvancedOptions) -> dict[str, Any]:
    return {
        "random_coil": options.random_coil,
        "neighbor_correction": options.neighbor_correction,
        "exclude_unassigned": options.exclude_unassigned,
        "terminal_correction": options.terminal_correction,
        "fill_small_gaps": options.fill_small_gaps,
        "correct_referencing": options.correct_referencing,
        "predict_secondary_structure": options.predict_secondary_structure,
    }


def build_cli_args(
    input_name: str,
    options: AdvancedOptions | None = None,
    entity_id: str | None = None,
) -> tuple[list[str], list[str]]:
    args = ["-b", input_name, "-mpl"]
    if entity_id:
        args.extend(["-entity", entity_id])
    warnings: list[str] = []
    if options is None:
        return args, warnings

    args.extend(["-r", RANDOM_COIL_FLAGS[options.random_coil]])
    args.extend(["-d", NEIGHBOR_FLAGS[options.neighbor_correction]])

    if options.exclude_unassigned:
        args.append("-no_i")

    if options.terminal_correction in {"end_corr2", "Yes"}:
        args.append("-end_corr2")
    elif options.terminal_correction in {"end_corr3"}:
        args.append("-end_corr3")
    elif options.terminal_correction in {"end_corr5"}:
        args.append("-end_corr5")
    elif options.terminal_correction in {"end_corr0", "No"}:
        args.append("-end_corr0")
    elif options.terminal_correction == "Exclude first 3 and last 3 aa":
        args.append("-no_fl6")

    args.append("-gapfill2" if options.fill_small_gaps else "-nogapfill")

    if options.correct_referencing:
        warnings.append(
            "Chemical shift re-referencing is not exposed by the current script. "
            "The GUI kept the script default for this setting."
        )

    if options.predict_secondary_structure:
        args.append("-dynamr")

    return args, warnings


def create_run(
    input_path: Path,
    options: AdvancedOptions | None = None,
    entity_id: str | None = None,
) -> Path:
    if options is None:
        options = AdvancedOptions()
    RUNS_DIR.mkdir(exist_ok=True)
    run_id = make_run_id()
    run_dir = RUNS_DIR / run_id
    work_dir = run_dir / "workspace"
    work_dir.mkdir(parents=True, exist_ok=True)
    copied_name = safe_filename(input_path.name)
    shutil.copy2(input_path, work_dir / copied_name)
    cli_args, warnings = build_cli_args(copied_name, options, entity_id)
    save_metadata(
        run_dir,
        {
            "advanced_options": advanced_options_to_dict(options),
            "command_args": cli_args,
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "entity_id": entity_id,
            "finished_at": "",
            "images": [],
            "input_name": copied_name,
            "logs": [],
            "return_code": None,
            "status": "running",
            "text_outputs": [],
            "warnings": warnings,
        },
    )
    return run_dir


def run_rci_job(
    run_dir: Path,
    process_registry: list[subprocess.Popen[str]] | None = None,
) -> dict:
    metadata = load_metadata(run_dir)
    input_name = metadata["input_name"]
    work_dir = run_dir / "workspace"
    command = [sys.executable, str(SCRIPT_PATH), *metadata.get("command_args", ["-b", input_name, "-mpl"])]
    (run_dir / "command.log").write_text(" ".join(command), encoding="utf-8")

    process = subprocess.Popen(
        command,
        cwd=work_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process_registry is not None:
        process_registry[:] = [process]
    try:
        stdout, stderr = process.communicate()
    finally:
        if process_registry is not None:
            process_registry.clear()

    (run_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(stderr, encoding="utf-8")

    if process.returncode is None:
        status = "cancelled"
    elif process.returncode == 0:
        status = "completed"
    else:
        status = "failed"

    metadata.update(
        {
            "finished_at": datetime.now().isoformat(timespec="seconds"),
            "return_code": process.returncode,
            "status": status,
        }
    )
    metadata.update(discover_artifacts(run_dir, input_name))
    save_metadata(run_dir, metadata)
    return metadata


def read_text_file(base_dir: Path, relative_path: str) -> str:
    path = base_dir / relative_path
    if not path.is_file():
        return f"(file not found: {relative_path})"
    return path.read_text(encoding="utf-8", errors="replace")


def can_delete_run(metadata: dict) -> bool:
    return metadata.get("status") in {"completed", "failed"}


def delete_run_directory(run_dir: Path) -> None:
    shutil.rmtree(run_dir)


QT_IMPORT_ERROR = None
try:
    from PySide6.QtCore import QObject, QThread, Qt, Signal
    from PySide6.QtGui import QAction, QMovie, QPixmap
    from PySide6.QtSvgWidgets import QSvgWidget
    from PySide6.QtWidgets import (
        QApplication,
        QComboBox,
        QFileDialog,
        QHBoxLayout,
        QLabel,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QPlainTextEdit,
        QSplitter,
        QStatusBar,
        QStackedLayout,
        QTabWidget,
        QToolBar,
        QVBoxLayout,
        QWidget,
    )
except ImportError as exc:
    QT_IMPORT_ERROR = exc


if QT_IMPORT_ERROR is None:
    class RunWorker(QObject):
        finished = Signal(str, dict)
        failed = Signal(str, str)

        def __init__(self, run_dir: Path) -> None:
            super().__init__()
            self.run_dir = run_dir
            self._process_registry: list[subprocess.Popen[str]] = []

        def cancel(self) -> None:
            if not self._process_registry:
                return
            process = self._process_registry[0]
            if process.poll() is not None:
                return
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()

        def run(self) -> None:
            try:
                metadata = run_rci_job(self.run_dir, self._process_registry)
            except Exception as exc:  # pragma: no cover - UI path
                self.failed.emit(self.run_dir.name, str(exc))
                return
            self.finished.emit(self.run_dir.name, metadata)


    class ImageViewer(QWidget):
        def __init__(self) -> None:
            super().__init__()
            self.setMinimumSize(420, 320)
            self.setStyleSheet("background:#f4f1ea; border:1px solid #d8d0c2; border-radius:12px;")
            self.current_path: Path | None = None
            self.movie: QMovie | None = None
            self.pixmap_source: QPixmap | None = None
            self.stack = QStackedLayout(self)
            self.stack.setContentsMargins(8, 8, 8, 8)

            self.message_label = QLabel("Run an analysis to preview output images.")
            self.message_label.setAlignment(Qt.AlignCenter)
            self.message_label.setWordWrap(True)

            self.raster_label = QLabel()
            self.raster_label.setAlignment(Qt.AlignCenter)

            self.svg_widget = QSvgWidget()

            self.stack.addWidget(self.message_label)
            self.stack.addWidget(self.raster_label)
            self.stack.addWidget(self.svg_widget)
            self.stack.setCurrentWidget(self.message_label)

        def show_image(self, image_path: Path) -> None:
            self.current_path = image_path
            self.movie = None
            self.pixmap_source = None
            lower_suffix = image_path.suffix.lower()
            if lower_suffix == ".svg":
                self.svg_widget.load(str(image_path))
                self.stack.setCurrentWidget(self.svg_widget)
            elif lower_suffix == ".gif":
                self.movie = QMovie(str(image_path))
                self.movie.setScaledSize(self._fit_size())
                self.raster_label.setMovie(self.movie)
                self.stack.setCurrentWidget(self.raster_label)
                self.movie.start()
            else:
                self.pixmap_source = QPixmap(str(image_path))
                if self.pixmap_source.isNull():
                    self.clear_view(f"Unable to load image:\n{image_path.name}")
                    return
                self.stack.setCurrentWidget(self.raster_label)
                self._refresh_pixmap()

        def clear_view(self, message: str) -> None:
            if self.movie is not None:
                self.movie.stop()
            self.movie = None
            self.pixmap_source = None
            self.current_path = None
            self.raster_label.clear()
            self.message_label.setText(message)
            self.stack.setCurrentWidget(self.message_label)

        def resizeEvent(self, event) -> None:  # pragma: no cover - UI path
            super().resizeEvent(event)
            if self.movie is not None:
                self.movie.setScaledSize(self._fit_size())
            elif self.pixmap_source is not None:
                self._refresh_pixmap()

        def _refresh_pixmap(self) -> None:
            if self.pixmap_source is None:
                return
            scaled = self.pixmap_source.scaled(
                self._fit_size(),
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation,
            )
            self.raster_label.setPixmap(scaled)

        def _fit_size(self):
            return self.rect().adjusted(16, 16, -16, -16).size()


    class RciMainWindow(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("RCI Qt Wrapper")
            self.resize(1340, 860)

            self.selected_input: Path | None = None
            self.conversion_input: Path | None = None
            self.selected_entity_id: str | None = None
            self.conversion_entity_id: str | None = None
            self.current_run_dir: Path | None = None
            self.current_metadata: dict = {}
            self.current_image_index = 0
            self.worker: RunWorker | None = None
            self.worker_thread: QThread | None = None
            self._exiting = False

            self._build_ui()
            self._load_runs()

        def _build_ui(self) -> None:
            central = QWidget()
            self.setCentralWidget(central)
            root_layout = QHBoxLayout(central)
            self.main_tabs = QTabWidget()
            root_layout.addWidget(self.main_tabs)

            analysis_tab = QWidget()
            analysis_layout = QHBoxLayout(analysis_tab)
            splitter = QSplitter()
            analysis_layout.addWidget(splitter)

            left_panel = QWidget()
            left_layout = QVBoxLayout(left_panel)

            self.file_label = QLabel("Input file: none selected")
            self.file_label.setWordWrap(True)
            self.status_label = QLabel("Status: idle")
            self.run_list = QListWidget()
            self.run_list.currentItemChanged.connect(self._on_run_selected)
            self.entity_label = QLabel("Entity: default")
            self.analysis_entity_combo = QComboBox()

            pick_button = QPushButton("Choose Input File")
            pick_button.clicked.connect(self._choose_input_file)
            self.run_button = QPushButton("Run Script")
            self.run_button.clicked.connect(self._start_run)
            self.delete_run_button = QPushButton("Delete Completed Run Output")
            self.delete_run_button.clicked.connect(self._delete_selected_run)

            left_layout.addWidget(self.file_label)
            left_layout.addWidget(self.entity_label)
            left_layout.addWidget(self.analysis_entity_combo)
            left_layout.addWidget(self.status_label)
            left_layout.addWidget(pick_button)
            left_layout.addWidget(self.run_button)
            left_layout.addWidget(self.delete_run_button)
            left_layout.addWidget(QLabel("Previous runs"))
            left_layout.addWidget(self.run_list, 1)

            right_panel = QWidget()
            right_layout = QVBoxLayout(right_panel)

            image_nav = QHBoxLayout()
            self.prev_button = QPushButton("Previous Image")
            self.prev_button.clicked.connect(lambda: self._step_image(-1))
            self.next_button = QPushButton("Next Image")
            self.next_button.clicked.connect(lambda: self._step_image(1))
            self.image_caption = QLabel("No images loaded")
            self.image_caption.setWordWrap(True)
            image_nav.addWidget(self.prev_button)
            image_nav.addWidget(self.next_button)
            image_nav.addWidget(self.image_caption, 1)

            self.image_viewer = ImageViewer()

            self.tabs = QTabWidget()
            self.logs_list = QListWidget()
            self.logs_list.currentItemChanged.connect(self._on_log_selected)
            self.log_text = QPlainTextEdit()
            self.log_text.setReadOnly(True)

            logs_tab = QWidget()
            logs_layout = QHBoxLayout(logs_tab)
            logs_layout.addWidget(self.logs_list, 1)
            logs_layout.addWidget(self.log_text, 3)

            self.outputs_list = QListWidget()
            self.outputs_list.currentItemChanged.connect(self._on_output_selected)
            self.output_text = QPlainTextEdit()
            self.output_text.setReadOnly(True)

            outputs_tab = QWidget()
            outputs_layout = QHBoxLayout(outputs_tab)
            outputs_layout.addWidget(self.outputs_list, 1)
            outputs_layout.addWidget(self.output_text, 3)

            self.tabs.addTab(logs_tab, "Logs")
            self.tabs.addTab(outputs_tab, "Text Outputs")

            right_layout.addLayout(image_nav)
            right_layout.addWidget(self.image_viewer, 3)
            right_layout.addWidget(self.tabs, 2)

            splitter.addWidget(left_panel)
            splitter.addWidget(right_panel)
            splitter.setSizes([320, 1020])

            self.main_tabs.addTab(analysis_tab, "Run Analysis")
            self.main_tabs.addTab(self._build_conversion_tab(), "Convert Formats")

            toolbar = QToolBar("Actions")
            toolbar.setMovable(False)
            self.addToolBar(toolbar)
            refresh_action = QAction("Refresh Runs", self)
            refresh_action.triggered.connect(self._load_runs)
            toolbar.addAction(refresh_action)

            exit_action = QAction("Exit", self)
            exit_action.triggered.connect(self._request_exit)
            toolbar.addAction(exit_action)

            self.setStatusBar(QStatusBar())
            self.delete_run_button.setEnabled(False)
            self._update_nav_buttons()

        def _build_conversion_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            self.convert_input_label = QLabel("Source file: none selected")
            self.convert_input_label.setWordWrap(True)
            self.convert_detected_label = QLabel("Detected format: unknown")
            self.convert_detected_label.setWordWrap(True)
            self.convert_entity_label = QLabel("Entity: default")
            self.convert_entity_combo = QComboBox()

            choose_input_button = QPushButton("Choose Source File")
            choose_input_button.clicked.connect(self._choose_conversion_input)

            self.convert_target_combo = QComboBox()
            self.convert_target_combo.addItem("NMR-STAR 3", "NMRSTAR3")
            self.convert_target_combo.addItem("NMR-STAR 2.1", "NMRSTAR21")
            self.convert_target_combo.addItem("SHIFTY", "SHIFTY")

            convert_button = QPushButton("Convert and Save As...")
            convert_button.clicked.connect(self._run_conversion)

            self.convert_status = QPlainTextEdit()
            self.convert_status.setReadOnly(True)
            self.convert_status.setPlainText(
                "Choose a supported input file, select a target format, and save the converted output."
            )

            layout.addWidget(self.convert_input_label)
            layout.addWidget(self.convert_detected_label)
            layout.addWidget(self.convert_entity_label)
            layout.addWidget(self.convert_entity_combo)
            layout.addWidget(choose_input_button)
            layout.addWidget(QLabel("Target format"))
            layout.addWidget(self.convert_target_combo)
            layout.addWidget(convert_button)
            layout.addWidget(self.convert_status, 1)
            return tab

        def _choose_input_file(self) -> None:
            start_dir = str(PROJECT_ROOT / "inputs")
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Choose NMR-STAR Input File",
                start_dir,
                "NMR-STAR files (*.str *.txt);;All files (*)",
            )
            if not file_path:
                return
            self.selected_input = Path(file_path)
            self.file_label.setText(f"Input file: {self.selected_input}")
            self._populate_entity_combo(
                self.analysis_entity_combo,
                self.selected_input,
                lambda entity_id: setattr(self, "selected_entity_id", entity_id),
                self.entity_label,
            )

        def _choose_conversion_input(self) -> None:
            start_dir = str(PROJECT_ROOT / "inputs")
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Choose Source File for Conversion",
                start_dir,
                "Supported files (*.str *.txt);;All files (*)",
            )
            if not file_path:
                return
            self.conversion_input = Path(file_path)
            detected = detect_supported_format(self.conversion_input)
            self.convert_input_label.setText(f"Source file: {self.conversion_input}")
            self.convert_detected_label.setText(f"Detected format: {detected}")
            self._populate_entity_combo(
                self.convert_entity_combo,
                self.conversion_input,
                lambda entity_id: setattr(self, "conversion_entity_id", entity_id),
                self.convert_entity_label,
            )

        def _populate_entity_combo(
            self,
            combo: QComboBox,
            source_path: Path,
            setter,
            label: QLabel,
        ) -> None:
            combo.blockSignals(True)
            combo.clear()
            entities = list_input_entities(source_path)
            for entity in entities:
                combo.addItem(
                    f"Entity {entity.entity_id} ({entity.residue_start}-{entity.residue_end}, {entity.residue_count} residues)",
                    entity.entity_id,
                )
            combo.blockSignals(False)

            if combo.count() == 0:
                setter(None)
                combo.setEnabled(False)
                label.setText("Entity: default")
                return

            combo.setEnabled(combo.count() > 1)
            selected_entity_id = combo.itemData(0)
            setter(selected_entity_id)
            label.setText(f"Entity: {selected_entity_id}")

            def handle_index_change(index: int) -> None:
                entity_id = combo.itemData(index) if index >= 0 else None
                setter(entity_id)
                label.setText(f"Entity: {entity_id}" if entity_id else "Entity: default")

            try:
                combo.currentIndexChanged.disconnect()
            except Exception:
                pass
            combo.currentIndexChanged.connect(handle_index_change)

        def _run_conversion(self) -> None:
            if self.conversion_input is None:
                QMessageBox.warning(self, "Missing Source File", "Choose a source file to convert first.")
                return
            if not self.conversion_input.exists():
                QMessageBox.warning(self, "Missing File", "The selected source file no longer exists.")
                return

            source_format = detect_supported_format(self.conversion_input)
            target_format = self.convert_target_combo.currentData()
            if source_format == target_format:
                QMessageBox.information(
                    self,
                    "Choose Another Format",
                    "Select a target format that differs from the detected input format.",
                )
                return

            suggested_name = suggested_output_name(self.conversion_input, target_format)
            output_path, _ = QFileDialog.getSaveFileName(
                self,
                "Save Converted File",
                str(PROJECT_ROOT / "inputs" / suggested_name),
                "Supported files (*.str *.txt);;All files (*)",
            )
            if not output_path:
                return

            output_file = Path(output_path)
            convert_file(self.conversion_input, target_format, output_file, self.conversion_entity_id)
            self.convert_status.setPlainText(
                "\n".join(
                    [
                        f"Converted: {self.conversion_input}",
                        f"Detected format: {source_format}",
                        f"Entity: {self.conversion_entity_id or 'default'}",
                        f"Target format: {target_format}",
                        f"Saved to: {output_file}",
                    ]
                )
            )
            self.statusBar().showMessage(f"Saved converted file to {output_file.name}.")

        def _load_runs(self) -> None:
            RUNS_DIR.mkdir(exist_ok=True)
            current_id = self.current_run_dir.name if self.current_run_dir else None
            self.run_list.clear()
            selected_row = None
            run_dirs = sorted([path for path in RUNS_DIR.iterdir() if path.is_dir()], reverse=True)
            for index, run_dir in enumerate(run_dirs):
                metadata = load_metadata(run_dir)
                item = QListWidgetItem(
                    f"{run_dir.name} [{metadata.get('status', 'unknown')}] {metadata.get('input_name', '')}"
                )
                item.setData(Qt.UserRole, str(run_dir))
                self.run_list.addItem(item)
                if run_dir.name == current_id:
                    selected_row = index

            if self.run_list.count() == 0:
                self.current_run_dir = None
                self.current_metadata = {}
                self.image_viewer.clear_view("No runs yet.")
                self.logs_list.clear()
                self.outputs_list.clear()
                self.log_text.clear()
                self.output_text.clear()
                self.image_caption.setText("No images loaded")
                self.delete_run_button.setEnabled(False)
                self._update_nav_buttons()
                return

            if selected_row is None:
                selected_row = 0
            self.run_list.setCurrentRow(selected_row)

        def _start_run(self) -> None:
            if self.selected_input is None:
                QMessageBox.warning(self, "Missing Input", "Choose an input NMR-STAR file first.")
                return
            if not self.selected_input.exists():
                QMessageBox.warning(self, "Missing File", "The selected input file no longer exists.")
                return
            if self.worker_thread is not None:
                QMessageBox.information(self, "Run In Progress", "Wait for the current run to finish.")
                return

            run_dir = create_run(self.selected_input, options=AdvancedOptions(), entity_id=self.selected_entity_id)
            self.current_run_dir = run_dir
            self.current_metadata = load_metadata(run_dir)
            self.status_label.setText(f"Status: running {run_dir.name}")
            self.statusBar().showMessage("Running script...")
            self.run_button.setEnabled(False)
            self._load_runs()

            self.worker_thread = QThread(self)
            self.worker = RunWorker(run_dir)
            self.worker.moveToThread(self.worker_thread)
            self.worker_thread.started.connect(self.worker.run)
            self.worker.finished.connect(self._on_run_finished)
            self.worker.failed.connect(self._on_run_failed)
            self.worker.finished.connect(self.worker_thread.quit)
            self.worker.failed.connect(self.worker_thread.quit)
            self.worker_thread.finished.connect(self._cleanup_worker)
            self.worker_thread.start()

        def _on_run_finished(self, run_id: str, metadata: dict) -> None:
            self.status_label.setText(f"Status: {metadata.get('status')} ({run_id})")
            self.statusBar().showMessage(f"Run finished with code {metadata.get('return_code')}.")
            self.current_run_dir = RUNS_DIR / run_id
            self.current_metadata = metadata
            self._load_runs()

        def _on_run_failed(self, run_id: str, message: str) -> None:
            self.status_label.setText(f"Status: failed ({run_id})")
            self.statusBar().showMessage("Run crashed before completion.")
            QMessageBox.critical(self, "Run Failed", message)
            self.current_run_dir = RUNS_DIR / run_id
            self._load_runs()

        def _cleanup_worker(self) -> None:
            self.run_button.setEnabled(True)
            self.delete_run_button.setEnabled(can_delete_run(self.current_metadata))
            if self.worker_thread is not None:
                self.worker_thread.deleteLater()
            self.worker_thread = None
            self.worker = None

        def _stop_worker(self) -> None:
            if self.worker is not None:
                self.worker.cancel()
            if self.worker_thread is not None:
                self.worker_thread.quit()
                if not self.worker_thread.wait(5000):
                    self.worker_thread.terminate()
                    self.worker_thread.wait(1000)
            self._cleanup_worker()

        def _confirm_stop_running_job(self) -> bool:
            if self.worker_thread is None:
                return True
            reply = QMessageBox.question(
                self,
                "Exit",
                "An analysis is still running. Stop it and exit?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            return reply == QMessageBox.Yes

        def _request_exit(self) -> None:
            if self._exiting:
                return
            if not self._confirm_stop_running_job():
                return
            self._exiting = True
            self._stop_worker()
            if self.image_viewer.movie is not None:
                self.image_viewer.movie.stop()
            app = QApplication.instance()
            if app is not None:
                app.quit()

        def closeEvent(self, event) -> None:  # pragma: no cover - UI path
            if self._exiting:
                event.accept()
                return
            if self.worker_thread is not None and not self._confirm_stop_running_job():
                event.ignore()
                return
            self._exiting = True
            self._stop_worker()
            if self.image_viewer.movie is not None:
                self.image_viewer.movie.stop()
            event.accept()

        def _on_run_selected(self, current: QListWidgetItem | None, previous: QListWidgetItem | None) -> None:
            del previous
            if current is None:
                return
            run_dir = Path(current.data(Qt.UserRole))
            self.current_run_dir = run_dir
            self.current_metadata = load_metadata(run_dir)
            self._refresh_run_views()

        def _refresh_run_views(self) -> None:
            metadata = self.current_metadata
            status = metadata.get("status", "unknown")
            self.status_label.setText(f"Status: {status} ({self.current_run_dir.name})")
            self.delete_run_button.setEnabled(can_delete_run(metadata) and self.worker_thread is None)

            images = metadata.get("images", [])
            self.current_image_index = 0
            if images:
                self._show_current_image()
            else:
                self.image_viewer.clear_view("No output images found for this run.")
                self.image_caption.setText("No images loaded")
            self._update_nav_buttons()

            self.logs_list.blockSignals(True)
            self.outputs_list.blockSignals(True)
            try:
                self.logs_list.clear()
                for relative_path in metadata.get("logs", []):
                    item = QListWidgetItem(Path(relative_path).name)
                    item.setData(Qt.UserRole, relative_path)
                    self.logs_list.addItem(item)

                self.outputs_list.clear()
                for relative_path in metadata.get("text_outputs", []):
                    item = QListWidgetItem(Path(relative_path).name)
                    item.setData(Qt.UserRole, relative_path)
                    self.outputs_list.addItem(item)

                if self.logs_list.count() > 0:
                    self.logs_list.setCurrentRow(0)
                    log_path = self.logs_list.currentItem().data(Qt.UserRole)
                    self.log_text.setPlainText(
                        read_text_file(self.current_run_dir, log_path)
                    )
                else:
                    self.log_text.setPlainText("")

                if self.outputs_list.count() > 0:
                    self.outputs_list.setCurrentRow(0)
                    output_path = self.outputs_list.currentItem().data(Qt.UserRole)
                    self.output_text.setPlainText(
                        read_text_file(self.current_run_dir, output_path)
                    )
                else:
                    self.output_text.setPlainText("")
            finally:
                self.logs_list.blockSignals(False)
                self.outputs_list.blockSignals(False)

        def _show_current_image(self) -> None:
            if self.current_run_dir is None:
                return
            images = self.current_metadata.get("images", [])
            if not images:
                return
            relative_path = images[self.current_image_index]
            image_path = self.current_run_dir / relative_path
            self.image_viewer.show_image(image_path)
            self.image_caption.setText(
                f"{self.current_image_index + 1} / {len(images)}  {Path(relative_path).name}"
            )

        def _step_image(self, delta: int) -> None:
            images = self.current_metadata.get("images", [])
            if not images:
                return
            self.current_image_index = (self.current_image_index + delta) % len(images)
            self._show_current_image()

        def _update_nav_buttons(self) -> None:
            image_count = len(self.current_metadata.get("images", []))
            enabled = image_count > 1
            self.prev_button.setEnabled(enabled)
            self.next_button.setEnabled(enabled)

        def _on_log_selected(self, current: QListWidgetItem | None, previous: QListWidgetItem | None) -> None:
            del previous
            if current is None or self.current_run_dir is None:
                return
            relative_path = current.data(Qt.UserRole)
            self.log_text.setPlainText(read_text_file(self.current_run_dir, relative_path))

        def _on_output_selected(self, current: QListWidgetItem | None, previous: QListWidgetItem | None) -> None:
            del previous
            if current is None or self.current_run_dir is None:
                return
            relative_path = current.data(Qt.UserRole)
            self.output_text.setPlainText(read_text_file(self.current_run_dir, relative_path))

        def _delete_selected_run(self) -> None:
            if self.current_run_dir is None:
                return
            if not can_delete_run(self.current_metadata):
                QMessageBox.information(
                    self,
                    "Run Not Deletable",
                    "Only completed or failed runs can be removed from the GUI.",
                )
                return
            reply = QMessageBox.question(
                self,
                "Delete Stored Run Output",
                f"Delete stored output for run {self.current_run_dir.name}?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return
            run_name = self.current_run_dir.name
            delete_run_directory(self.current_run_dir)
            self.current_run_dir = None
            self.current_metadata = {}
            self.statusBar().showMessage(f"Deleted run output for {run_name}.")
            self._load_runs()


def main() -> None:
    if QT_IMPORT_ERROR is not None:
        raise SystemExit(
            "PySide6 is required for the Qt GUI. Install it in the project environment and try again.\n"
            f"Original import error: {QT_IMPORT_ERROR}"
        )
    if not SCRIPT_PATH.exists():
        raise SystemExit(f"Missing script: {SCRIPT_PATH}")
    app = QApplication(sys.argv)
    window = RciMainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
