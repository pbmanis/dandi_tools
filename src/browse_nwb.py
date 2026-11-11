"""Browse the NWB files in a directory, grouped by recording day, and display each
protocol with display_nwb.

The files are found in the directory and its subdirectories (so a dandiset directory,
with one subdirectory per subject, also works). The day is taken from the file name,
e.g. DCN_IC_inj~2021.04.05~S00C00.nwb or sub-X_ses-DCN-IC-inj~2025-07-24~S00C00_icephys.nwb.

Usage:
    uv run browse_nwb [directory]

Keys:
    Ctrl+Right / Ctrl+Left: next / previous protocol (continues into the next / previous file)
    Ctrl+Down / Ctrl+Up: next / previous file
    PgDown / PgUp: next / previous day
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Union

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt6 import QtCore, QtGui, QtWidgets

from src.display_nwb import ReadNWB

re_day = re.compile(r"~(?P<year>\d{4})[.-](?P<month>\d{2})[.-](?P<day>\d{2})~")
re_cell = re.compile(r"~(?P<cell>S\d+C\d+)")

FILE_ROLE = QtCore.Qt.ItemDataRole.UserRole


def day_of(path: Path) -> str:
    """The recording day (yyyy.mm.dd) from an NWB file name"""
    daymatch = re_day.search(path.name)
    if daymatch is None:
        return "unknown day"
    return f"{daymatch['year']}.{daymatch['month']}.{daymatch['day']}"


def find_nwb_files(directory: Union[str, Path]) -> dict:
    """Find the NWB files in a directory and its subdirectories

    Returns:
        dict: {day: [files]}, sorted by day and file name
        int: the number of entries that are not readable files (e.g., broken links)
    """
    days = {}
    missing = 0
    for f in sorted(Path(directory).rglob("*.nwb"), key=lambda f: f.name):
        if not f.is_file():  # e.g., a link to a file that no longer exists
            missing += 1
            continue
        days.setdefault(day_of(f), []).append(f)
    return dict(sorted(days.items())), missing


class BrowseNWB(QtWidgets.QMainWindow):
    def __init__(self, directory: Union[str, Path, None] = None):
        super().__init__()
        self.setWindowTitle("NWB browser")
        self.directory = None
        self.days = {}  # {day: [files]}
        self.files = []  # all files, in display order
        self.reader = None  # ReadNWB for the current file
        self.current_file = None
        self.protocol_names = []
        self.protocol_index = 0

        # left side: directory chooser, the days and files, and file information
        choose_button = QtWidgets.QPushButton("Choose directory…")
        choose_button.clicked.connect(self.choose_directory)
        self.dir_label = QtWidgets.QLabel("No directory chosen")
        self.dir_label.setWordWrap(True)
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(["Day / cell"])
        self.tree.currentItemChanged.connect(self.on_tree_selection)
        self.info = QtWidgets.QPlainTextEdit()
        self.info.setReadOnly(True)
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.addWidget(choose_button)
        left_layout.addWidget(self.dir_label)
        left_layout.addWidget(self.tree, stretch=3)
        left_layout.addWidget(self.info, stretch=1)

        # right side: navigation, and the plot
        self.protocol_combo = QtWidgets.QComboBox()
        self.protocol_combo.currentIndexChanged.connect(self.show_protocol)
        nav = QtWidgets.QHBoxLayout()
        nav.addWidget(QtWidgets.QLabel("Protocol:"))
        nav.addWidget(self.protocol_combo, stretch=1)
        for label, tip, shortcut, slot in [
            ("◀ Day", "Previous day (PgUp)", "PgUp", lambda: self.step_day(-1)),
            ("◀ File", "Previous file (Ctrl+Up)", "Ctrl+Up", lambda: self.step_file(-1)),
            ("◀ Protocol", "Previous protocol (Ctrl+Left)", "Ctrl+Left", lambda: self.step_protocol(-1)),
            ("Protocol ▶", "Next protocol (Ctrl+Right)", "Ctrl+Right", lambda: self.step_protocol(1)),
            ("File ▶", "Next file (Ctrl+Down)", "Ctrl+Down", lambda: self.step_file(1)),
            ("Day ▶", "Next day (PgDown)", "PgDown", lambda: self.step_day(1)),
        ]:
            button = QtWidgets.QPushButton(label)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            nav.addWidget(button)
            QtGui.QShortcut(QtGui.QKeySequence(shortcut), self, activated=slot)
        self.figure = Figure(figsize=(8, 10), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        right = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.addLayout(nav)
        right_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        right_layout.addWidget(self.canvas, stretch=1)

        splitter = QtWidgets.QSplitter()
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([320, 980])
        self.setCentralWidget(splitter)

        if directory is not None:
            self.load_directory(directory)

    def choose_directory(self):
        start = str(self.directory) if self.directory is not None else str(Path.cwd())
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Choose a directory of NWB files", start
        )
        if directory:
            self.load_directory(directory)

    def load_directory(self, directory: Union[str, Path]):
        """List the NWB files in the directory by day, and show the first one"""
        self.directory = Path(directory)
        self.days, missing = find_nwb_files(self.directory)
        self.files = [f for files in self.days.values() for f in files]
        self.dir_label.setText(str(self.directory))
        self.current_file = None
        self.tree.blockSignals(True)
        self.tree.clear()
        for day, files in self.days.items():
            day_item = QtWidgets.QTreeWidgetItem([f"{day:s}  ({len(files):d})"])
            for f in files:
                cellmatch = re_cell.search(f.name)
                label = cellmatch["cell"] if cellmatch is not None else f.stem
                file_item = QtWidgets.QTreeWidgetItem([label])
                file_item.setData(0, FILE_ROLE, str(f))
                file_item.setToolTip(0, str(f))
                day_item.addChild(file_item)
            self.tree.addTopLevelItem(day_item)
        self.tree.expandAll()
        self.tree.blockSignals(False)
        message = f"{len(self.files):d} NWB files on {len(self.days):d} days in {self.directory!s}"
        if missing > 0:
            message += f"; {missing:d} .nwb entries skipped (not found: broken links?)"
        self.statusBar().showMessage(message)
        if len(self.files) > 0:
            self.select_file(self.files[0])
        else:
            self.info.setPlainText(
                "No NWB files found."
                + (f"\n{missing:d} .nwb entries are links to files that do not exist." if missing > 0 else "")
            )
            self.figure.clear()
            self.canvas.draw_idle()

    def file_item(self, f: Path):
        """The tree item for a file"""
        for i in range(self.tree.topLevelItemCount()):
            day_item = self.tree.topLevelItem(i)
            for j in range(day_item.childCount()):
                if day_item.child(j).data(0, FILE_ROLE) == str(f):
                    return day_item.child(j)
        return None

    def select_file(self, f: Path):
        item = self.file_item(f)
        if item is not None:
            self.tree.setCurrentItem(item)

    def on_tree_selection(self, current, previous):
        if current is None:
            return
        filename = current.data(0, FILE_ROLE)
        if filename is None:  # a day: show its first file
            if current.childCount() > 0:
                self.tree.setCurrentItem(current.child(0))
            return
        self.show_file(Path(filename))

    def show_file(self, f: Path, protocol_index: int = 0):
        """Read a file and show one of its protocols"""
        if f != self.current_file:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
            try:
                reader = ReadNWB()
                reader.readfile(f)
            except Exception as e:
                self.info.setPlainText(f"Could not read {f!s}:\n{e!r}")
                self.statusBar().showMessage(f"Could not read {f.name:s}")
                self.figure.clear()
                self.canvas.draw_idle()
                return
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()
            if self.reader is not None and self.reader.pick_cid is not None:
                self.canvas.mpl_disconnect(self.reader.pick_cid)
            self.reader = reader
            self.current_file = f
            self.protocol_names = list(reader.protocols.keys())
            self.protocol_combo.blockSignals(True)
            self.protocol_combo.clear()
            self.protocol_combo.addItems(
                [f"{name:s}  ({p.mode:s}, {p.nsweeps:d} sweeps)" for name, p in reader.protocols.items()]
            )
            self.protocol_combo.blockSignals(False)
            self.info.setPlainText(self.file_info_text(f))
        if protocol_index < 0:
            protocol_index = len(self.protocol_names) + protocol_index
        self.show_protocol(protocol_index)

    def file_info_text(self, f: Path) -> str:
        lines = [f.name, ""]
        lines += [f"{key:s}: {value!s}" for key, value in self.reader.info.items()]
        lines += ["", f"{len(self.protocol_names):d} protocols:"]
        lines += [
            f"  {name:s} ({p.mode:s}, {p.nsweeps:d} sweeps{', optical' if p.ostim is not None else '':s})"
            for name, p in self.reader.protocols.items()
        ]
        return "\n".join(lines)

    def show_protocol(self, index: int):
        """Plot one protocol of the current file"""
        if self.reader is None or len(self.protocol_names) == 0:
            self.figure.clear()
            self.canvas.draw_idle()
            return
        index = min(max(index, 0), len(self.protocol_names) - 1)
        self.protocol_index = index
        self.protocol_combo.blockSignals(True)
        self.protocol_combo.setCurrentIndex(index)
        self.protocol_combo.blockSignals(False)
        protocol = self.reader.protocols[self.protocol_names[index]]
        try:
            self.reader.plot_traces(
                protocol,
                fig=self.figure,
                nprotocols=len(self.protocol_names),
                protocol_number=index + 1,
            )
        except Exception as e:
            self.figure.clear()
            self.statusBar().showMessage(f"Could not plot {protocol.protocol:s}: {e!r}")
        else:
            self.statusBar().showMessage(
                f"{self.current_file.name:s}: {protocol.protocol:s} "
                f"({index + 1:d} of {len(self.protocol_names):d})"
            )
        self.canvas.draw_idle()

    def step_protocol(self, delta: int):
        """Go to the next/previous protocol, continuing into the next/previous file"""
        if self.current_file is None:
            return
        index = self.protocol_index + delta
        if 0 <= index < len(self.protocol_names):
            self.show_protocol(index)
        else:
            self.step_file(delta, protocol_index=0 if delta > 0 else -1)

    def step_file(self, delta: int, protocol_index: int = 0):
        if self.current_file is None:
            return
        index = self.files.index(self.current_file) + delta
        if 0 <= index < len(self.files):
            self.tree.blockSignals(True)
            self.select_file(self.files[index])
            self.tree.blockSignals(False)
            self.show_file(self.files[index], protocol_index=protocol_index)

    def step_day(self, delta: int):
        if self.current_file is None:
            return
        days = list(self.days.keys())
        index = days.index(day_of(self.current_file)) + delta
        if 0 <= index < len(days):
            self.select_file(self.days[days[index]][0])


def main():
    parser = argparse.ArgumentParser(description="Browse the NWB files in a directory")
    parser.add_argument(
        dest="directory", type=str, nargs="?", default=None, help="directory of NWB files"
    )
    args = parser.parse_args()
    app = QtWidgets.QApplication(sys.argv[:1])
    window = BrowseNWB(args.directory)
    window.resize(1300, 950)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
