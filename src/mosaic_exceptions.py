"""Exceptions found by acq4_mosaic, written as a workbook to check by hand, and the decisions
made in it, read back to steer the next run.

Every run of ``acq4_mosaic`` writes ``exceptions_<date>-<time>.xlsx`` to its output directory,
with a sheet for each kind of exception and a README sheet:

``image kind``
    Images whose appearance (see :mod:`src.image_kind`) is ambiguous, or whose kind was set in
    the decision table. Decisions:
    BF, fluorescence, green, red, yellow, mixed, exclude. "fluorescence" keeps the colour of a
    fluorescence folder and uses yellow otherwise.
``image kind (automatic)``
    Images whose channel was set from a confident verdict of the classifier, overriding the
    folder name (e.g. a dye fill in a cell folder, which the folder rule calls brightfield).
    Listed for checking; the same decisions apply.
``refused matches``
    Groups whose registration would have moved them more than 250 um (usually the slice was
    moved between images). Decisions: stage (the default: leave at the stage position), accept
    (use the move), exclude (leave the group out), manual (leave at the stage position and list
    the group for a manual mosaic).
``annotation mismatches``
    Cells whose position, placed by hand in a Mosaic Editor file, is outside all of the cell's
    own higher-magnification images. Decisions: annotation (the default: keep it), inferred (use
    the centre of the cell's highest-magnification image), exclude (no position).

Fill in the ``decision`` column (leave it empty for the default) and give the file to the next
run with ``--decisions``. The workbook written by that run lists the same exceptions with the
decisions filled in, so one file can be kept up to date.
"""

import io
from pathlib import Path
from typing import Optional

import numpy as np

choices = {
    "image kind": ["BF", "fluorescence", "green", "red", "yellow", "mixed", "exclude"],
    "image kind (automatic)": ["BF", "fluorescence", "green", "red", "yellow", "mixed", "exclude"],
    "refused matches": ["stage", "accept", "exclude", "manual"],
    "annotation mismatches": ["annotation", "inferred", "exclude"],
}

# columns that identify a row (with the sheet), for matching decisions to images, groups and cells
keys = {
    "image kind": ("slice", "image"),
    "image kind (automatic)": ("slice", "image"),
    "refused matches": ("slice", "group", "channel", "mag", "run"),
    "annotation mismatches": ("slice", "cell"),
}

_readme = [
    ("Sheet", "What is listed", "decision (empty = the default, first listed)"),
    ("image kind", "Images whose appearance is ambiguous (or set in this table): need a decision", ", ".join(choices["image kind"])),
    ("image kind (automatic)", "Images whose channel was set from the classifier's confident verdict, overriding the folder name: check, decide only to change", ", ".join(choices["image kind (automatic)"])),
    ("refused matches", "Group moves of more than 250 um that were not used", ", ".join(choices["refused matches"])),
    ("annotation mismatches", "Hand-placed cell positions outside all of the cell's own images", ", ".join(choices["annotation mismatches"])),
    ("", "", ""),
    ("How to use", "Fill in the decision column, save, and run acq4_mosaic again with --decisions <this file>.", ""),
    ("", "image kind: fluorescence keeps the colour of a fluorescence folder, otherwise yellow. mixed = fluorescence with brightfield for context.", ""),
    ("", "refused matches: stage = leave at the stage position, accept = use the move, exclude = leave the group out, manual = stage, listed for a manual mosaic.", ""),
    ("", "annotation mismatches: annotation = keep the hand-placed position, inferred = centre of the cell's highest-magnification image, exclude = no position.", ""),
]


def _key(sheet: str, row: dict) -> tuple:
    """The identifying values of a row, as strings.

    Parameters
    ----------
    sheet : str
        Sheet name.
    row : dict
        The row.

    Returns
    -------
    tuple of str
        The values of the sheet's key columns.
    """
    out = []
    for k in keys[sheet]:
        v = row.get(k, "")
        if isinstance(v, float) and v.is_integer():
            v = int(v)
        out.append(str(v))
    return tuple(out)


class Decisions:
    """Decisions read from an exceptions workbook (empty if there is none).

    Parameters
    ----------
    path : str or pathlib.Path, optional
        The workbook.

    Raises
    ------
    ValueError
        If a decision is not one of the choices for its sheet.
    """

    def __init__(self, path: Optional[Path] = None):
        """Read the decisions from a workbook (none if path is None).

        Parameters
        ----------
        path : str or pathlib.Path, optional
            The workbook.
        """
        self.path = Path(path) if path is not None else None
        self.table = {sheet: {} for sheet in choices}
        if self.path is None:
            return
        import openpyxl

        wb = openpyxl.load_workbook(self.path, read_only=True, data_only=True)
        for sheet in choices:
            if sheet not in wb.sheetnames:
                continue
            rows = list(wb[sheet].iter_rows(values_only=True))
            if not rows:
                continue
            header = [str(h) if h is not None else "" for h in rows[0]]
            for values in rows[1:]:
                row = dict(zip(header, values))
                decision = row.get("decision")
                if decision is None or str(decision).strip() == "":
                    continue
                decision = str(decision).strip()
                allowed = {c.lower(): c for c in choices[sheet]}
                if decision.lower() not in allowed:
                    raise ValueError(f"{self.path.name}, sheet {sheet}: decision {decision!r} is not one of {choices[sheet]}")
                self.table[sheet][_key(sheet, row)] = allowed[decision.lower()]

    def get(self, sheet: str, **row) -> Optional[str]:
        """The decision for an image, group or cell, if one was made.

        Parameters
        ----------
        sheet : str
            Sheet name (the two image-kind sheets are looked up together: a decision made in
            either applies to the image).
        **row
            The key columns of the sheet (e.g. ``slice=..., image=...``).

        Returns
        -------
        str or None
            The decision, or None for the default.
        """
        if sheet.startswith("image kind"):  # the two image-kind sheets share their decisions
            k = _key("image kind", row)
            return self.table["image kind"].get(k) or self.table["image kind (automatic)"].get(k)
        return self.table[sheet].get(_key(sheet, row))

    def __len__(self) -> int:
        """Number of decisions read.

        Returns
        -------
        int
            The number of rows with a decision, over all sheets.
        """
        return sum(len(v) for v in self.table.values())


def write_workbook(filename: Path, exceptions: dict, decisions: Optional[Decisions] = None):
    """Write the exceptions of one or more slices as a workbook, one sheet per kind.

    Parameters
    ----------
    filename : pathlib.Path
        The .xlsx file.
    exceptions : dict
        Sheet name -> list of rows (dicts). A row may hold "_thumbnail", a 2D image shown in the
        sheet.
    decisions : Decisions, optional
        Decisions to fill in (from the workbook given to this run).
    """
    import openpyxl
    from openpyxl.drawing.image import Image as XLImage
    from openpyxl.styles import Font, PatternFill
    from openpyxl.worksheet.datavalidation import DataValidation

    wb = openpyxl.Workbook()
    readme = wb.active
    readme.title = "README"
    for r in _readme:
        readme.append(list(r))
    readme["A1"].font = readme["B1"].font = readme["C1"].font = Font(bold=True)
    readme.column_dimensions["A"].width = 24
    readme.column_dimensions["B"].width = 110
    readme.column_dimensions["C"].width = 60
    for sheet in choices:
        rows = exceptions.get(sheet, [])
        ws = wb.create_sheet(sheet)
        has_thumbs = any("_thumbnail" in r for r in rows)
        columns = ["decision"] + (["thumbnail"] if has_thumbs else [])
        for r in rows:
            for k in r:
                if k not in columns and not k.startswith("_"):
                    columns.append(k)
        if not rows:
            columns.append("(none in this run)")
        ws.append(columns)
        for c in ws[1]:
            c.font = Font(bold=True)
        fill = PatternFill("solid", fgColor="FFF2CC")
        dv = DataValidation(type="list", formula1='"' + ",".join(choices[sheet]) + '"', allow_blank=True)
        ws.add_data_validation(dv)
        for i, r in enumerate(rows, start=2):
            decided = decisions.get(sheet, **r) if decisions is not None else None
            values = [decided or r.get("decision")] + ([""] if has_thumbs else [])
            values += [_cell(r.get(k)) for k in columns[len(values):]]
            ws.append(values)
            ws.cell(row=i, column=1).fill = fill
            dv.add(ws.cell(row=i, column=1))
            if "_thumbnail" in r:
                img = XLImage(io.BytesIO(r["_thumbnail"]))
                img.width, img.height = 120, int(120 * r.get("_aspect", 0.75))
                ws.add_image(img, ws.cell(row=i, column=2).coordinate)
                ws.row_dimensions[i].height = img.height * 0.78
        ws.column_dimensions["A"].width = 14
        if has_thumbs:
            ws.column_dimensions["B"].width = 18
        ws.freeze_panes = "B2"
        for j, name in enumerate(columns, start=1):
            if name in ("slice", "image", "notes", "slice_dir"):
                ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 40
    filename.parent.mkdir(parents=True, exist_ok=True)
    wb.save(filename)


def _cell(v):
    """A value as Excel can hold it (numbers rounded, lists joined).

    Parameters
    ----------
    v : object
        The value.

    Returns
    -------
    object
        A number, string or None.
    """
    if isinstance(v, (np.floating, float)):
        return round(float(v), 3)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v)
    return v
