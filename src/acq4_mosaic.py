"""Experimental: assemble slice mosaics from the images and videos in acq4 data directories.

For one slice directory, this:

1. Finds the images (image_NNN.tif) and videos (video_NNN.ma) in the slice directory and its
   subdirectories (fluorescence folders, cell folders, ...), reading the position, pixel size and
   objective of each from the .index files. Protocol directories are not searched.
2. Assigns a channel to each image from the name of its folder (FITC, GFP, YFP -> green;
   TRITC, mCherry, tdTomato, DiI -> red; Lucifer yellow, Alexa -> yellow). Images in other
   folders (the slice and cell folders) are treated as brightfield: some are fluorescence, but
   the intensities do not tell them apart reliably (brightfield images can be mostly dark bath,
   and tissue autofluorescence can make fluorescence images look like brightfield).
   Folders named garbage, bad, etc. are left out. Videos that are too large (usually a recording
   that was not stopped) are left out; the others (z-stacks taken while focusing by hand) are
   projected to one image by a multi-scale focus stack, which keeps what is in focus in each
   frame and suppresses the out-of-focus flare of bright somata (see src/stack_projection.py).
3. Stitches the lowest-magnification brightfield images of the slice (the "base"). Neighbouring
   tiles usually overlap by only a 40-200 um strip, and the stage error can be larger than that,
   so every pair that shares an edge is registered (normalized cross-correlation of band-passed
   images, within 200 um of the stage positions). Strings running along a strip fix the shift
   across the strip but hardly along it, so each pair is weighted, direction by direction, by
   how sharply its correlation peak is defined (the curvature of the peak). The translations
   are then solved for together by least squares, with a weak pull toward the stage positions,
   dropping the most inconsistent pair until the rest agree.
4. Registers the other images as groups (one folder, channel and magnification at a time),
   lowest magnification first, against everything placed so far: first the whole group within a
   large window (the manual corrections in the .index files are mostly whole-folder offsets, up to
   ~1.5 mm), then each image within a small window. Images of different kinds (fluorescence vs
   brightfield) are compared by gradient magnitude, which does not depend on the sign of the
   contrast: nylon strings and tissue edges are dark in brightfield but bright (autofluorescent)
   in fluorescence. If a group cannot be registered at its recorded magnification, the other
   objectives on the rig are tried (the objective setting is sometimes not changed when the
   objective is). Placements that are not clearly better than the alternatives, or that would
   move the images more than 250 um (usually the slice was moved between images), are not used:
   the image stays at its stage position, and this is reported.
5. Renders the slice (brightfield in grey, fluorescence as transparent coloured overlays),
   before and after registration, and a close-up of each cell with its higher-magnification
   images, and writes a JSON report with every image's placement and the problems found.
   The report carries each image's .index metadata forward, with its placement (corrected
   magnification and position) in the form of an acq4 transform. The placements are also
   written as an acq4 Mosaic Editor file (.mosaic), and the video projections as TIFFs.

Only translations are used (no rotation, scaling or warping), as the stage only translates.

Usage:
    uv run acq4_mosaic <slice_dir> [<slice_dir> ...] -o <output_dir>
    uv run acq4_mosaic --random 4 --root <acq4 data root> -o <output_dir>
"""

import argparse
import itertools
import json
import os
import random
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as mpl
import numpy as np
import scipy.fft
import scipy.ndimage as ndi
import tifffile
from matplotlib.patches import Rectangle
from pyqtgraph import configfile

import MetaArray as EM
from src.image_kind import classify, measure, thumbnail_jpeg
from src.mosaic_exceptions import Decisions, write_workbook
from src.stack_projection import project_stack

re_mag = re.compile(r"^\s*(?P<mag>\d+(?:\.\d+)?)\s*x", re.IGNORECASE)
re_image = re.compile(r"^(image|video)_\d+.*\.(tif|tiff|ma)$", re.IGNORECASE)
re_cell = re.compile(r"^cell_\d+$")
re_slice = re.compile(r"^slice_\d+$")
re_sweep = re.compile(r"^\d{3}$")

# a group registration that would move images further than this (um) is not used: such moves are
# usually the slice having been moved between images, not stage error
max_move_um = 250.0

# folder name -> channel. The innermost folder with a match wins.
channel_patterns = [
    ("yellow", re.compile(r"lucifer|(^|[^a-z])ly([^a-z]|$)|alexa|biocytin", re.IGNORECASE)),
    ("red", re.compile(r"tritc|trittc|mcherry|tdtom|dii|rfp|cy3", re.IGNORECASE)),
    ("green", re.compile(r"fitc|ftic|gfp|yfp|chr2", re.IGNORECASE)),
    ("BF", re.compile(r"(^|[^a-z])bf([^a-z]|$)|bright", re.IGNORECASE)),
]
exclude_pattern = re.compile(r"garbage|junk|trash|(^|[^a-z])bad([^a-z]|$)", re.IGNORECASE)

# overlay colours (RGB) for fluorescence channels
channel_colors = {
    "green": (0.1, 1.0, 0.1),
    "red": (1.0, 0.15, 0.1),
    "yellow": (1.0, 0.9, 0.0),
}


_tile_serial = itertools.count()  # unique, never reused (unlike id(), which is reused after a copy is freed)


@dataclass
class Tile:
    """One image (or the projection of one video), and where it goes

    For a video, the image is a projection of its frames (see :func:`load_tile`). Positions are
    world (stage) coordinates in um, with +x to the right and +y up; the image data are stored in
    that orientation (row 0 at the top).

    Attributes
    ----------
    path : pathlib.Path
        The image (.tif) or video (.ma) file.
    group : str
        Folder, relative to the slice directory ("." for the slice directory).
    kind : str
        "image" or "video".
    objective : str
        The objective named in the .index file.
    mag : float
        Magnification in use (may differ from the .index after a correction).
    pixel : float
        Pixel size in use, um.
    xmin, ymax : float
        World position of the left and top edges at the stage position, um.
    shape : tuple
        (rows, cols).
    flip_x, flip_y : bool
        Whether the stored image had to be flipped to put +x to the right and +y at the top.
    z : float
        Focus (z) position, um.
    timestamp : float
        Acquisition time (s since the epoch), from the .index file.
    user_offset : tuple
        The userTransform position in the .index file (manual correction), um.
    cell : str or None
        The cell folder the image is in, if any.
    nbytes : int
        File size.
    channel : str
        "BF", "green", "red" or "yellow".
    channel_from : str
        "folder" (from the folder name) or "default".
    excluded : str
        The reason the tile is not used, or "" if it is used.
    dx, dy : float
        Correction added to the stage position, um.
    placement : str
        How the tile was placed: "stage", "stitched", "group", "refined" or "stage+stitched"
        (aligned with the rest of its group, which as a whole is at the stage position).
    quality : float
        NCC peak of the registration that placed the tile.
    notes : list of str
        Notes for the report.
    data : numpy.ndarray or None
        The image (float32), once loaded.
    frames : int
        Number of frames (videos).
    projection : dict or None
        How a video's frames were projected to one image.
    appearance : dict or None
        Brightfield vs fluorescence measurements and the rule's verdict (src.image_kind).
    mixed : bool
        Fluorescence with brightfield for context (set in the decision table).
    serial : int
        Unique number, the key of the resampling cache (``id()`` is reused after an object is
        freed, which made the temporary copies tried at other magnifications pick up each
        other's cached images).
    index_info : dict
        The image's whole .index entry, as recorded.
    index_mag, index_pixel : float
        Magnification and pixel size (um) according to the .index file.
    stage_xmin, stage_ymax : float
        Left and top edges (um) according to the .index file.
    """

    path: Path
    group: str  # folder, relative to the slice directory ("." for the slice directory)
    kind: str  # "image" or "video"
    objective: str
    mag: float
    pixel: float  # um per pixel
    xmin: float  # world position of the left edge, um
    ymax: float  # world position of the top edge, um
    shape: tuple  # (rows, cols)
    flip_x: bool  # the image must be flipped to put +x to the right
    flip_y: bool  # ... and +y at the top
    z: float  # um
    timestamp: float
    user_offset: tuple  # the userTransform in the .index file, um
    cell: Optional[str]  # the cell folder it is in, if any
    nbytes: int
    channel: str = "BF"  # BF, green, red or yellow
    channel_from: str = "default"  # folder or default
    excluded: str = ""  # the reason, if the tile is not used
    dx: float = 0.0  # correction to the stage position, um
    dy: float = 0.0
    placement: str = "stage"  # stage, stitched, group, refined, stage+stitched
    quality: float = float("nan")
    notes: list = field(default_factory=list)
    data: Optional[np.ndarray] = field(default=None, repr=False)
    frames: int = 0  # videos: number of frames
    appearance: Optional[dict] = None  # brightfield vs fluorescence measurements (src.image_kind)
    mixed: bool = False  # fluorescence with brightfield for context (from the decision table)
    serial: int = field(default_factory=lambda: next(_tile_serial), repr=False)  # key for caches
    projection: Optional[dict] = None  # videos: how the frames were projected to one image
    # as recorded in the .index file, before any correction
    index_info: dict = field(default_factory=dict, repr=False)  # the image's whole .index entry
    index_mag: float = float("nan")
    index_pixel: float = float("nan")
    stage_xmin: float = float("nan")
    stage_ymax: float = float("nan")

    @property
    def name(self) -> str:
        """The file name, with its folder relative to the slice directory.

        Returns
        -------
        str
            e.g. ``"image_003.tif"`` or ``"cell_001/video_000.ma"``.
        """
        return self.path.name if self.group == "." else f"{self.group}/{self.path.name}"

    @property
    def width(self) -> float:
        """Width of the image in um.

        Returns
        -------
        float
            Columns times the pixel size.
        """
        return self.shape[1] * self.pixel

    @property
    def height(self) -> float:
        """Height of the image in um.

        Returns
        -------
        float
            Rows times the pixel size.
        """
        return self.shape[0] * self.pixel

    def bounds(self, placed: bool = True) -> tuple:
        """(xmin, xmax, ymin, ymax) in um

        Parameters
        ----------
        placed : bool, optional
            Include the correction (dx, dy) (default True); False gives the stage position.

        Returns
        -------
        tuple of float
            (xmin, xmax, ymin, ymax), um.
        """
        dx, dy = (self.dx, self.dy) if placed else (0.0, 0.0)
        x0 = self.xmin + dx
        y1 = self.ymax + dy
        return (x0, x0 + self.width, y1 - self.height, y1)

    def center(self, placed: bool = True) -> tuple:
        """Centre of the image in world coordinates.

        Parameters
        ----------
        placed : bool, optional
            Include the correction (dx, dy) (default True).

        Returns
        -------
        tuple of float
            (x, y), um.
        """
        b = self.bounds(placed)
        return ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2)

    def with_mag(self, mag: float) -> "Tile":
        """A copy of the tile as if it had been taken with another objective (same center)

        The pixel size is scaled by the ratio of the magnifications. The image data are shared, not
        copied; notes are not copied.

        Parameters
        ----------
        mag : float
            The other magnification.

        Returns
        -------
        Tile
            The copy.
        """
        cx, cy = self.center(placed=False)
        pixel = self.pixel * self.mag / mag
        t = Tile(**{k: getattr(self, k) for k in self.__dataclass_fields__ if k not in ("notes", "data", "serial")})
        t.mag, t.pixel = mag, pixel
        t.xmin = cx - t.width / 2
        t.ymax = cy + t.height / 2
        t.data = self.data
        return t

    def mosaic_user_transform(self) -> dict:
        """The acq4 userTransform that moves the image from its .index position to its placement here.

        acq4 places an image by its .index ``transform`` (stage position and pixel size), then
        applies the item's ``userTransform`` (in world coordinates, m). This is the userTransform
        that puts the image where it was placed here, including a corrected magnification (as a
        scale about the world origin). It replaces, not adds to, the userTransform in the .index
        file, as in the acq4 Mosaic Editor.

        Returns
        -------
        dict
            ``{"pos": [x, y], "scale": [k, k], "angle": 0.0}``, pos in m.
        """
        k = self.pixel / self.index_pixel
        x0, _, _, y1 = self.bounds()
        pos = ((x0 - k * self.stage_xmin) * 1e-6, (y1 - k * self.stage_ymax) * 1e-6)
        return {"pos": [float(pos[0]), float(pos[1])], "scale": [float(k), float(k)], "angle": 0.0}

    def corrected_transform(self) -> Optional[dict]:
        """The image's .index ``transform`` with the corrections made here applied.

        The same form as the .index transform (``pos`` in m, ``scale`` in m per pixel with the
        signs of the original, ``angle`` in degrees), so it places the image directly, with no
        userTransform.

        Returns
        -------
        dict or None
            The corrected transform, or None if the .index entry has no transform.
        """
        tr = self.index_info.get("transform")
        if tr is None:
            return None
        u = self.mosaic_user_transform()
        k = u["scale"][0]
        pos = [float(v) for v in tr["pos"]]
        scale = [float(v) for v in tr["scale"]]
        pos[0], pos[1] = k * pos[0] + u["pos"][0], k * pos[1] + u["pos"][1]
        scale[0], scale[1] = k * scale[0], k * scale[1]
        return {"pos": pos, "scale": scale, "angle": float(tr.get("angle", 0.0))}


def json_safe(value):
    """Convert .index values (OrderedDicts, tuples, numpy scalars and arrays) to JSON types.

    Parameters
    ----------
    value : object
        A value read from a .index file.

    Returns
    -------
    object
        The same value as dicts, lists, floats, ints, strings, bools and None.
    """
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def file_number(path: Path) -> int:
    """The number in an acq4 file name (``image_012.tif`` -> 12).

    Parameters
    ----------
    path : pathlib.Path
        The file.

    Returns
    -------
    int
        The number, or -1 if the name has none.
    """
    m = re.search(r"(\d+)", path.stem)
    return int(m.group(1)) if m else -1


def compass(dx: float, dy: float) -> str:
    """The direction of a displacement in words (+x right, +y up).

    Parameters
    ----------
    dx, dy : float
        The displacement.

    Returns
    -------
    str
        e.g. "upper left" or "right".
    """
    angle = np.degrees(np.arctan2(dy, dx))
    names = ["right", "upper right", "top", "upper left", "left", "lower left", "bottom", "lower right"]
    return names[int(np.round(angle / 45.0)) % 8]


_provenance = None


def provenance() -> str:
    """Repository name, version and git commit, for stamping figures.

    The commit is marked "+ uncommitted changes" when files under ``src/`` differ from it (or
    are not yet committed), since then the commit alone does not identify the code.

    Returns
    -------
    str
        e.g. "dandi_acq4 0.2.0, git 64bf21f + uncommitted changes".
    """
    global _provenance
    if _provenance is None:
        import importlib.metadata
        import subprocess

        try:
            version = importlib.metadata.version("dandi_acq4")
        except importlib.metadata.PackageNotFoundError:
            version = "unknown version"
        root = Path(__file__).resolve().parent.parent
        try:
            commit = subprocess.run(
                ["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, check=True
            ).stdout.strip()
            changed = subprocess.run(
                ["git", "status", "--porcelain", "--", "src"], cwd=root, capture_output=True, text=True, check=True
            ).stdout.strip()
            commit = f"git {commit}" + (" + uncommitted changes" if changed else "")
        except (OSError, subprocess.CalledProcessError):
            commit = "git commit unknown"
        _provenance = f"dandi_acq4 {version}, {commit}"
    return _provenance


def read_index(directory: Path) -> dict:
    """Read the .index file of an acq4 directory.

    Parameters
    ----------
    directory : pathlib.Path
        The directory.

    Returns
    -------
    dict
        The .index contents, or ``{"__error__": message}`` if the file could not be read (a few
        .index files are damaged).
    """
    try:
        return configfile.readConfigFile(str(directory / ".index"))
    except Exception as e:  # a few .index files are damaged
        return {"__error__": repr(e)}


def magnification(objective: str) -> float:
    """The magnification in an objective name.

    Parameters
    ----------
    objective : str
        The objective as written in the .index file, e.g. "40x 0.8na ACHROPLAN".

    Returns
    -------
    float
        The magnification (40.0 in the example), or NaN if the name does not start with one.
    """
    m = re_mag.match(objective or "")
    return float(m["mag"]) if m is not None else float("nan")


def folder_channel(group: str) -> Optional[str]:
    """The channel given by the innermost folder name that has one

    Parameters
    ----------
    group : str
        Folder path relative to the slice directory.

    Returns
    -------
    str or None
        "yellow", "red", "green" or "BF" (see ``channel_patterns``), or None if no folder name
        names a channel.
    """
    for part in reversed(Path(group).parts):
        for channel, pattern in channel_patterns:
            if pattern.search(part):
                return channel
    return None


def image_directories(slice_dir: Path):
    """The slice directory and its subdirectories, except protocol (sweep) directories

    A directory is taken to hold a protocol if it has subdirectories named like sweeps (``000``,
    ``001``, ...); it is skipped, with everything under it.

    Parameters
    ----------
    slice_dir : pathlib.Path
        The slice directory.

    Yields
    ------
    pathlib.Path
        Each directory to search for images.
    """
    for dirpath, dirnames, filenames in os.walk(slice_dir):
        d = Path(dirpath)
        has_sweeps = any(re_sweep.match(n) for n in dirnames)
        dirnames[:] = sorted(n for n in dirnames if not n.startswith(".") and not has_sweeps)
        if not has_sweeps:
            yield d


def scan_slice(slice_dir: Path, max_video_mb: float = 500.0) -> tuple:
    """Find the images and videos of a slice. Returns (tiles, notes).

    Positions, pixel sizes, objectives and the manual userTransform come from the .index files.
    Images in a directory that are not in its .index have no position and are not used.

    Parameters
    ----------
    slice_dir : pathlib.Path
        The slice directory.
    max_video_mb : float, optional
        Videos larger than this (MB) are excluded (default 500): usually a recording that was not
        stopped.

    Returns
    -------
    tiles : list of Tile
        One per image or video listed in a .index file (excluded tiles included, with the reason).
    notes : list of str
        Problems found (unreadable .index files, missing files, images without positions).
    """
    tiles = []
    notes = []
    for d in image_directories(slice_dir):
        group = d.relative_to(slice_dir).as_posix()
        index = read_index(d) if (d / ".index").exists() else {}
        if "__error__" in index:
            notes.append(f"{group}: .index could not be read ({index['__error__'][:120]})")
            continue
        listed = set()
        for name, info in index.items():
            if not isinstance(info, dict) or not re_image.match(name):
                continue
            listed.add(name)
            path = d / name
            if not path.exists():
                notes.append(f"{group}/{name}: in .index, but the file is missing")
                continue
            transform = info.get("transform")
            if transform is None:
                notes.append(f"{group}/{name}: no transform in .index")
                continue
            sx, sy = (s * 1e6 for s in transform["scale"][:2])
            x0, y0 = (p * 1e6 for p in transform["pos"][:2])
            z = transform["pos"][2] * 1e6 if len(transform["pos"]) > 2 else float("nan")
            region = info.get("region", [0, 0, 0, 0])
            binning = info.get("binning", [1, 1])
            shape = (int(region[3] // binning[1]), int(region[2] // binning[0]))
            user = info.get("userTransform", {}).get("pos", (0.0, 0.0))
            parts = Path(group).parts
            cell = next((p for p in parts if re_cell.match(p)), None)
            t = Tile(
                path=path,
                group=group,
                kind="video" if name.lower().startswith("video") else "image",
                objective=str(info.get("objective", "")),
                mag=magnification(info.get("objective", "")),
                pixel=abs(sx),
                xmin=x0 if sx > 0 else x0 + shape[1] * sx,
                ymax=y0 if sy < 0 else y0 + shape[0] * sy,
                shape=shape,
                flip_x=sx < 0,
                flip_y=sy > 0,
                z=z,
                timestamp=float(info.get("__timestamp__", 0.0)),
                user_offset=(user[0] * 1e6, user[1] * 1e6),
                cell=cell,
                nbytes=path.stat().st_size,
                index_info=info,
            )
            t.index_mag, t.index_pixel, t.stage_xmin, t.stage_ymax = t.mag, t.pixel, t.xmin, t.ymax
            if abs(abs(sy) - abs(sx)) > 0.01 * abs(sx):
                t.notes.append("pixels are not square")
            channel = folder_channel(group)
            if channel is not None:
                t.channel, t.channel_from = channel, "folder"
            if any(exclude_pattern.search(p) for p in parts):
                t.excluded = "folder marked as bad"
            if t.kind == "video":
                t.frames = int(info.get("frames", 0))
                if t.nbytes > max_video_mb * 1e6:
                    t.excluded = f"video is {t.nbytes / 1e6:.0f} MB (more than {max_video_mb:.0f} MB)"
            tiles.append(t)
        # images that are in the directory but not in the .index (no position information)
        for f in d.iterdir():
            if re_image.match(f.name) and f.name not in listed:
                notes.append(f"{group}/{f.name}: not in .index (no position); not used")
    return tiles, notes


def load_tile(t: Tile, video_projection: str = "focus"):
    """Read the image (or the projection of a video), oriented with +x right, +y up

    A video (z-stack) is projected to one image with
    :func:`src.stack_projection.project_stack`; how is recorded in ``t.projection`` and in the
    notes. The result is stored in ``t.data``; ``t.shape`` is corrected if the file does not match
    the .index region.

    Parameters
    ----------
    t : Tile
        The tile to load.
    video_projection : {"focus", "mip", "best"}, optional
        Projection for videos (default "focus": multi-scale focus stack).
    """
    if t.kind == "image":
        a = tifffile.imread(t.path).astype(np.float32)
        if a.ndim == 3:  # a stack saved as one file: use the first frame
            a = a[0]
    else:
        m = EM.MetaArray(file=str(t.path), readAllData=False)
        t.frames = m.shape[0]

        def frame(i):
            """Frame i of the video as (rows, cols).

            Parameters
            ----------
            i : int
                Frame index.

            Returns
            -------
            numpy.ndarray
                The frame (float32), transposed from the MetaArray's [x, y] order.
            """
            return np.asarray(m[int(i)], dtype=np.float32).T  # MetaArray frames are [x, y]

        a, info = project_stack(frame, t.frames, t.pixel, method=video_projection)
        sweep = max(info.frame_sharpness) / max(np.median(info.frame_sharpness), 1e-12)
        t.projection = {
            "method": info.method,
            "frames": t.frames,
            "used": info.used,
            "dark": info.dark,
            "sharpest_frame": info.best_frame,
            "frame_sharpness": info.frame_sharpness,
            "pyramid_levels": info.levels,
        }
        t.notes.append(
            f"video: {info.summary()} ({'focus varies ' + format(sweep, '.1f') + 'x' if sweep > 1.5 else 'focus constant'})"
        )
    if a.shape != tuple(t.shape):
        t.notes.append(f"image is {a.shape}, .index region gives {t.shape}")
        t.shape = a.shape
    if t.flip_x:
        a = a[:, ::-1]
    if t.flip_y:
        a = a[::-1, :]
    t.data = np.ascontiguousarray(a)


def flatten(tiles: list):
    """Divide out the illumination pattern, estimated as the median of the (normalized) images
    with the same channel, magnification and size; needs at least 5 images.

    For fluorescence the pattern is the 30th percentile instead of the median (the background is
    the lower part of the distribution). Fluorescence images then also have their own background
    (20th percentile) subtracted. ``t.data`` is changed in place.

    Parameters
    ----------
    tiles : list of Tile
        Loaded tiles.
    """
    sets = {}
    for t in tiles:
        if not t.excluded:
            sets.setdefault((t.channel, t.mag, t.shape), []).append(t)
    for (channel, mag, shape), members in sets.items():
        if len(members) < 5:
            continue
        small = [ndi.zoom(t.data, 0.25, order=1) for t in members]
        stack = np.array([s / max(np.median(s), 1.0) for s in small])
        if channel == "BF":
            flat = np.median(stack, axis=0)
        else:  # fluorescence: the background is the lower part of the distribution
            flat = np.percentile(stack, 30, axis=0)
        flat = ndi.gaussian_filter(flat, 4.0)
        flat = flat / flat.mean()
        flat = ndi.zoom(flat, (shape[0] / flat.shape[0], shape[1] / flat.shape[1]), order=1)
        flat = np.clip(flat, 0.2, None)
        for t in members:
            t.data = t.data / flat
    for t in tiles:  # fluorescence: subtract each image's background (autofluorescence, camera offset)
        if not t.excluded and t.channel != "BF":
            t.data = t.data - np.percentile(t.data[::4, ::4], 20)


def check_objectives(tiles: list, log=print, min_ncc: float = 0.5) -> list:
    """Images of the same place, taken with the same objective, look alike at their own pixel
    size; images taken with different objectives do not. So if an image matches images in its
    folder that the .index says were taken with another objective (and matches fewer images with
    its own objective), its objective is taken to be wrong. It is then given the objective of the
    images it matches, and placed relative to the best match (not from the stage position, which
    includes the offset of the wrong objective). Returns a list of the changes.

    Parameters
    ----------
    tiles : list of Tile
        Loaded tiles. Changed tiles get a new mag, pixel, xmin and ymax, and a note.
    log : callable, optional
        Function for progress messages (default print).
    min_ncc : float, optional
        Smallest NCC (on 4x-reduced gradient features, with at least 60% overlap) that counts as a
        match (default 0.5).

    Returns
    -------
    list of dict
        One entry per changed image: image, index_mag, used_mag, matched (the best match), ncc.
    """
    changes = []
    by_folder = {}
    for t in tiles:
        by_folder.setdefault(t.group, []).append(t)
    for group, members in by_folder.items():
        if len({t.mag for t in members}) < 2:
            continue
        small = {}
        for t in members:
            a = ndi.zoom(t.data, 0.25, order=1)
            small[id(t)] = (features(a, "grad"), np.ones_like(a, dtype=np.float32))
        matches = {id(t): [] for t in members}  # (other tile, ncc, dx, dy in the other's pixels)
        for i, ta in enumerate(members):
            for tb in members[i + 1 :]:
                fa, ma = small[id(ta)]
                fb, mb = small[id(tb)]
                if fa.shape != fb.shape:
                    continue
                ncc, overlap = masked_ncc(fa, ma, fb, mb)
                h, w = fb.shape
                ncc = np.where(overlap >= 0.6 * h * w, ncc, -1.0)
                pi, pj = np.unravel_index(np.argmax(ncc), ncc.shape)
                if ncc[pi, pj] >= min_ncc:
                    # tb's top-left corner is at (r, c) in ta's (4x reduced) pixels
                    r, c = (pi - h + 1) * 4, (pj - w + 1) * 4
                    matches[id(ta)].append((tb, float(ncc[pi, pj]), c, r))
                    matches[id(tb)].append((ta, float(ncc[pi, pj]), -c, -r))
        for t in members:
            votes = {}
            for other, ncc, c, r in matches[id(t)]:
                votes[other.mag] = votes.get(other.mag, 0) + 1
            others = {m: n for m, n in votes.items() if m != t.mag}
            if not others:
                continue
            mag = max(others, key=others.get)
            if others[mag] <= votes.get(t.mag, 0):
                t.notes.append(f"objective: matches images labelled {mag:g}x and {t.mag:g}x; not changed")
                continue
            best = max((m for m in matches[id(t)] if m[0].mag == mag), key=lambda m: m[1])
            other, ncc, c, r = best
            old = t.mag
            t.mag, t.pixel = mag, other.pixel
            # other's pixel (c, r) holds t's top-left corner -> t's top-left in world coordinates
            t.xmin = other.xmin - c * other.pixel
            t.ymax = other.ymax + r * other.pixel
            t.notes.append(
                f"objective: .index says {old:g}x, but the image matches {others[mag]:d} image(s) taken at "
                f"{mag:g}x (best {other.name}, NCC {ncc:.2f}); used {mag:g}x, placed from {other.name}"
            )
            changes.append({"image": t.name, "index_mag": old, "used_mag": mag, "matched": other.name, "ncc": round(ncc, 2)})
            log(f"   {t.name}: .index says {old:g}x, matches {mag:g}x images (e.g. {other.name}, NCC {ncc:.2f})")
    return changes


# --- resampling and compositing ----------------------------------------------------------------


class Canvas:
    """A world-coordinate region sampled at res um/pixel; row 0 is the top (largest y)

    Parameters
    ----------
    xmin, xmax, ymin, ymax : float
        The region, um. The canvas covers at least this region (whole pixels, from xmin and ymax).
    res : float
        Pixel size, um.

    Attributes
    ----------
    x0, y1 : float
        World position of the left and top edges, um.
    res : float
        Pixel size, um.
    shape : tuple of int
        (rows, cols).
    """

    def __init__(self, xmin: float, xmax: float, ymin: float, ymax: float, res: float):
        self.x0, self.y1, self.res = xmin, ymax, res
        self.shape = (max(int(np.ceil((ymax - ymin) / res)), 1), max(int(np.ceil((xmax - xmin) / res)), 1))

    @property
    def bounds(self):
        """The region the canvas covers.

        Returns
        -------
        tuple of float
            (xmin, xmax, ymin, ymax), um.
        """
        return (self.x0, self.x0 + self.shape[1] * self.res, self.y1 - self.shape[0] * self.res, self.y1)

    def to_pixel(self, x: float, y: float) -> tuple:
        """(row, col) of a world point

        Parameters
        ----------
        x, y : float
            World position, um.

        Returns
        -------
        tuple of float
            (row, col), fractional.
        """
        return ((self.y1 - y) / self.res, (x - self.x0) / self.res)


_resample_cache = {}


def resampled(t: Tile, res: float, mode: str) -> tuple:
    """The tile's image (mode 'raw') or registration features ('dog' or 'grad') at res um/pixel,
    and its mask

    Images are smoothed before they are shrunk (anti-aliasing) and interpolated linearly. In the
    feature modes, a 3-pixel border is masked out (filter edge effects). Results are cached in
    ``_resample_cache`` (cleared after each slice). See "Registration methods" in the
    documentation for the constants.

    Parameters
    ----------
    t : Tile
        A loaded tile.
    res : float
        Pixel size of the result, um.
    mode : {"raw", "dog", "grad"}
        The image itself, or features for registration (see :func:`features`).

    Returns
    -------
    image : numpy.ndarray
        float32.
    mask : numpy.ndarray
        float32, 1 where the image is valid.
    """
    key = (t.serial, t.mag, round(res, 6), mode)
    if key in _resample_cache:
        return _resample_cache[key]
    zoom = t.pixel / res
    a = t.data
    if zoom < 1:
        a = ndi.gaussian_filter(a, 0.5 / zoom)
    out_shape = (max(int(round(a.shape[0] * zoom)), 1), max(int(round(a.shape[1] * zoom)), 1))
    a = ndi.zoom(a, (out_shape[0] / a.shape[0], out_shape[1] / a.shape[1]), order=1)
    mask = np.ones_like(a, dtype=np.float32)
    if mode != "raw":
        a = features(a, mode)
        border = 3
        mask[:border, :] = mask[-border:, :] = mask[:, :border] = mask[:, -border:] = 0
    _resample_cache[key] = (a.astype(np.float32), mask)
    return _resample_cache[key]


def features(a: np.ndarray, mode: str) -> np.ndarray:
    """Band-passed image for registration, robustly normalized.
    'dog': difference of Gaussians (images of the same kind);
    'grad': gradient magnitude (does not depend on the sign of the contrast)

    The image is log-compressed above its dark floor, band-passed, then standardized with the
    median and the MAD and clipped at +-6 robust SD. See "Registration methods" in the
    documentation for the constants.

    Parameters
    ----------
    a : numpy.ndarray
        2D image.
    mode : {"dog", "grad"}
        Feature type.

    Returns
    -------
    numpy.ndarray
        The features, in robust SD units, between -6 and 6.
    """
    lo = np.percentile(a, 0.5)
    a = np.log1p(np.maximum(a - lo, 0.0) / max(np.percentile(a, 50) - lo, 1.0))
    if mode == "dog":
        f = ndi.gaussian_filter(a, 1.0) - ndi.gaussian_filter(a, 10.0)
    else:
        g = ndi.gaussian_gradient_magnitude(a, 1.5)
        f = g - ndi.gaussian_filter(g, 10.0)
    med = np.median(f)
    mad = np.median(np.abs(f - med)) * 1.4826 + 1e-9  # MAD * 1.4826 is the SD estimate. 1e-9 avoids divide by zero.
    return np.clip((f - med) / mad, -6.0, 6.0)  # Gaussian noise clipping (6 SDs) is robust to outliers


def feather(shape: tuple) -> np.ndarray:
    """Blending weight: 1 in the middle, falling to ~0 at the edges

    The weight rises linearly over the outer 15% of each dimension.

    Parameters
    ----------
    shape : tuple of int
        (rows, cols).

    Returns
    -------
    numpy.ndarray
        float32 weights between 0 and 1.
    """
    r = np.minimum(np.arange(shape[0]), np.arange(shape[0])[::-1]) + 1.0
    c = np.minimum(np.arange(shape[1]), np.arange(shape[1])[::-1]) + 1.0
    w = np.minimum.outer(r / max(shape[0] * 0.15, 1.0), c / max(shape[1] * 0.15, 1.0))
    return np.clip(w, 0.0, 1.0).astype(np.float32)


def composite(tiles: list, canvas: Canvas, mode: str, placed: bool = True, blend: str = "feather") -> tuple:
    """Paste tiles into the canvas. Returns (image, weight).

    Overlapping tiles are averaged, weighted by their masks (and, with ``blend="feather"``, by a
    weight that falls toward each tile's edges, which hides the seams).

    Parameters
    ----------
    tiles : list of Tile
        Loaded tiles.
    canvas : Canvas
        Region and resolution of the result.
    mode : {"raw", "dog", "grad"}
        What to paste (see :func:`resampled`).
    placed : bool, optional
        Use the corrected positions (default True) or the stage positions.
    blend : {"feather", "flat"}, optional
        Blending of overlapping tiles (default "feather"); anything but "feather" gives equal
        weights.

    Returns
    -------
    image : numpy.ndarray
        The composite (0 where no tile covers the canvas).
    weight : numpy.ndarray
        Sum of the weights at each pixel (0 where no tile covers the canvas).
    """
    acc = np.zeros(canvas.shape, np.float32)
    wsum = np.zeros(canvas.shape, np.float32)
    for t in tiles:
        a, mask = resampled(t, canvas.res, mode)
        b = t.bounds(placed)
        r0, c0 = canvas.to_pixel(b[0], b[3])
        r0, c0 = int(round(r0)), int(round(c0))
        rs = slice(max(r0, 0), min(r0 + a.shape[0], canvas.shape[0]))
        cs = slice(max(c0, 0), min(c0 + a.shape[1], canvas.shape[1]))
        if rs.start >= rs.stop or cs.start >= cs.stop:
            continue
        ar = slice(rs.start - r0, rs.stop - r0)
        ac = slice(cs.start - c0, cs.stop - c0)
        w = mask[ar, ac] * (feather(a.shape)[ar, ac] if blend == "feather" else 1.0)
        acc[rs, cs] += a[ar, ac] * w
        wsum[rs, cs] += w
    image = np.where(wsum > 0, acc / np.maximum(wsum, 1e-9), 0.0)
    return image, wsum


# --- registration -------------------------------------------------------------------------------


def masked_ncc(fixed, fmask, moving, mmask):
    """Masked normalized cross-correlation (Padfield, 2012) for all translations. This method
    uses the Fourier transform to compute the NCC efficiently. The masks are 1 where the image
    is valid.

    Element [i, j] is for the moving image's top-left corner at (i - h + 1, j - w + 1) in the
    fixed image, where (h, w) is the moving image's shape. Returns (ncc, overlap in pixels).

    Reference: D. Padfield, "Masked Object Registration in the Fourier Domain," in IEEE
    Transactions on Image Processing, vol. 21, no. 5, pp. 2706-2718, May 2012,
    doi: 10.1109/TIP.2011.2181402.

    Parameters
    ----------
    fixed : numpy.ndarray
        The fixed image.
    fmask : numpy.ndarray
        Its mask.
    moving : numpy.ndarray
        The moving image.
    mmask : numpy.ndarray
        Its mask.

    Returns
    -------
    ncc : numpy.ndarray
        NCC for each translation, between -1 and 1 (0 where it is undefined), shape
        (fixed rows + moving rows - 1, fixed cols + moving cols - 1).
    overlap : numpy.ndarray
        Number of valid pixels that overlap at each translation.
    """
    shape = (fixed.shape[0] + moving.shape[0] - 1, fixed.shape[1] + moving.shape[1] - 1)
    fshape = [scipy.fft.next_fast_len(s, real=True) for s in shape]
    fixed = fixed * fmask
    moving = (moving * mmask)[::-1, ::-1]
    mmask_r = mmask[::-1, ::-1]

    def fft(a):
        """Real 2D FFT, padded to the common fast size.

        Parameters
        ----------
        a : numpy.ndarray
            2D array.

        Returns
        -------
        numpy.ndarray
            Its transform.
        """
        return scipy.fft.rfft2(a, fshape)

    def xc(fa, fb):
        """Cross-correlation from two transforms, cropped to the full-overlap size.

        Parameters
        ----------
        fa, fb : numpy.ndarray
            Transforms from ``fft`` (fb of a flipped image, so the product is a correlation).

        Returns
        -------
        numpy.ndarray
            The correlation for each translation.
        """
        return scipy.fft.irfft2(fa * fb, fshape)[: shape[0], : shape[1]]

    F, F2, Fm = fft(fixed), fft(fixed * fixed), fft(fmask)
    M, M2, Mm = fft(moving), fft(moving * moving), fft(mmask_r)
    overlap = np.round(xc(Fm, Mm))
    sF, sM = xc(F, Mm), xc(Fm, M)
    with np.errstate(divide="ignore", invalid="ignore"):
        num = xc(F, M) - sF * sM / overlap
        vF = xc(F2, Mm) - sF**2 / overlap
        vM = xc(Fm, M2) - sM**2 / overlap
        ncc = num / np.sqrt(np.maximum(vF * vM, 1e-12))
    ncc[~np.isfinite(ncc)] = 0.0
    return np.clip(ncc, -1.0, 1.0), overlap


@dataclass
class Match:
    """The result of :func:`register`: the best translation and how clearly it is defined.

    Attributes
    ----------
    dx, dy : float
        Correction to add to the moving tiles' positions, um.
    peak : float
        NCC at the best shift.
    second : float
        Best NCC away from the peak (outside an exclusion radius).
    z : float
        (peak - median) / robust SD of the NCC over the search window.
    overlap : float
        Fraction of the moving image that overlaps the fixed images at the peak.
    info : numpy.ndarray or None
        2x2 curvature (information matrix) of the NCC peak in (x, y), per um^2: large along
        directions in which the shift is well defined.
    """
    dx: float  # um, correction to add to the moving tiles' positions
    dy: float
    peak: float  # NCC at the best shift
    second: float  # best NCC away from the peak
    z: float  # (peak - median) / robust SD over the search window
    overlap: float  # fraction of the moving image that overlaps at the peak
    info: Optional[np.ndarray] = None  # 2x2 curvature of the NCC peak in (x, y), per um^2

    def ok(self, min_peak: float, min_z: float = 7.0, min_margin: float = 0.04) -> bool:
        """Whether the match is clear enough to use.

        Parameters
        ----------
        min_peak : float
            Smallest NCC peak.
        min_z : float, optional
            Smallest z-score of the peak (default 7).
        min_margin : float, optional
            Smallest margin of the peak over the next-best peak (default 0.04).

        Returns
        -------
        bool
            True if all three criteria are met.
        """
        return self.peak >= min_peak and self.z >= min_z and self.peak - self.second >= min_margin

    @property
    def score(self) -> float:
        """Margin of the peak over the next-best peak.

        Returns
        -------
        float
            peak - second.
        """
        return self.peak - self.second


def register(
    fixed_tiles: list,
    moving_tiles: list,
    res: float,
    radius: float,
    mode: str,
    min_overlap: float = 0.3,
    min_overlap_px: int = 0,
):
    """Find the translation of the moving tiles (as a group) that best matches the fixed tiles,
    within radius um of their current positions. The overlap must be at least min_overlap of the
    moving image and at least min_overlap_px pixels. Returns a Match or None.

    Both sides are composited on canvases at ``res`` and compared by masked NCC. A peak near the
    edge of the search window must be up to 0.1 better than one at the centre (a mild preference
    for small shifts). The peak is located to sub-pixel precision (parabola fits), and its
    curvature gives the information matrix ``Match.info``.

    Parameters
    ----------
    fixed_tiles : list of Tile
        Tiles already placed (only those near the moving tiles are used).
    moving_tiles : list of Tile
        Tiles to place, moved together.
    res : float
        Resolution of the comparison, um per pixel.
    radius : float
        Largest shift searched for, um.
    mode : {"dog", "grad"}
        Features compared (see :func:`features`).
    min_overlap : float, optional
        Smallest overlap as a fraction of the moving image (default 0.3).
    min_overlap_px : int, optional
        Smallest overlap in pixels (default 0).

    Returns
    -------
    Match or None
        None if there are no fixed tiles nearby, too few valid pixels, or too few allowed shifts.
    """
    mb = np.array([t.bounds() for t in moving_tiles])
    mcanvas = Canvas(mb[:, 0].min(), mb[:, 1].max(), mb[:, 2].min(), mb[:, 3].max(), res)
    fcanvas = Canvas(
        mcanvas.bounds[0] - radius, mcanvas.bounds[1] + radius, mcanvas.bounds[2] - radius, mcanvas.bounds[3] + radius, res
    )
    near = [t for t in fixed_tiles if overlaps(t.bounds(), fcanvas.bounds)]
    if len(near) == 0:
        return None
    m_img, m_w = composite(moving_tiles, mcanvas, mode)
    f_img, f_w = composite(near, fcanvas, mode)
    m_mask = (m_w > 0).astype(np.float32)
    f_mask = (f_w > 0).astype(np.float32)
    if m_mask.sum() < 50 or f_mask.sum() < 50:
        return None
    ncc, overlap = masked_ncc(f_img, f_mask, m_img, m_mask)
    h, w = m_img.shape
    r_nom, c_nom = fcanvas.to_pixel(mcanvas.x0, mcanvas.y1)
    i = np.arange(ncc.shape[0])[:, None] - h + 1 - r_nom  # row shift for each element
    j = np.arange(ncc.shape[1])[None, :] - w + 1 - c_nom
    rpx = radius / res
    window = (i**2 + j**2 <= rpx**2) & (overlap >= max(min_overlap * m_mask.sum(), min_overlap_px))
    if window.sum() < 10:
        return None
    # a mild preference for small shifts: a peak at the edge of the window must be 0.1 better
    vals = np.where(window, ncc - 0.1 * (i**2 + j**2) / rpx**2, -np.inf)
    pi, pj = np.unravel_index(np.argmax(vals), vals.shape)
    peak = ncc[pi, pj]
    # sub-pixel peak (parabola in each direction)
    di = dj = 0.0
    if 0 < pi < ncc.shape[0] - 1 and window[pi - 1, pj] and window[pi + 1, pj]:
        a, b, c = ncc[pi - 1, pj], peak, ncc[pi + 1, pj]
        di = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) < 0 else 0.0
    if 0 < pj < ncc.shape[1] - 1 and window[pi, pj - 1] and window[pi, pj + 1]:
        a, b, c = ncc[pi, pj - 1], peak, ncc[pi, pj + 1]
        dj = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) < 0 else 0.0
    # how sharply the peak is defined in each direction: along a straight feature (a string
    # running along the overlap) the NCC hardly changes, and that direction gets little weight
    step = 2
    info = None
    if step <= pi < ncc.shape[0] - step and step <= pj < ncc.shape[1] - step:
        f = ncc
        hjj = (f[pi, pj + step] - 2 * peak + f[pi, pj - step]) / step**2
        hii = (f[pi + step, pj] - 2 * peak + f[pi - step, pj]) / step**2
        hij = (f[pi + step, pj + step] - f[pi + step, pj - step] - f[pi - step, pj + step] + f[pi - step, pj - step]) / (
            4 * step**2
        )
        # rows run toward -y, so the cross term changes sign in (x, y)
        hess = np.array([[hjj, -hij], [-hij, hii]]) / res**2
        ev, vec = np.linalg.eigh(-hess)
        ev = np.clip(ev, 1e-4 * max(ev.max(), 1e-12), None)
        info = (vec * ev) @ vec.T
    exclusion = max(4.0, 0.08 * min(h, w))
    away = window & (((i - i[pi, 0]) ** 2 + (j - j[0, pj]) ** 2) > exclusion**2)
    second = ncc[away].max() if away.any() else -1.0
    wv = ncc[window]
    med = np.median(wv)
    mad = np.median(np.abs(wv - med)) * 1.4826 + 1e-6
    return Match(
        dx=(j[0, pj] + dj) * res,
        dy=-(i[pi, 0] + di) * res,
        peak=float(peak),
        second=float(second),
        z=float((peak - med) / mad),
        overlap=float(overlap[pi, pj] / m_mask.sum()),
        info=info,
    )


def overlaps(a: tuple, b: tuple) -> bool:
    """Whether two rectangles overlap (touching edges do not count).

    Parameters
    ----------
    a, b : tuple of float
        Rectangles as (xmin, xmax, ymin, ymax).

    Returns
    -------
    bool
        True if their interiors intersect.
    """
    return a[0] < b[1] and b[0] < a[1] and a[2] < b[3] and b[2] < a[3]


def overlap_area(a: tuple, b: tuple) -> float:
    """Area of the intersection of two rectangles.

    Parameters
    ----------
    a, b : tuple of float
        Rectangles as (xmin, xmax, ymin, ymax).

    Returns
    -------
    float
        The area (0 if they do not overlap), in the units of the rectangles squared.
    """
    return max(0.0, min(a[1], b[1]) - max(a[0], b[0])) * max(0.0, min(a[3], b[3]) - max(a[2], b[2]))


def stitch(tiles: list, radius: float = 200.0, log=print) -> dict:
    """Register neighbouring pairs of tiles and solve for all of the translations together,
    starting from the tiles' current placements.

    The tiles usually overlap by only a strip 40-200 um wide along one edge (a few % of the
    area), and the stage error can be larger than the overlap, so every pair that shares an
    edge (or would, within the search radius) is tried, and a match needs only an overlap of
    min_strip um x 1/4 of the tile edge. Pairs of the same channel are compared by band-passed
    intensity; pairs of different channels (fluorescence and brightfield) by gradient magnitude.

    The pairs are combined by :func:`solve_pairs`; the most inconsistent pair is dropped and the
    solution repeated until every pair agrees within max(2 pixels, 3 x the median residual). Tiles
    that are connected to the others by at least one pair get the solved translation and placement
    "stitched"; the others are left as they were.

    Parameters
    ----------
    tiles : list of Tile
        Loaded tiles; their dx, dy and placement are changed.
    radius : float, optional
        Largest correction searched for in each pair, um (default 200).
    log : callable, optional
        Function for progress messages (default print).

    Returns
    -------
    dict
        Report: tiles, pairs_tried, pairs_registered, pairs_rejected, pairs_used, unconnected
        (names), residual_um (median residual along each pair's best-defined direction).
    """
    res = max(t.pixel for t in tiles)
    min_strip = 25.0  # um
    d0 = np.array([(t.dx, t.dy) for t in tiles])
    pairs = []
    tried = 0
    for a in range(len(tiles)):
        for b in range(a + 1, len(tiles)):
            A, B = tiles[a].bounds(), tiles[b].bounds()
            ox = min(A[1], B[1]) - max(A[0], B[0])  # overlap in x, um (negative: gap)
            oy = min(A[3], B[3]) - max(A[2], B[2])
            edge = min(tiles[a].width, tiles[a].height, tiles[b].width, tiles[b].height)
            # neighbours: a shared edge at least a third of a tile long, and overlapping (or
            # separated by less than half the search radius) in the other direction
            if not ((ox > edge / 3 and oy > -radius / 2) or (oy > edge / 3 and ox > -radius / 2)):
                continue
            tried += 1
            same = tiles[a].channel == tiles[b].channel
            mode = "dog" if same else "grad"
            min_px = int(min_strip * edge / 4 / res**2)
            m = register([tiles[a]], [tiles[b]], res, radius, mode, min_overlap=0.0, min_overlap_px=min_px)
            # small windows and smooth features make the z-score small even for clear peaks,
            # so for pairs the margin over the next-best peak is what matters
            min_peak, min_margin = (0.25, 0.08) if same else (0.15, 0.05)
            if m is not None and m.info is not None and m.ok(min_peak=min_peak, min_z=0.0, min_margin=min_margin):
                # measured: (d_b - d_a) should be the current difference plus the correction
                pairs.append([a, b, d0[b] - d0[a] + np.array([m.dx, m.dy]), m.info])
    n = len(tiles)
    report = {"tiles": n, "pairs_tried": tried, "pairs_registered": len(pairs), "pairs_rejected": 0}
    if len(pairs) == 0:
        return report
    for iteration in range(3 * n):
        sol, whitened = solve_pairs(pairs, n, d0)
        limit = max(2.0 * res, 3.0 * np.median(whitened))
        if not (whitened > limit).any():
            break
        # drop the worst pair only, then solve again (one bad pair spreads error to its neighbours)
        worst = int(np.argmax(whitened))
        report["pairs_rejected"] += 1
        pairs = pairs[:worst] + pairs[worst + 1 :]
        if len(pairs) == 0:
            sol = d0
            whitened = np.zeros(0)
            break
    connected = set()
    for a, b, *_ in pairs:
        connected.update((a, b))
    for k, t in enumerate(tiles):
        if k in connected:
            t.dx, t.dy = float(sol[k, 0]), float(sol[k, 1])
            t.placement = "stitched"
    report["pairs_used"] = len(pairs)
    report["unconnected"] = [t.name for k, t in enumerate(tiles) if k not in connected]
    report["residual_um"] = round(float(np.median(whitened)), 2) if len(pairs) else None
    return report


def solve_pairs(pairs: list, n: int, d0: np.ndarray, prior: float = 0.01) -> tuple:
    """Least squares for the tile translations d (n x 2): each pair (a, b, m, W) contributes
    (d_b - d_a - m)' W (d_b - d_a - m), where W is the pair's 2x2 information matrix, and each
    tile a weak pull toward its starting translation d0. Returns (d, the residual of each pair
    along its best-defined direction, um).

    Parameters
    ----------
    pairs : list
        [a, b, m, W]: tile indices, measured d_b - d_a (um, length 2) and information matrix.
    n : int
        Number of tiles.
    d0 : numpy.ndarray
        Starting translations (n x 2), um.
    prior : float, optional
        Weight of the pull toward d0, relative to a typical pair (default 0.01).

    Returns
    -------
    d : numpy.ndarray
        Translations (n x 2), um.
    residuals : numpy.ndarray
        For each pair, the residual along its best-defined direction, um.
    """
    scale = np.median([np.trace(W) for *_, W in pairs]) / 2.0
    rows, rhs = [], []
    for a, b, m, W in pairs:
        L = np.linalg.cholesky(W / scale + 1e-9 * np.eye(2))
        block = np.zeros((2, 2 * n))
        block[:, 2 * b : 2 * b + 2] = L.T
        block[:, 2 * a : 2 * a + 2] = -L.T
        rows.append(block)
        rhs.append(L.T @ m)
    rows.append(np.sqrt(prior) * np.eye(2 * n))
    rhs.append(np.sqrt(prior) * np.asarray(d0, dtype=float).ravel())
    sol, *_ = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)
    sol = sol.reshape(n, 2)
    whitened = []
    for a, b, m, W in pairs:
        r = sol[b] - sol[a] - m
        ev, vec = np.linalg.eigh(W)
        whitened.append(abs(vec[:, -1] @ r))  # along the best-defined direction
    return sol, np.array(whitened)


def hysteresis(base: list) -> Optional[dict]:
    """Do the stitching corrections depend on the direction the stage last moved?
    (Least-squares fit of correction = a + h * sign(step), for x and y.)

    Parameters
    ----------
    base : list of Tile
        The stitched base tiles (at least 6 stitched tiles are needed).

    Returns
    -------
    dict or None
        For each axis with enough moves in both directions: half_backlash_um (h) and n; None if
        there are too few tiles.
    """
    ordered = sorted((t for t in base if t.placement == "stitched"), key=lambda t: t.timestamp)
    if len(ordered) < 6:
        return None
    out = {}
    for axis, pos, corr in (("x", "xmin", "dx"), ("y", "ymax", "dy")):
        steps = np.sign(np.diff([getattr(t, pos) for t in ordered]))
        c = np.array([getattr(t, corr) for t in ordered[1:]])
        use = steps != 0
        if use.sum() < 4 or len(set(steps[use])) < 2:
            continue
        A = np.column_stack([np.ones(use.sum()), steps[use]])
        (offset, h), *_ = np.linalg.lstsq(A, c[use], rcond=None)
        out[axis] = {"half_backlash_um": round(float(h), 1), "n": int(use.sum())}
    return out or None


def register_group(
    placed: list, group: list, base_res: float, rig_mags: list, log=print, decision: Optional[str] = None
) -> dict:
    """Place a group of tiles (one folder, channel and magnification) against the placed tiles

    The group is first registered as a whole: in a large window (3x the group's extent, 300 to
    1500 um) at a coarse resolution, then in a small window at a finer one. If that fails and the
    rig has other objectives, the group is tried at each of their magnifications (the objective
    setting is sometimes not changed when the objective is). If the group is placed, each tile is
    then refined on its own within 30 um. A match that would move the group more than
    ``max_move_um`` (250 um) is not used and is reported (``refused_move_um``): such moves are
    usually the slice having been moved between images. Groups that cannot be placed stay at
    the stage position, but their images are stitched to each other (placement
    "stage+stitched"), so that, for example, the overlapping z-stacks of a filled cell are
    consistent with each other.

    Parameters
    ----------
    placed : list of Tile
        Tiles placed so far.
    group : list of Tile
        The tiles to place; their dx, dy, placement, quality (and possibly mag) are changed.
    base_res : float
        Pixel size of the base tiles, um.
    rig_mags : list of float
        Magnifications used on this rig (from all of the slice's tiles).
    log : callable, optional
        Function for progress messages (default print).
    decision : str, optional
        Decision from the exceptions workbook for a refused move: "stage" (default), "accept",
        "exclude" or "manual".

    Returns
    -------
    dict
        Report: group, channel, mag, n, mode, match (see :func:`summary`), placement, and
        magnification_corrected and refined when they apply, and within_group (the report of
        :func:`stitch`) when a group that could not be placed was aligned within itself.
    """
    g0 = group[0]
    same_kind = all(t.channel == "BF" and not t.mixed for t in group + placed)
    mode = "dog" if same_kind else "grad"
    min_peak = 0.25 if mode == "dog" else 0.15
    result = {"group": g0.group, "channel": g0.channel, "mag": g0.mag, "n": len(group), "mode": mode}

    def attempt(tiles: list) -> Optional[Match]:
        """Register the tiles as a group: coarse, then fine.

        Parameters
        ----------
        tiles : list of Tile
            The group (or a copy at another magnification).

        Returns
        -------
        Match or None
            The fine match if it is clear, otherwise the coarse match (which may not be ok).
        """
        pixel = tiles[0].pixel
        near = [t for t in placed if overlaps(t.bounds(), expand(group_bounds(tiles), 1500.0))]
        finest = min((t.pixel for t in near), default=base_res)
        fine_res = max(pixel * 2.0, finest, base_res if pixel >= base_res else finest)
        extent = max(group_bounds(tiles)[1] - group_bounds(tiles)[0], group_bounds(tiles)[3] - group_bounds(tiles)[2])
        radius = float(np.clip(3.0 * extent, 300.0, 1500.0))
        coarse_res = max(fine_res, radius / 150.0)
        m = register(placed, tiles, coarse_res, radius, mode)
        if m is None or not m.ok(min_peak):
            return m
        shift_tiles(tiles, m.dx, m.dy)
        m2 = register(placed, tiles, fine_res, 4.0 * coarse_res, mode)
        shift_tiles(tiles, -m.dx, -m.dy)
        if m2 is not None and m2.ok(min_peak, min_z=5.0, min_margin=0.0):
            m2.dx += m.dx
            m2.dy += m.dy
            return m2
        return m

    def usable(m: Optional[Match]) -> bool:
        """A clear match that does not move the group more than max_move_um."""
        return m is not None and m.ok(min_peak) and np.hypot(m.dx, m.dy) <= max_move_um

    def too_far(m: Optional[Match]) -> bool:
        """A clear match that would move the group more than max_move_um."""
        return m is not None and m.ok(min_peak) and np.hypot(m.dx, m.dy) > max_move_um

    match = attempt(group)
    result["match"] = summary(match)
    refused = [match] if too_far(match) else []
    if not usable(match) and len(rig_mags) > 1:
        # the objective may not be what the .index says
        for mag in rig_mags:
            if mag == g0.mag:
                continue
            alt = [t.with_mag(mag) for t in group]
            m = attempt(alt)
            if too_far(m):
                refused.append(m)
            if usable(m) and (match is None or not usable(match) or m.peak > match.peak + 0.05):
                log(f"      {g0.group} ({g0.channel}): matches as {mag:g}x, not {g0.mag:g}x (NCC {m.peak:.2f})")
                for t, a in zip(group, alt):
                    t.notes.append(f"objective: .index says {t.mag:g}x, images match as {mag:g}x")
                    t.mag, t.pixel, t.xmin, t.ymax = a.mag, a.pixel, a.xmin, a.ymax
                result["magnification_corrected"] = [g0.objective, mag]
                match = m
                result["match"] = summary(match)
                break
    if refused and not usable(match):
        far = max(refused, key=lambda m: m.peak)
        result["refused_move_um"] = [round(far.dx, 1), round(far.dy, 1)]
        result["refused_ncc"] = round(far.peak, 3)
        result["decision"] = decision or "stage"
    if refused and not usable(match) and decision == "accept":
        match = max(refused, key=lambda m: m.peak)
        for t in group:
            t.notes.append(f"move of {np.hypot(match.dx, match.dy):.0f} um accepted in the decision table")
        result["accepted_move_um"] = [round(match.dx, 1), round(match.dy, 1)]
    elif refused and not usable(match) and decision == "exclude":
        for t in group:
            t.excluded = "group excluded in the decision table (refused move)"
        result["placement"] = "excluded"
        return result
    if not usable(match) and "accepted_move_um" not in result:
        for t in group:
            t.placement = "stage"
            if refused:
                far = max(refused, key=lambda m: m.peak)
                t.notes.append(
                    f"group registration gave a move of {np.hypot(far.dx, far.dy):.0f} um (NCC {far.peak:.2f}), more than "
                    f"{max_move_um:.0f} um (the slice was probably moved between images): not used, left at the stage position"
                )
            else:
                t.notes.append("group registration not confident: left at the stage position")
            if decision == "manual":
                t.notes.append("manual mosaic requested in the decision table")
        if refused:
            far = max(refused, key=lambda m: m.peak)
            log(f"      {g0.group} ({g0.channel} {g0.mag:g}x): move of {np.hypot(far.dx, far.dy):.0f} um refused (> {max_move_um:.0f} um)")
        result["placement"] = "stage"
        if len(group) > 1:
            # the group cannot be placed on the slice (e.g. 40x stacks of a filled cell on a 4x
            # base), but its images can still be made consistent with each other
            result["within_group"] = stitch(group, radius=100.0, log=log)
            for t in group:
                if t.placement == "stitched":
                    t.placement = "stage+stitched"
                    t.notes.append("aligned with the other images of its group (the group as a whole is at the stage position)")
        return result
    for t in group:
        t.dx, t.dy, t.placement, t.quality = match.dx, match.dy, "group", match.peak
    result["placement"] = "group"
    # refine each tile in a small window
    if len(group) > 1:
        res = max(group[0].pixel * 2.0, base_res if group[0].pixel >= base_res else group[0].pixel * 2.0)
        refined = 0
        for t in group:
            m = register(placed, [t], res, 30.0, mode)
            if m is not None and m.ok(min_peak, min_z=6.0) and np.hypot(m.dx, m.dy) < 30.0:
                t.dx += m.dx
                t.dy += m.dy
                t.placement, t.quality = "refined", m.peak
                refined += 1
        result["refined"] = refined
    return result


def summary(m: Optional[Match]) -> Optional[dict]:
    """A Match as a dictionary for the report.

    Parameters
    ----------
    m : Match or None
        The match.

    Returns
    -------
    dict or None
        Its fields (except info), rounded to 3 decimals; None if there is no match.
    """
    if m is None:
        return None
    return {k: round(v, 3) for k, v in m.__dict__.items() if k != "info"}


def group_bounds(tiles: list) -> tuple:
    """The rectangle that holds all of the tiles (at their corrected positions).

    Parameters
    ----------
    tiles : list of Tile
        The tiles.

    Returns
    -------
    tuple of float
        (xmin, xmax, ymin, ymax), um.
    """
    b = np.array([t.bounds() for t in tiles])
    return (b[:, 0].min(), b[:, 1].max(), b[:, 2].min(), b[:, 3].max())


def expand(b: tuple, d: float) -> tuple:
    """A rectangle grown by the same distance on every side.

    Parameters
    ----------
    b : tuple of float
        (xmin, xmax, ymin, ymax).
    d : float
        Distance to add on each side.

    Returns
    -------
    tuple of float
        The grown rectangle.
    """
    return (b[0] - d, b[1] + d, b[2] - d, b[3] + d)


def shift_tiles(tiles: list, dx: float, dy: float):
    """Add a translation to the corrections of the tiles.

    Parameters
    ----------
    tiles : list of Tile
        The tiles (changed in place).
    dx, dy : float
        Translation, um.
    """
    for t in tiles:
        t.dx += dx
        t.dy += dy


# --- the slice ----------------------------------------------------------------------------------


# cell markers: positions placed by hand in a Mosaic Editor file (magenta circles) and positions
# inferred here (blue squares: the centre of the cell's highest-magnification image). Saturated
# colours: white markers disappear on grey images.
annotated_style = {"marker": "o", "color": "magenta"}
inferred_style = {"marker": "s", "color": "deepskyblue"}

# parts of an acq4 Mosaic Editor file that are not about the images (as the editor saves them)
_terminal = {"multi": False, "optional": False, "renamable": False, "removable": False, "multiable": False}
_mosaic_filter_state = {
    "flowchart": {
        "pos": [0.0, 0.0],
        "bypass": False,
        "terminals": {"dataIn": {"io": "in", **_terminal}, "dataOut": {"io": "out", **_terminal}},
        "nodes": [],
        "connects": [],
        "inputNode": {"pos": [-150.0, 0.0], "bypass": False, "terminals": {"dataIn": {"io": "out", **_terminal}}},
        "outputNode": {"pos": [300.0, 0.0], "bypass": False, "terminals": {"dataOut": {"io": "in", **_terminal}}},
    }
}
_mosaic_view_state = {
    "targetRange": [[0.0, 1.0], [0.0, 1.0]],
    "viewRange": [[0.0, 1.0], [0.0, 1.0]],
    "yInverted": False,
    "xInverted": False,
    "aspectLocked": 1,
    "autoRange": [False, False],
    "autoPan": [False, False],
    "autoVisibleOnly": [False, False],
    "linkedViews": [None, None],
    "defaultPadding": 0.02,
    "mouseEnabled": [True, True],
    "mouseMode": 3,
    "enableMenu": True,
    "wheelScaleFactor": -0.125,
    "background": None,
    "logMode": [False, False],
    "limits": {"xLimits": [-1e307, 1e307], "yLimits": [-1e307, 1e307], "xRange": [None, None], "yRange": [None, None]},
}


def cell_style(c: dict) -> dict:
    """Marker style for a cell: placed by hand (annotated) or inferred.

    Parameters
    ----------
    c : dict
        A cell entry from ``SliceMosaic.cells`` (with "source").

    Returns
    -------
    dict
        ``annotated_style`` or ``inferred_style`` (marker and colour).
    """
    return annotated_style if c["source"].startswith("annotated") else inferred_style


max_gap_s = 600.0  # images in one folder more than this far apart in time are registered separately


class SliceMosaic:
    """The images of one slice, their placement, cells, figures and reports.

    Typical use::

        mosaic = SliceMosaic(slice_dir).build()
        mosaic.save_overview(path, title)

    Parameters
    ----------
    slice_dir : str or pathlib.Path
        The slice directory (``.../yyyy.mm.dd_nnn/slice_nnn``).
    max_video_mb : float, optional
        Videos larger than this (MB) are excluded (default 500).
    video_projection : {"focus", "mip", "best"}, optional
        How videos are projected to one image (default "focus"; see :mod:`src.stack_projection`).
    decisions : Decisions, optional
        Decisions from an exceptions workbook (see :mod:`src.mosaic_exceptions`).
    auto_kind : bool, optional
        Let confident brightfield/fluorescence verdicts (10x and above) override the folder
        name (default True).
    log : callable, optional
        Function for progress messages (default print).

    Attributes
    ----------
    tiles : list of Tile
        All images and videos found, including excluded ones.
    notes : list of str
        Problems found while scanning.
    report : dict
        Everything the build found and did, for the JSON report.
    cells : dict
        Cell name -> {"x", "y", "source"} (um), from :meth:`find_cells`.
    base : list of Tile
        The stitched base tiles (after :meth:`build`).
    base_res : float
        Pixel size of the base tiles, um (after :meth:`build`).
    """
    def __init__(
        self,
        slice_dir: Union[str, Path],
        max_video_mb: float = 500.0,
        video_projection: str = "focus",
        decisions: Optional[Decisions] = None,
        auto_kind: bool = True,
        log=print,
    ):
        """Scan the slice directory for images and videos (nothing is read yet).

        Parameters
        ----------
        slice_dir : str or pathlib.Path
            The slice directory.
        max_video_mb : float, optional
            Videos larger than this (MB) are excluded (default 500).
        video_projection : {"focus", "mip", "best"}, optional
            How videos are projected to one image (default "focus").
        decisions : Decisions, optional
            Decisions from an exceptions workbook (see :mod:`src.mosaic_exceptions`).
        auto_kind : bool, optional
            Let confident brightfield/fluorescence verdicts of the classifier (images at 10x and
            above) override the folder name (default True).
        log : callable, optional
            Function for progress messages (default print).
        """
        self.slice_dir = Path(slice_dir)
        self.log = log
        self.video_projection = video_projection
        self.decisions = decisions if decisions is not None else Decisions()
        self.auto_kind = auto_kind
        self.title = slice_title(self.slice_dir)
        self.tiles, self.notes = scan_slice(self.slice_dir, max_video_mb=max_video_mb)
        self.report = {"slice_dir": str(self.slice_dir), "notes": self.notes, "groups": []}
        self.cells = {}
        self.apply_kind_decisions()

    def apply_kind_decisions(self):
        """Apply the image-kind decisions of the decision table (channel, mixed, exclude)."""
        for t in self.tiles:
            d = self.decisions.get("image kind", slice=self.title, image=t.name)
            if d is None:
                continue
            if d == "exclude":
                t.excluded = "excluded in the decision table"
            elif d == "BF":
                t.channel = "BF"
            elif d == "mixed":
                t.channel, t.mixed = "BF", True
            elif d == "fluorescence":
                t.channel = t.channel if t.channel != "BF" else "yellow"
            else:
                t.channel = d
            t.channel_from = "decision"
            t.notes.append(f"decision table: {d}")

    @property
    def used(self) -> list:
        """The tiles that are not excluded.

        Returns
        -------
        list of Tile
            Tiles with an empty ``excluded``.
        """
        return [t for t in self.tiles if not t.excluded]

    def build(self):
        """Load, flatten and place all of the tiles, then find the cells.

        Steps: load the images and video projections; divide out the illumination pattern; check the
        objectives; stitch the base (the slice-level folder with the most brightfield images at the
        lowest magnification); register the other groups, lowest magnification first; stitch all of
        the slice-level low-magnification images together; find the cells; compare with the manual
        corrections. Results go to the tiles and to ``report``.

        Returns
        -------
        SliceMosaic
            self, for chaining.
        """
        t0 = time.time()
        self.processed = time.strftime("%Y-%m-%d %H:%M:%S")
        for t in self.used:
            try:
                load_tile(t, video_projection=self.video_projection)
            except Exception as e:
                t.excluded = f"could not be read: {e!r}"[:200]
        tiles = self.used
        if len(tiles) == 0:
            self.report["error"] = "no usable images"
            return self
        for t in tiles:  # brightfield or fluorescence, by appearance (before any correction)
            m = measure(t.data, t.pixel)
            kind, p = classify(m)
            t.appearance = {"rule": kind, "p_fluorescence": p, **m.as_dict()}
            # a confident verdict overrides the folder name (the classifier was trained at 10x and
            # above); decisions in the decision table override both
            if not self.auto_kind or t.channel_from == "decision" or not t.mag >= 10:
                continue
            if kind == "fluorescence" and t.channel == "BF":
                t.channel, t.channel_from = "yellow", "appearance"
                t.notes.append(f"appearance: fluorescence (p {p:.2f}), although the folder name says brightfield")
            elif kind == "BF" and t.channel != "BF":
                t.notes.append(f"appearance: brightfield (p fluorescence {p:.2f}), although the folder name says {t.channel}")
                t.channel, t.channel_from = "BF", "appearance"
        flatten(tiles)
        self.report["objectives_changed"] = check_objectives(tiles, log=self.log)
        # the base: the slice-level folder with the most brightfield images at the lowest
        # magnification; if that has fewer than 3, the folder with the most images of any channel
        low = min(t.mag for t in tiles)
        candidates = {}
        for t in tiles:
            if t.mag == low and t.cell is None:
                candidates.setdefault((t.group, t.channel), []).append(t)
        if not candidates:
            candidates = {(t.group, t.channel): [] for t in tiles if t.mag == low}
            for t in tiles:
                if t.mag == low:
                    candidates[(t.group, t.channel)].append(t)
        bf = {k: v for k, v in candidates.items() if k[1] == "BF" and len(v) >= 3}
        pool = bf if bf else candidates
        base_key = max(pool, key=lambda k: (len(pool[k]), k[0] == "."))
        base_group = base_key[0]
        base = pool[base_key]
        self.base_res = base[0].pixel
        self.log(f"   base: {len(base):d} {base[0].channel} images at {low:g}x in {base_group}")
        self.report["base"] = {"group": base_group, "mag": low, "channel": base[0].channel}
        self.base = base
        self.report["base"].update(stitch(base, log=self.log))
        self.report["base"]["hysteresis"] = hysteresis(base)
        placed = list(base)
        # everything else, a group at a time: lowest magnification first, big folders first
        groups = {}
        last = {}
        for t in sorted(tiles, key=lambda t: t.timestamp):
            if t in base:
                continue
            key = (t.group, t.channel, t.mag)
            run = last.get(key, (0, t.timestamp))
            if t.timestamp - run[1] > max_gap_s:  # a new run of images in the same folder
                run = (run[0] + 1, t.timestamp)
            last[key] = (run[0], t.timestamp)
            groups.setdefault(key + (run[0],), []).append(t)
        rig_mags = sorted({t.mag for t in self.tiles if np.isfinite(t.mag)})
        for key in sorted(groups, key=lambda k: (k[2], k[0] != ".", -len(groups[k]))):
            decision = self.decisions.get("refused matches", slice=self.title, group=key[0], channel=key[1], mag=key[2], run=key[3])
            result = register_group(placed, groups[key], self.base_res, rig_mags, log=self.log, decision=decision)
            result["run"] = key[3]
            self.log(
                f"   {key[0]:s}{'' if key[3] == 0 else ' run ' + str(key[3] + 1)} {key[1]:s} {key[2]:g}x ({len(groups[key]):d}): {result['placement']:s}"
                + (
                    f" offset ({groups[key][0].dx:+.0f}, {groups[key][0].dy:+.0f}) um, NCC {result['match']['peak']:.2f}"
                    if result["placement"] == "group"
                    else ""
                )
            )
            self.report["groups"].append(result)
            if result["placement"] not in ("stage", "excluded"):
                placed.extend(groups[key])
        # finally, all of the slice-level images at the base magnification together: the group
        # offsets above are only a starting point, as the tiles in a group are not placed
        # consistently with each other (each was taken at a different stage position)
        joint = [t for t in tiles if t.cell is None and t.mag <= low * 2.6 and not t.excluded]
        if len(joint) > len(base):
            self.report["joint"] = stitch(joint, radius=150.0, log=self.log)
            self.log(
                f"   all {len(joint):d} slice-level {low:g}x images together: {self.report['joint']['pairs_used']:d} pairs, "
                f"residual {self.report['joint']['residual_um']} um, {len(self.report['joint']['unconnected']):d} not connected"
            )
        self.find_cells()
        self.compare_user()
        self.report["seconds"] = round(time.time() - t0, 1)
        return self

    def find_cells(self):
        """Cell positions: from a Mosaic Editor file in the slice directory if there is one,
        otherwise the centre of the cell's highest-magnification image

        Positions in a Mosaic Editor file are used as they are (they were placed by hand on the
        images as aligned in the editor). Results go to ``cells`` and to ``report["cells"]``.
        """
        annotated = {}
        for mf in sorted(self.slice_dir.glob("*.mosaic")):
            try:
                items = json.load(open(mf))["items"]
            except Exception:
                continue
            for item in items:
                if item.get("type") == "CellCanvasItem":
                    pos = item["userTransform"]["pos"]
                    annotated[item["name"]] = (pos[0] * 1e6, pos[1] * 1e6, mf.name)
        cells = sorted({t.cell for t in self.tiles if t.cell is not None} | set(annotated))
        cells += sorted(d.name for d in self.slice_dir.iterdir() if d.is_dir() and re_cell.match(d.name) and d.name not in cells)
        for cell in cells:
            own = [t for t in self.used if t.cell == cell]
            decision = self.decisions.get("annotation mismatches", slice=self.title, cell=cell)
            if decision == "exclude":
                self.cells[cell] = {"x": None, "y": None, "source": "excluded in the decision table"}
                if cell in annotated:
                    self.cells[cell]["annotated"] = annotated[cell]
            elif cell in annotated and decision != "inferred":
                x, y, source = annotated[cell]
                self.cells[cell] = {"x": x, "y": y, "source": f"annotated in {source}"}
            elif len(own) > 0:
                top = max(own, key=lambda t: (t.mag, t.timestamp))
                x, y = top.center()
                self.cells[cell] = {
                    "x": x,
                    "y": y,
                    "source": f"centre of {top.name} ({top.mag:g}x)" + (" (decision table)" if decision == "inferred" else ""),
                }
                if cell in annotated:
                    self.cells[cell]["annotated"] = annotated[cell]
            else:
                self.cells[cell] = {"x": None, "y": None, "source": "no images and no annotation"}
        self.report["cells"] = self.cells

    def compare_user(self):
        """Compare the corrections with the manual ones (userTransform) in the .index files,
        relative to the base (the manual ones are relative to the stage positions, ours float)

        For each folder with a manual correction, the median manual and automatic corrections and the
        median distance between them go to ``report["manual_vs_automatic"]``.
        """
        base = getattr(self, "base", [])
        if not base:
            return
        b_dx = np.median([t.dx for t in base])
        b_dy = np.median([t.dy for t in base])
        u_bx = np.median([t.user_offset[0] for t in base])
        u_by = np.median([t.user_offset[1] for t in base])
        rows = []
        groups = {}
        for t in self.used:
            groups.setdefault(t.group, []).append(t)
        for g, members in groups.items():
            user = np.array([(t.user_offset[0] - u_bx, t.user_offset[1] - u_by) for t in members])
            ours = np.array([(t.dx - b_dx, t.dy - b_dy) for t in members])
            if np.abs(user).max() < 1.0:
                continue
            rows.append(
                {
                    "group": g,
                    "n": len(members),
                    "manual_um": [round(float(v), 1) for v in np.median(user, axis=0)],
                    "automatic_um": [round(float(v), 1) for v in np.median(ours, axis=0)],
                    "difference_um": round(float(np.median(np.hypot(*(ours - user).T))), 1),
                    "placement": sorted({t.placement for t in members}),
                }
            )
        self.report["manual_vs_automatic"] = rows

    # --- rendering ---

    def channel_ranges(self) -> dict:
        """Display range of each channel, from all of its tiles.

        Returns
        -------
        dict
            Channel -> (lo, hi): the 0.5-99.8 percentiles for brightfield; the 90-99.95 percentiles
            for fluorescence (to show the brightest structures, not the autofluorescence).
        """
        ranges = {}
        for channel in {t.channel for t in self.used}:
            sample = np.concatenate([t.data[::8, ::8].ravel() for t in self.used if t.channel == channel])
            if channel == "BF":
                ranges[channel] = tuple(np.percentile(sample, [0.5, 99.8]))
            else:  # show the brightest structures, not the autofluorescence
                lo, hi = np.percentile(sample, [90, 99.95])
                ranges[channel] = (lo, max(hi, lo + 1.0))
        return ranges

    def render(
        self, canvas: Canvas, max_mag: float, placed: bool = True, ranges: Optional[dict] = None, cells: bool = True
    ) -> np.ndarray:
        """RGB image: brightfield in grey, fluorescence as transparent overlays.
        cells: include the images in the cell folders

        Brightfield magnifications are composited separately and laid over each other, higher
        magnifications on top. Each fluorescence channel is laid over in its colour, with an opacity
        that follows its intensity.

        Parameters
        ----------
        canvas : Canvas
            Region and resolution.
        max_mag : float
            Leave out tiles at higher magnifications.
        placed : bool, optional
            Corrected positions (default True) or stage positions.
        ranges : dict, optional
            Display ranges (default :meth:`channel_ranges`).
        cells : bool, optional
            Include the images in the cell folders (default True).

        Returns
        -------
        numpy.ndarray
            (rows, cols, 3) RGB image, values 0-1.
        """
        ranges = ranges or self.channel_ranges()
        tiles = sorted((t for t in self.used if t.mag <= max_mag and (cells or t.cell is None)), key=lambda t: t.mag)
        rgb = np.zeros(canvas.shape + (3,), np.float32)
        bf = [t for t in tiles if t.channel == "BF"]
        # higher magnifications on top: composite each magnification, then lay them over
        for mag in sorted({t.mag for t in bf}):
            img, w = composite([t for t in bf if t.mag == mag], canvas, "raw", placed)
            lo, hi = ranges["BF"]
            g = np.clip((img - lo) / (hi - lo), 0, 1)
            cover = (w > 0)[..., None]
            rgb = np.where(cover, g[..., None], rgb)
        for channel, color in channel_colors.items():
            ft = [t for t in tiles if t.channel == channel]
            if not ft:
                continue
            img, w = composite(ft, canvas, "raw", placed, blend="feather")
            lo, hi = ranges[channel]
            alpha = np.clip((img - lo) / (hi - lo), 0, 1) ** 0.8 * 0.9 * (w > 0)
            rgb = rgb * (1 - alpha[..., None]) + np.array(color)[None, None, :] * alpha[..., None]
        return np.clip(rgb, 0, 1)

    def overview_canvas(self, max_pixels: int = 2400) -> Canvas:
        """A canvas covering the slice-level low-magnification images.

        Parameters
        ----------
        max_pixels : int, optional
            Largest size of the canvas, pixels (default 2400); the resolution is never finer than
            the base.

        Returns
        -------
        Canvas
            The canvas, with a 50 um margin.
        """
        low = min(t.mag for t in self.used)
        b = group_bounds([t for t in self.used if t.mag <= low * 2.6 and t.cell is None] or self.used)
        res = max(self.base_res, max(b[1] - b[0], b[3] - b[2]) / max_pixels)
        return Canvas(*expand(b, 50.0), res)

    def save_overview(self, filename: Path, title: str):
        """Figure of the slice at the stage positions and as registered, with the cells.

        Higher-magnification images are outlined, and cell positions are marked (annotated and
        inferred positions differently).

        Parameters
        ----------
        filename : pathlib.Path
            PNG file to write.
        title : str
            Figure title.
        """
        canvas = self.overview_canvas()
        ranges = self.channel_ranges()
        low = min(t.mag for t in self.used)
        before = self.render(canvas, max_mag=low * 2.6, placed=False, ranges=ranges, cells=False)
        after = self.render(canvas, max_mag=low * 2.6, placed=True, ranges=ranges, cells=False)
        fig, axes = mpl.subplots(1, 2, figsize=(16, 8.6), layout="constrained")
        ext = (canvas.bounds[0], canvas.bounds[1], canvas.bounds[2], canvas.bounds[3])
        for ax, img, label in ((axes[0], before, "stage positions"), (axes[1], after, "registered")):
            ax.imshow(img, extent=ext, origin="upper", interpolation="antialiased")
            ax.set_title(label)
            ax.set_xlabel("x (um)")
            ax.set_aspect("equal")
        axes[0].set_ylabel("y (um)")
        ax = axes[1]
        for t in self.used:
            if t.mag > low * 2.6:
                b = t.bounds()
                ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fill=False, ec="cyan", lw=0.6, alpha=0.7))
        for cell, c in self.cells.items():
            if c["x"] is None:
                continue
            st = cell_style(c)
            ax.plot(c["x"], c["y"], st["marker"], ms=9, mfc="none", mec=st["color"], mew=1.5)
            ax.annotate(cell.replace("cell_", "C"), (c["x"], c["y"]), xytext=(6, 6), textcoords="offset points", color=st["color"], fontsize=9)
        ax.text(
            0.01, 0.01, "Magenta circles: cell positions placed by hand in the Mosaic Editor (.mosaic) file.\n"
            "Blue squares: cell positions inferred here (centre of the cell's highest-magnification image).\n"
            "C followed by a number labels each cell.   Cyan boxes: higher-magnification images.",
            transform=ax.transAxes, color="white", fontsize=7, va="bottom",
            bbox={"facecolor": "black", "alpha": 0.5, "edgecolor": "none"},
        )
        legend = sorted({(t.channel, t.group) for t in self.used if t.mag <= low * 2.6 and t.cell is None})
        text = "   ".join(f"{ch}: {g}" for ch, g in legend)
        fig.suptitle(f"{title}\n{text}", fontsize=10)
        self.stamp(fig)
        fig.savefig(filename, dpi=110)
        mpl.close(fig)

    def save_qc(self, filename: Path, title: str, max_groups: int = 6):
        """For each group that was moved: its edges (magenta) on the base's edges (green), at the
        stage positions and as registered. White where they coincide.

        An earlier figure of the same name is always removed first.

        Parameters
        ----------
        filename : pathlib.Path
            PNG file to write.
        title : str
            Figure title.
        max_groups : int, optional
            Largest number of groups shown (the largest groups; default 6).

        Returns
        -------
        bool
            True if a figure was written (False if no group was moved).
        """
        base = self.base
        low = min(t.mag for t in self.used)
        moved = {}
        for t in self.used:
            if t not in base and t.cell is None and t.mag <= low * 2.6:
                moved.setdefault((t.group, t.channel, t.mag), []).append(t)
        keys = sorted(moved, key=lambda k: -len(moved[k]))[:max_groups]
        if filename.exists():
            filename.unlink()  # never leave a figure from an earlier run
        if not keys:
            return False
        fig, axes = mpl.subplots(len(keys), 2, figsize=(12, 5.2 * len(keys)), layout="constrained", squeeze=False)
        for row, key in zip(axes, keys):
            grp = moved[key]
            b = expand(group_bounds(grp), 100.0)
            canvas = Canvas(*b, max(self.base_res, (b[1] - b[0]) / 800.0, grp[0].pixel * 2.0))
            bimg, bw = composite(base, canvas, "grad")
            for ax, placed in zip(row, (False, True)):
                saved = [(t.dx, t.dy) for t in grp]
                if not placed:  # stage position, relative to the (stitched) base
                    bdx = np.median([t.dx for t in base])
                    bdy = np.median([t.dy for t in base])
                    for t in grp:
                        t.dx, t.dy = bdx, bdy
                gimg, gw = composite(grp, canvas, "grad")
                for t, (dx, dy) in zip(grp, saved):
                    t.dx, t.dy = dx, dy
                rgb = np.zeros(canvas.shape + (3,), np.float32)
                rgb[..., 1] = np.clip(bimg / 4, 0, 1) * (bw > 0)
                rgb[..., 0] = rgb[..., 2] = np.clip(gimg / 4, 0, 1) * (gw > 0)
                ax.imshow(rgb, extent=canvas.bounds)
                ax.set_title(
                    f"{key[0]} {key[1]} {key[2]:g}x: "
                    + (
                        "stage positions"
                        if not placed
                        else f"registered ({sorted({t.placement for t in grp})}, median "
                        f"{np.median([t.dx for t in grp]):+.0f}, {np.median([t.dy for t in grp]):+.0f} um)"
                    ),
                    fontsize=9,
                )
        fig.suptitle(f"{title}: edges of each group (magenta) on the base (green), white where they are aligned", fontsize=10)
        self.stamp(fig)
        fig.savefig(filename, dpi=80)
        mpl.close(fig)
        return True

    def save_cell(self, cell: str, filename: Path, title: str, size: float = 1500.0, max_images: int = 12):
        """The slice around the cell (slice-level images, with the outlines of the
        higher-magnification images), and the cell's own highest-magnification images

        The cell's images at a higher magnification than the lowest (slice) magnification (all of
        its images if it has none) are shown at the right in file-number order (images, then videos), in 1 column
        (1-2 images), 2 columns (3-4) or 3 columns (more). Video projections are shown on a
        square-root scale, so that dim processes are visible next to a bright soma. Where the
        cell's position is outside an image, the panel says how far away and in which direction.

        Parameters
        ----------
        cell : str
            Cell name (a key of ``cells``).
        filename : pathlib.Path
            PNG file to write.
        title : str
            Figure title.
        size : float, optional
            Width of the region around the cell, um (default 1500).
        max_images : int, optional
            Largest number of the cell's images shown (default 12); the number left out is noted.

        Returns
        -------
        bool
            True if a figure was written (False if the cell has no position).
        """
        c = self.cells[cell]
        if c["x"] is None:
            return False
        ranges = self.channel_ranges()
        b = (c["x"] - size / 2, c["x"] + size / 2, c["y"] - size / 2, c["y"] + size / 2)
        low = min(t.mag for t in self.used)
        canvas = Canvas(*b, max(self.base_res / 2.0, size / 1200.0))
        img = self.render(canvas, max_mag=low * 2.6, ranges=ranges, cells=False)
        own = [t for t in self.used if t.cell == cell]
        own = [t for t in own if t.mag > low] or own
        own = sorted(own, key=lambda t: (t.group, t.kind != "image", file_number(t.path), t.path.name))
        shown = own[:max_images]
        n = len(shown)
        ncols = 1 if n <= 2 else (2 if n <= 4 else 3)
        nrows = max(int(np.ceil(n / ncols)), 1)
        side = max(7.0, 3.9 * nrows)  # the slice panel is square, as tall as the image panels
        fig = mpl.figure(figsize=(side + 5.0 * ncols * (n > 0), side + 1.6), layout="constrained")
        gs = fig.add_gridspec(nrows, 1 + ncols * (n > 0), width_ratios=[side] + [5.0] * ncols * (n > 0))
        ax = fig.add_subplot(gs[:, 0])
        ax.imshow(img, extent=canvas.bounds, origin="upper", interpolation="antialiased")
        for t in self.used:
            if (t.mag > low * 2.6 or t in shown) and overlaps(t.bounds(), b):
                tb = t.bounds()
                color = "yellow" if t in shown else ("cyan" if t.cell == cell else "0.6")
                ax.add_patch(Rectangle((tb[0], tb[2]), tb[1] - tb[0], tb[3] - tb[2], fill=False, ec=color, lw=0.8))
        for other, oc in self.cells.items():
            if oc["x"] is not None and overlaps((oc["x"],) * 2 + (oc["y"],) * 2, b):
                st = cell_style(oc)
                this = other == cell
                ax.plot(oc["x"], oc["y"], st["marker"], ms=12 if this else 7, mfc="none", mec=st["color"], mew=1.5 if this else 0.8)
                ax.annotate(
                    other.replace("cell_", "C"), (oc["x"], oc["y"]), xytext=(6, 6), textcoords="offset points",
                    color=st["color"], fontsize=10 if this else 8, fontweight="bold" if this else "normal",
                )
        ax.set_xlim(b[0], b[1])
        ax.set_ylim(b[2], b[3])
        how = "placed by hand" if cell_style(c) is annotated_style else "INFERRED"
        ax.set_title(f"{size:.0f} um around {cell} (C{cell[5:]}): position {how} ({c['source']})", fontsize=9)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        inferred = cell_style(c) is inferred_style
        mark = "x" if inferred else "+"
        for k, t in enumerate(shown):
            ax = fig.add_subplot(gs[k // ncols, 1 + k % ncols])
            tb = t.bounds()
            if t.channel == "BF":
                lo, hi = np.percentile(t.data, [0.5, 99.8])
                g = np.clip((t.data - lo) / max(hi - lo, 1e-6), 0, 1)
                rgb = np.repeat(g[..., None], 3, axis=2)
            else:
                lo, hi = np.percentile(t.data, [5, 99.9])
                g = np.clip((t.data - lo) / max(hi - lo, 1e-6), 0, 1)
                rgb = g[..., None] * np.array(channel_colors[t.channel])[None, None, :]
            if t.kind == "video":  # dye fills: a bright soma and dim processes; compress the range
                rgb = rgb**0.5
            ax.imshow(rgb, extent=tb, origin="upper")
            if overlaps((c["x"],) * 2 + (c["y"],) * 2, tb):
                ax.plot(c["x"], c["y"], mark, ms=14, color=cell_style(c)["color"], mew=1.2)
            else:  # say where the cell is, rather than leave the mark out
                cx, cy = (tb[0] + tb[1]) / 2, (tb[2] + tb[3]) / 2
                dx, dy = c["x"] - cx, c["y"] - cy
                ax.text(
                    0.5, 0.04, f"{cell} position is outside this image: {np.hypot(dx, dy):.0f} um from its centre, "
                    f"toward {compass(dx, dy)}",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=cell_style(c)["color"],
                    bbox={"facecolor": "black", "alpha": 0.6, "edgecolor": "none"},
                )
            ax.set_xlim(tb[0], tb[1])
            ax.set_ylim(tb[2], tb[3])
            shown_as = f"\n{t.projection['method']} projection of {t.frames:d} frames, square-root scale" if t.projection else ""
            ax.set_title(f"{t.name} ({t.mag:g}x, {t.channel}, {t.placement}){shown_as}", fontsize=8)
            ax.tick_params(labelsize=7)
        legend = (
            "Magenta circle (left) and magenta + (right): cell position placed by hand in the Mosaic Editor (.mosaic) file.\n"
            "Blue square (left) and blue x (right): cell position inferred here (centre of the cell's highest-magnification "
            "image), used only when there is no hand-placed position.\n"
            "C followed by a number labels each cell (bold: this cell).   Boxes on the left: yellow = the images at the right, "
            "cyan = other images of this cell, grey = images of other cells."
            + (f"\n{len(own) - n:d} more images of this cell are not shown." if len(own) > n else "")
        )
        fig.suptitle(f"{title}\n{legend}", fontsize=8)
        self.stamp(fig)
        fig.savefig(filename, dpi=100)
        mpl.close(fig)
        return True

    def stamp(self, fig):
        """Add the provenance stamps to a figure.

        Top left: the repository, its version and git commit. Bottom right: when the slice was
        processed.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure.
        """
        fig.text(0.003, 0.997, provenance(), ha="left", va="top", fontsize=8, color="0.3")
        fig.text(0.997, 0.003, f"processed {self.processed}", ha="right", va="bottom", fontsize=8, color="0.3")

    def exceptions(self, min_mag: float = 10.0) -> dict:
        """The exceptions of this slice, to be checked by hand (see :mod:`src.mosaic_exceptions`).

        Parameters
        ----------
        min_mag : float, optional
            Smallest magnification listed in "image kind" (default 10: the appearance rule is
            meant for 10x, 40x and 63x).

        Returns
        -------
        dict
            Sheet name -> rows: "image kind" (appearance ambiguous, or set in the decision
            table), "image kind (automatic)" (channel set from a confident verdict of the
            classifier, overriding the folder name), "refused matches" (group moves of more than 250 um) and "annotation
            mismatches" (hand-placed cell positions outside all of the cell's own images).
        """
        out = {"image kind": [], "image kind (automatic)": [], "refused matches": [], "annotation mismatches": []}
        for t in self.tiles:
            a = t.appearance
            if t.excluded == "excluded in the decision table":  # keep decided rows in the table
                out["image kind"].append(
                    {"slice": self.title, "image": t.name, "mag": t.mag, "kind": t.kind, "now treated as": "excluded (decision table)", "slice_dir": str(self.slice_dir)}
                )
                continue
            if a is None or not np.isfinite(t.mag) or t.mag < min_mag:
                continue
            now_says = "BF" if t.channel == "BF" else "fluorescence"
            automatic = t.channel_from == "appearance"
            if not automatic and t.channel_from != "decision" and a["rule"] == now_says:
                continue
            sheet = "image kind (automatic)" if automatic else "image kind"
            out[sheet].append(
                {
                    "slice": self.title,
                    "image": t.name,
                    "mag": t.mag,
                    "kind": t.kind,
                    "now treated as": ("mixed" if t.mixed else t.channel) + (" (decision table)" if t.channel_from == "decision" else f" ({t.channel_from})"),
                    "appearance rule": a["rule"],
                    "p(fluorescence)": a["p_fluorescence"],
                    "structured_skew": a["structured_skew"],
                    "dark_excursion": a["dark_excursion"],
                    "tail_ratio": a["tail_ratio"],
                    "raw_tail_ratio": a["raw_tail_ratio"],
                    "texture_spread": a["texture_spread"],
                    "bright_area": a["bright_area"],
                    "n_bright": a["n_bright"],
                    "excluded": t.excluded,
                    "slice_dir": str(self.slice_dir),
                    "_thumbnail": thumbnail_jpeg(t.data, width=200, gamma=0.5) if t.data is not None else None,
                    "_aspect": t.data.shape[0] / t.data.shape[1] if t.data is not None else 0.75,
                }
            )
            if out[sheet][-1]["_thumbnail"] is None:
                del out[sheet][-1]["_thumbnail"]
        for g in self.report.get("groups", []):
            if "refused_move_um" not in g:
                continue
            dx, dy = g["refused_move_um"]
            out["refused matches"].append(
                {
                    "slice": self.title,
                    "group": g["group"],
                    "channel": g["channel"],
                    "mag": g["mag"],
                    "run": g.get("run", 0),
                    "images": g["n"],
                    "move_dx_um": dx,
                    "move_dy_um": dy,
                    "move_um": round(float(np.hypot(dx, dy)), 1),
                    "ncc": g.get("refused_ncc"),
                    "applied": g.get("decision", "stage"),
                    "slice_dir": str(self.slice_dir),
                }
            )
        low = min((t.mag for t in self.used), default=float("nan"))
        for cell, c in self.cells.items():
            if c["x"] is not None and c["source"].startswith("annotated"):
                ax, ay, source = c["x"], c["y"], c["source"]
            elif "annotated" in c:  # overridden in the decision table: still listed
                ax, ay, source = c["annotated"][0], c["annotated"][1], f"annotated in {c['annotated'][2]}"
            else:
                continue
            own = [t for t in self.used if t.cell == cell and t.mag > low]
            if not own or any(overlaps((ax,) * 2 + (ay,) * 2, t.bounds()) for t in own):
                continue
            c = {"x": ax, "y": ay, "source": source}
            near = min(own, key=lambda t: np.hypot(c["x"] - t.center()[0], c["y"] - t.center()[1]))
            cx, cy = near.center()
            top = max(own, key=lambda t: (t.mag, t.timestamp))
            out["annotation mismatches"].append(
                {
                    "slice": self.title,
                    "cell": cell,
                    "annotated_x_um": round(c["x"], 1),
                    "annotated_y_um": round(c["y"], 1),
                    "source": c["source"],
                    "nearest own image": near.name,
                    "distance_um": round(float(np.hypot(c["x"] - cx, c["y"] - cy)), 1),
                    "annotation is toward": compass(c["x"] - cx, c["y"] - cy),
                    "image placement": near.placement,
                    "inferred_x_um": round(top.center()[0], 1),
                    "inferred_y_um": round(top.center()[1], 1),
                    "slice_dir": str(self.slice_dir),
                }
            )
        return out

    def tile_table(self) -> list:
        """One entry per tile for the JSON report, including the metadata carried forward.

        Returns
        -------
        list of dict
            For each tile: name, kind, channel, objective, magnification, pixel size, stage position,
            automatic and manual corrections, placement, exclusion, size, notes, projection (videos),
            ``corrected`` (the placement as an acq4 transform; see :meth:`Tile.corrected_transform` and
            :meth:`Tile.mosaic_user_transform`) and ``index`` (the .index entry as recorded).
        """
        return [
            {
                "name": t.name,
                "kind": t.kind,
                "channel": t.channel,
                "channel_from": t.channel_from,
                "objective": t.objective,
                "mag": t.mag,
                "pixel_um": round(t.pixel, 4),
                "stage_xy_um": [round(t.xmin, 1), round(t.ymax, 1)],
                "correction_um": [round(t.dx, 1), round(t.dy, 1)],
                "manual_correction_um": [round(t.user_offset[0], 1), round(t.user_offset[1], 1)],
                "placement": t.placement,
                "excluded": t.excluded,
                "MB": round(t.nbytes / 1e6, 1),
                "notes": t.notes,
                "projection": t.projection,
                "appearance": t.appearance,
                # carried forward for atlas alignment: the .index entry as recorded, and the
                # placement made here in the same form (see Tile.corrected_transform)
                "corrected": None
                if t.excluded
                else {
                    "objective_index": t.objective,
                    "mag_index": t.index_mag,
                    "mag_used": t.mag,
                    "pixel_m": [t.pixel * 1e-6, t.pixel * 1e-6],
                    "bounds_um": [round(v, 2) for v in t.bounds()],
                    "transform": t.corrected_transform(),
                    "mosaic_userTransform": t.mosaic_user_transform(),
                },
                "index": json_safe(t.index_info),
            }
            for t in self.tiles
        ]

    def save_mosaic_file(self, filename: Path, levels: Optional[dict] = None):
        """Write the placements as an acq4 Mosaic Editor file (.mosaic).

        Each image (and video) is an ImageCanvasItem whose userTransform places it as registered here
        (see :meth:`Tile.mosaic_user_transform`). Cells are CellCanvasItems; inferred positions are
        named "cell_nnn (inferred)". Landmark (MarkersCanvasItem) and ruler items from the Mosaic
        Editor files in the slice directory are copied unchanged. ``rootPath`` is the slice
        directory and file names are relative to it. The acq4 Mosaic Editor resolves file names
        relative to the directory holding the .mosaic file, so to open it there, copy it into the
        slice directory.

        Parameters
        ----------
        filename : pathlib.Path
            The file to write.
        levels : dict, optional
            Display range for each channel (default: each image's 0.5-99.8 percentiles).
        """
        items = []
        low = min(t.mag for t in self.used)
        for z, t in enumerate(sorted(self.used, key=lambda t: (t.mag, t.timestamp))):
            rel = t.path.relative_to(self.slice_dir).as_posix()
            lo, hi = np.percentile(t.data[::4, ::4], [0.5, 99.8]) if levels is None else levels[t.channel]
            items.append(
                {
                    "type": "ImageCanvasItem",
                    "name": f"{self.slice_dir.name}/{t.path.name}" if t.group == "." else rel,
                    "visible": t.mag <= low * 2.6 and t.cell is None,
                    "alpha": 1.0 if t.channel == "BF" else 0.6,
                    "userTransform": t.mosaic_user_transform(),
                    "z": z,
                    "scalable": False,
                    "rotatable": True,
                    "movable": True,
                    "filename": rel,
                    "imagestate": {
                        "gradient": {"mode": "rgb", "ticks": [[0.0, [0, 0, 0, 255]], [1.0, [255, 255, 255, 255]]], "ticksVisible": True},
                        "levels": [float(lo), float(hi)],
                        "mode": "mono",
                    },
                    "filter": _mosaic_filter_state,
                    "composition": "SourceOver" if t.channel == "BF" else "Plus",
                }
            )
        z = len(items)
        for cell, c in self.cells.items():
            if c["x"] is None:
                continue
            inferred = cell_style(c) is inferred_style
            items.append(
                {
                    "type": "CellCanvasItem",
                    "name": f"{cell} (inferred)" if inferred else cell,
                    "visible": True,
                    "alpha": 1.0,
                    "userTransform": {"pos": [c["x"] * 1e-6, c["y"] * 1e-6], "scale": [1.0, 1.0], "angle": 0.0},
                    "z": z,
                    "scalable": False,
                    "rotatable": False,
                    "movable": True,
                    "filename": None,
                }
            )
            z += 1
        # landmarks and rulers from the Mosaic Editor files in the slice directory, unchanged
        for mf in sorted(self.slice_dir.glob("*.mosaic")):
            try:
                kept = [i for i in json.load(open(mf))["items"] if i.get("type") in ("MarkersCanvasItem", "RulerCanvasItem")]
            except Exception:
                continue
            items.extend(kept)
        b = expand(group_bounds([t for t in self.used if t.mag <= low * 2.6 and t.cell is None] or self.used), 50.0)
        view = dict(_mosaic_view_state)
        view["targetRange"] = view["viewRange"] = [[b[0] * 1e-6, b[1] * 1e-6], [b[2] * 1e-6, b[3] * 1e-6]]
        state = {
            "contents": "MosaicEditor_save",
            "version": [2, 0],
            "rootPath": str(self.slice_dir),
            "items": items,
            "view": view,
        }
        with open(filename, "w") as f:
            json.dump(state, f, indent=1)

    def save_projections(self, directory: Path) -> list:
        """Write the projection of each video as a float32 TIFF; returns the files written.

        The TIFFs are oriented with row 0 at the top (largest y) and column 0 at the left (smallest
        x); their metadata give the source file, the projection, the pixel size and the placement
        (bounds in um).

        Parameters
        ----------
        directory : pathlib.Path
            Output directory; the slice's folder structure is kept below it.

        Returns
        -------
        list of pathlib.Path
            The files written.
        """
        written = []
        for t in self.used:
            if t.kind != "video" or t.data is None:
                continue
            out = directory / t.group / (t.path.stem + "_projection.tif")
            out.parent.mkdir(parents=True, exist_ok=True)
            b = t.bounds()
            tifffile.imwrite(
                out,
                t.data.astype(np.float32),
                metadata={
                    "source": str(t.path),
                    "projection": t.projection,
                    "orientation": "row 0 is the top (largest y), column 0 the left (smallest x)",
                    "pixel_um": t.pixel,
                    "bounds_um": [round(v, 2) for v in b],
                },
            )
            written.append(out)
        return written


def slice_title(slice_dir: Path) -> str:
    """A short name for a slice: its last three path components.

    Parameters
    ----------
    slice_dir : pathlib.Path
        The slice directory.

    Returns
    -------
    str
        e.g. ``"Rig2/2023.10.30_000/slice_003"``.
    """
    parts = slice_dir.parts
    return "/".join(parts[-3:]) if len(parts) >= 3 else str(slice_dir)


def run_slice(
    slice_dir: Path,
    output: Path,
    cells: Optional[list] = None,
    max_video_mb: float = 500.0,
    video_projection: str = "focus",
    cell_region_um: float = 1500.0,
    decisions: Optional[Decisions] = None,
    auto_kind: bool = True,
):
    """Build the mosaic of one slice and write its figures, files and report.

    Writes ``<name>~mosaic.png``, ``~qc.png``, ``~<cell>.png``, ``~aligned.mosaic``,
    ``~projections/`` and ``~report.json`` to the output directory, where name is the slice title
    with "/" replaced by "~".

    Parameters
    ----------
    slice_dir : pathlib.Path
        The slice directory.
    output : pathlib.Path
        Output directory (created if needed).
    cells : list of str, optional
        Cells to draw (default: all).
    max_video_mb : float, optional
        Videos larger than this (MB) are excluded (default 500).
    video_projection : {"focus", "mip", "best"}, optional
        How videos are projected to one image (default "focus").
    cell_region_um : float, optional
        Width of the region around each cell in the cell figures, um (default 1500).
    decisions : Decisions, optional
        Decisions from an exceptions workbook.
    auto_kind : bool, optional
        Let confident brightfield/fluorescence verdicts override the folder name (default True).

    Returns
    -------
    SliceMosaic
        The built mosaic.
    """
    title = slice_title(slice_dir)
    name = title.replace("/", "~")
    print(f"{title}:")
    mosaic = SliceMosaic(
        slice_dir, max_video_mb=max_video_mb, video_projection=video_projection, decisions=decisions, auto_kind=auto_kind
    )
    excluded = [t for t in mosaic.tiles if t.excluded]
    print(f"   {len(mosaic.tiles):d} images/videos, {len(excluded):d} left out")
    mosaic.build()
    output.mkdir(parents=True, exist_ok=True)
    report = dict(mosaic.report)
    report["tiles"] = mosaic.tile_table()
    if "error" in report:
        print(f"   {report['error']}")
    else:
        mosaic.save_overview(output / f"{name}~mosaic.png", title)
        mosaic.save_mosaic_file(output / f"{name}~aligned.mosaic", levels=mosaic.channel_ranges())
        mosaic.save_projections(output / f"{name}~projections")
        mosaic.save_qc(output / f"{name}~qc.png", title)
        for cell in cells if cells is not None else list(mosaic.cells):
            if cell in mosaic.cells and mosaic.save_cell(
                cell, output / f"{name}~{cell}.png", f"{title}/{cell}", size=cell_region_um
            ):
                pass
        for row in report.get("manual_vs_automatic", []):
            print(
                f"   manual vs automatic, {row['group']}: manual {row['manual_um']}, "
                f"automatic {row['automatic_um']} um (difference {row['difference_um']} um)"
            )
    with open(output / f"{name}~report.json", "w") as f:
        json.dump(report, f, indent=1, default=str)
    _resample_cache.clear()
    return mosaic


def random_cells(root: Path, n: int, seed: int) -> list:
    """n cell directories (with recordings) chosen at random under root

    Parameters
    ----------
    root : pathlib.Path
        acq4 data root (searched for day directories ``yyyy.mm.dd_nnn``).
    n : int
        Number of cells.
    seed : int
        Random seed.

    Returns
    -------
    list of pathlib.Path
        The cell directories.
    """
    cells = []
    for day in sorted(root.rglob("20[0-9][0-9].[0-9][0-9].[0-9][0-9]_[0-9][0-9][0-9]")):
        for cell in sorted(day.glob("slice_*/cell_*")):
            if cell.is_dir():
                cells.append(cell)
    return random.Random(seed).sample(cells, min(n, len(cells)))


def main():
    """Command line entry point (``acq4_mosaic``); see the module docstring for usage."""
    parser = argparse.ArgumentParser(description="Experimental: build slice mosaics from acq4 images")
    parser.add_argument("slices", nargs="*", type=Path, help="slice directories")
    parser.add_argument("-o", "--output", type=Path, default=Path("mosaic_test"), help="output directory")
    parser.add_argument("--root", type=Path, action="append", default=[], help="acq4 data root(s) for --random")
    parser.add_argument("--random", type=int, default=0, help="number of cells to choose at random")
    parser.add_argument("--seed", type=int, default=None, help="random seed")
    parser.add_argument("--max-video-mb", type=float, default=500.0, help="leave out larger videos")
    parser.add_argument(
        "--video-projection",
        choices=["focus", "mip", "best"],
        default="focus",
        help="how a video (z-stack) becomes one image: multi-scale focus stack (default), "
        "maximum projection, or the sharpest frame",
    )
    parser.add_argument(
        "--cell-region-um",
        type=float,
        default=None,
        help="width of the region around each cell in the cell figures, um (default: the dataset's "
        "mosaic_cell_region_um if --dataset is given, otherwise 1500)",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="dataset entry in the configuration file to take settings from (mosaic_cell_region_um)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "config" / "dandi_configuration.cfg",
        help="configuration file for --dataset",
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=None,
        help="exceptions workbook (from an earlier run) with decisions filled in",
    )
    parser.add_argument(
        "--no-auto-kind",
        action="store_true",
        help="do not let the brightfield/fluorescence classifier override folder names",
    )
    args = parser.parse_args()
    decisions = Decisions(args.decisions) if args.decisions is not None else Decisions()
    if args.decisions is not None:
        print(f"{len(decisions):d} decisions read from {args.decisions}")
    cell_region_um = 1500.0
    if args.dataset is not None:
        config = configfile.readConfigFile(str(args.config))
        if args.dataset not in config:
            raise SystemExit(f"dataset {args.dataset} is not in {args.config}")
        cell_region_um = float(config[args.dataset].get("mosaic_cell_region_um", cell_region_um))
    if args.cell_region_um is not None:
        cell_region_um = args.cell_region_um
    jobs = [(Path(s), None) for s in args.slices]
    if args.random > 0:
        seed = args.seed if args.seed is not None else random.randrange(10000)
        print(f"random seed {seed:d}")
        chosen = []
        for root in args.root:
            chosen += random_cells(root, args.random, seed)
        chosen = random.Random(seed).sample(chosen, min(args.random, len(chosen)))
        for cell in chosen:
            jobs.append((cell.parent, [cell.name]))
    exceptions = {}
    for slice_dir, cells in jobs:
        mosaic = run_slice(
            slice_dir,
            args.output,
            cells=cells,
            max_video_mb=args.max_video_mb,
            video_projection=args.video_projection,
            cell_region_um=cell_region_um,
            decisions=decisions,
            auto_kind=not args.no_auto_kind,
        )
        for sheet, rows in mosaic.exceptions().items():
            exceptions.setdefault(sheet, []).extend(rows)
        mosaic.tiles = []  # free the images before the next slice
    workbook = args.output / f"exceptions_{time.strftime('%Y%m%d-%H%M%S')}.xlsx"
    write_workbook(workbook, exceptions, decisions)
    print(
        f"exceptions: {len(exceptions.get('image kind', [])):d} image kind, {len(exceptions.get('image kind (automatic)', [])):d} "
        f"set automatically, {len(exceptions.get('refused matches', [])):d} refused "
        f"matches, {len(exceptions.get('annotation mismatches', [])):d} annotation mismatches -> {workbook}"
    )


if __name__ == "__main__":
    main()
