"""Projections of acq4 video stacks (z-series taken while focusing by hand) to one image.

A video stack in acq4 (``video_NNN.ma``) is usually a series of frames taken while the
experimenter turned the focus knob. In most stacks the focus position is not recorded (the
per-frame ``translation`` in the MetaArray is constant); where it is recorded, the steps are
irregular and can reverse. So which frames are in focus, and where, has to be found from the
images themselves.

A plain maximum-intensity projection (MIP) has two problems with these stacks, most of all for
dye-filled cells:

1. Out-of-focus light. A bright soma, out of focus, spreads into a wide halo ("flare") that is
   brighter than the in-focus fine processes under it, so the MIP shows the halo there, not the
   processes.
2. Noise. The maximum of n noisy frames is biased upward, more so as n grows.

:func:`project_stack` offers three methods:

``"focus"`` (default)
    Multi-scale focus stacking (extended depth of field). Each frame, lightly
    smoothed, is split into a Laplacian pyramid: band-pass images from fine detail (processes)
    to the size of a soma, plus a smooth remainder. In each band, at each pixel, the value is
    taken from the frame with the most local energy in that band, i.e. the frame in which the
    structure of that size is most in focus. Out-of-focus light is smooth, so it carries little
    band energy and is not picked. The smooth remainder (the coarsest band, larger than a soma)
    is where the flare is; it is taken as the minimum over the frames, because at any place the
    flare from a structure is least in the frame where that structure is in focus. Where no
    frame stands out in a band (its best energy is not well above the mean energy over the
    frames: nothing of that size is in focus there, only noise), the band is averaged over the
    frames instead of selected, so the noise is averaged down rather than picked at its largest
    (the "match measure" of Burt and Kolczynski, 1993). No pixel is allowed below its value in
    the darkest frame. This is the pyramid fusion of Burt and Adelson (1985), as used for
    extended depth of field in microscopy (e.g. Forster et al., 2004).
``"mip"``
    Maximum-intensity projection of the lightly smoothed frames (what the acq4 Mosaic
    Editor does), for comparison.
``"best"``
    The single sharpest frame (overall).

None of the methods uses the focal step, so irregular or reversed steps, and several frames at the
same focus, need no special treatment. Frames much darker than the rest (shutter closed, lamp
off) are left out.

References
----------
P. J. Burt and R. J. Kolczynski, "Enhanced image capture through fusion," Proc. Fourth
International Conference on Computer Vision (1993), 173-182, doi:10.1109/ICCV.1993.378222.

P. J. Burt and E. H. Adelson, "Merging images through pattern decomposition," Proc. SPIE 575,
Applications of Digital Image Processing VIII (1985), doi:10.1117/12.966501.

B. Forster, D. Van De Ville, J. Berent, D. Sage and M. Unser, "Complex wavelets for extended
depth-of-field: a new method for the fusion of multichannel microscopy images," Microscopy
Research and Technique 65 (2004), 33-42, doi:10.1002/jemt.20092.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import scipy.ndimage as ndi


@dataclass
class ProjectionInfo:
    """What :func:`project_stack` did, for reports.

    Attributes
    ----------
    method : str
        "focus", "mip" or "best".
    n_frames : int
        Frames in the stack.
    used : list of int
        Indices of the frames that were used.
    dark : list of int
        Indices of the frames left out as dark.
    frame_sharpness : list of float
        Overall sharpness of each used frame, relative to the sharpest.
    best_frame : int
        Index of the sharpest frame.
    levels : int
        Number of band-pass levels in the pyramid ("focus" only).
    focus_index : numpy.ndarray or None
        For each pixel, the index of the frame its finest detail was taken from, or -1 where no
        frame stood out and the detail was averaged ("focus" only); shows which parts of the
        image came from which focal plane.
    """

    method: str = "focus"
    n_frames: int = 0
    used: list = field(default_factory=list)
    dark: list = field(default_factory=list)
    frame_sharpness: list = field(default_factory=list)
    best_frame: int = -1
    levels: int = 0
    focus_index: Optional[np.ndarray] = field(default=None, repr=False)

    def summary(self) -> str:
        """One line describing the projection, for the notes of a tile.

        Returns
        -------
        str
            The method, the number of frames used and left out, and the sharpest frame.
        """
        what = {"focus": "multi-scale focus stack", "mip": "maximum projection", "best": "sharpest frame"}[self.method]
        return (
            f"{what} of {len(self.used):d}/{self.n_frames:d} frames"
            + (f" ({len(self.dark):d} dark frames left out)" if self.dark else "")
            + f"; sharpest frame {self.best_frame:d}"
            + (f"; {self.levels:d} pyramid levels" if self.method == "focus" else "")
        )


def _up(a: np.ndarray, shape: tuple) -> np.ndarray:
    """Expand a pyramid level to a given shape (bilinear).

    Parameters
    ----------
    a : numpy.ndarray
        2D image.
    shape : tuple of int
        Target (rows, cols).

    Returns
    -------
    numpy.ndarray
        ``a`` resampled to ``shape``.
    """
    return ndi.zoom(a, (shape[0] / a.shape[0], shape[1] / a.shape[1]), order=1)


def laplacian_pyramid(a: np.ndarray, levels: int) -> tuple:
    """Split an image into band-pass levels and a smooth remainder.

    Level 0 holds the finest detail, at full resolution; each further level is half the size and
    holds detail twice as large. ``reconstruct(*laplacian_pyramid(a, n))`` gives back ``a``
    exactly (to rounding).

    Parameters
    ----------
    a : numpy.ndarray
        2D image.
    levels : int
        Number of band-pass levels.

    Returns
    -------
    bands : list of numpy.ndarray
        The band-pass images, finest first.
    residual : numpy.ndarray
        The smooth remainder (the coarsest level).
    """
    bands = []
    g = a
    for _ in range(levels):
        down = ndi.gaussian_filter(g, 1.0)[::2, ::2]
        bands.append(g - _up(down, g.shape))
        g = down
    return bands, g


def reconstruct(bands: list, residual: np.ndarray) -> np.ndarray:
    """Inverse of :func:`laplacian_pyramid`.

    Parameters
    ----------
    bands : list of numpy.ndarray
        Band-pass images, finest first.
    residual : numpy.ndarray
        The smooth remainder.

    Returns
    -------
    numpy.ndarray
        The image.
    """
    a = residual
    for b in reversed(bands):
        a = _up(a, b.shape) + b
    return a


def frame_sharpness(frame: np.ndarray, detail_px: float = 1.0) -> float:
    """Overall sharpness of a frame: variance of the Laplacian of Gaussian, relative to brightness.

    Parameters
    ----------
    frame : numpy.ndarray
        2D image.
    detail_px : float, optional
        Scale (sigma, pixels) of the Laplacian of Gaussian (default 1).

    Returns
    -------
    float
        The sharpness (larger is sharper). Dividing by the squared mean makes frames of different
        brightness comparable.
    """
    return float(ndi.gaussian_laplace(frame, detail_px).var() / max(float(frame.mean()), 1.0) ** 2)


def project_stack(
    get_frame: Callable[[int], np.ndarray],
    n_frames: int,
    pixel_um: float,
    method: str = "focus",
    smooth_px: float = 0.7,
    largest_um: float = 20.0,
    salience: tuple = (2.0, 4.0),
    dark_fraction: float = 0.5,
) -> tuple:
    """Project a z-stack to one image (see the module docstring for the methods).

    The frames are read one at a time and never all held in memory: the focus stack keeps, for
    each pyramid level, the best energy so far and the value that goes with it.

    Parameters
    ----------
    get_frame : callable
        ``get_frame(i)`` returns frame i as a 2D array (rows, cols).
    n_frames : int
        Number of frames.
    pixel_um : float
        Pixel size, um. It sets the number of pyramid levels.
    method : {"focus", "mip", "best"}, optional
        Projection method (default "focus").
    smooth_px : float, optional
        Sigma (pixels) of the light Gaussian smoothing of each frame (default 0.7).
    largest_um : float, optional
        Size (um) of the largest structure that is focus-selected as a whole (default 20 um,
        about a soma). The pyramid has enough levels that its coarsest band is at least this
        large; anything larger (the flare) goes to the remainder, which is the minimum over the
        frames.
    salience : tuple of float, optional
        (low, high) limits on the ratio of the best band energy at a pixel to the mean band
        energy over the frames (default (2, 4)). Below low, the band is averaged over the frames
        (nothing in focus; for pure noise the ratio is about 1.5-2 with tens of frames); above
        high, it is taken from the best frame; in between, a blend.
    dark_fraction : float, optional
        Frames whose median intensity is less than this fraction of the median over the frames
        are left out as dark (default 0.5).

    Returns
    -------
    image : numpy.ndarray
        The projection (float32, the shape of a frame), in the frames' intensity units.
    info : ProjectionInfo
        What was done.
    """
    if method not in ("focus", "mip", "best"):
        raise ValueError(f"unknown projection method {method!r}")
    info = ProjectionInfo(method=method, n_frames=n_frames)

    def frame(i):
        """Frame i, lightly smoothed.

        Parameters
        ----------
        i : int
            Frame index.

        Returns
        -------
        numpy.ndarray
            The frame (float32), smoothed with sigma ``smooth_px``.
        """
        return ndi.gaussian_filter(np.asarray(get_frame(i), dtype=np.float32), smooth_px)

    level = np.array([float(np.median(np.asarray(get_frame(i))[::4, ::4])) for i in range(n_frames)])
    info.used = [i for i in range(n_frames) if level[i] >= dark_fraction * np.median(level)]
    info.dark = [i for i in range(n_frames) if i not in info.used]
    levels = max(1, int(np.ceil(np.log2(max(largest_um / pixel_um, 2.0)))))
    info.levels = levels if method == "focus" else 0
    sharp = []
    best_img = best_e = best_v = lowest = mip = floor = None
    for i in info.used:
        f = frame(i)
        sharp.append(frame_sharpness(f))
        floor = f if floor is None else np.minimum(floor, f)
        if method == "mip":
            mip = f if mip is None else np.maximum(mip, f)
        elif method == "best":
            if sharp[-1] >= max(sharp):
                best_img = f
        else:
            bands, residual = laplacian_pyramid(f, levels)
            energy = [ndi.gaussian_filter(b * b, 1.0) for b in bands]
            if best_e is None:
                best_e, best_v, lowest = energy, bands, residual
                sum_e = [e.copy() for e in energy]
                sum_v = [b.copy() for b in bands]
                index = np.full(bands[0].shape, i, dtype=np.int32)
            else:
                for k in range(levels):
                    better = energy[k] > best_e[k]
                    best_e[k] = np.where(better, energy[k], best_e[k])
                    best_v[k] = np.where(better, bands[k], best_v[k])
                    sum_e[k] += energy[k]
                    sum_v[k] += bands[k]
                    if k == 0:
                        index = np.where(better, i, index)
                lowest = np.minimum(lowest, residual)
    rel = np.array(sharp) / max(max(sharp), 1e-12)
    info.frame_sharpness = [round(float(v), 3) for v in rel]
    info.best_frame = info.used[int(np.argmax(rel))]
    if method == "mip":
        image = mip
    elif method == "best":
        image = best_img
    else:
        # bands taken from different frames do not always add up: where a strong negative lobe
        # (the dark surround of an in-focus soma) is not matched by its own frame's remainder,
        # the sum dips below the background, so no pixel is allowed below its darkest frame
        n = len(info.used)
        lo, hi = salience
        fused = []
        for k in range(levels):
            ratio = best_e[k] / (sum_e[k] / n + 1e-12)
            w = np.clip((ratio - lo) / (hi - lo), 0.0, 1.0)
            fused.append(w * best_v[k] + (1.0 - w) * sum_v[k] / n)
            if k == 0:
                info.focus_index = np.where(w > 0.5, index, -1)
        image = np.maximum(reconstruct(fused, lowest), floor)
    return image.astype(np.float32), info
