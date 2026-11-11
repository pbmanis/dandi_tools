# Video stacks: projecting to one image

acq4 videos (`video_NNN.ma`) in these data are mostly z-stacks taken while turning the focus by
hand, often of dye-filled cells or fluorescent labels. Most do not record the focus position (the
per-frame `translation` in the MetaArray is constant). Where it is recorded, the steps are
irregular and can reverse. So which frames are in focus, and where, is found from the images.
The code is in `src/stack_projection.py` ({func}`src.stack_projection.project_stack`).

## Why not a maximum projection

A maximum-intensity projection (MIP, what the acq4 Mosaic Editor shows) has two problems with
these stacks:

1. **Flare.** A bright soma out of focus spreads into a wide halo that is brighter than the
   in-focus fine processes under it, so the MIP shows the halo there.
2. **Noise.** The maximum of many noisy frames is biased upward, more so as the number of frames
   grows.

## The focus stack (default)

Each frame is smoothed lightly (Gaussian, σ = 0.7 pixel) and split into a Laplacian pyramid:
band-pass images from fine detail (processes) up to the size of a soma (20 µm), plus a smooth
remainder. Then, pixel by pixel:

* **In each band**, the value is taken from the frame with the most local energy in that band,
  which is the frame in which structure of that size is most in focus. Out-of-focus light is
  smooth, carries little band energy, and is not picked.
* **Where no frame stands out** in a band (its best energy is less than 2–4 times the mean
  energy over the frames: nothing of that size is in focus there, only noise), the band is
  averaged over the frames instead, so the noise is averaged down rather than picked at its
  largest.
* **The smooth remainder** (larger than a soma), where the flare is, is the minimum over the
  frames: at any place, the flare from a structure is least in the frame where that structure
  is in focus.
* **No pixel is allowed below its darkest frame.** Bands taken from different frames do not
  always add up: the dark surround of an in-focus soma, selected from one frame and not matched
  by that frame's remainder, would otherwise dip below the background.

None of this uses the focal step, so irregular or reversed steps, and repeated frames at one
focus, need no special treatment. Frames much darker than the rest (median below half the
median frame: shutter closed) are left out. The frames are read one at a time, so large stacks
need little memory.

The method is the pyramid fusion of Burt and Adelson (1985) with the selection rule of Burt and
Kolczynski (1993), as used for extended depth of field in microscopy (Forster et al., 2004).

## Measured effect

On three stacks, comparing the focus stack with a MIP of the same smoothed frames (soma = the
brightest blob; halo = the median 15–35 µm from it, above the background more than 80 µm away;
noise = robust SD of a high-pass image far from the cell):

| Stack | Halo, % of soma: MIP → focus stack | Noise: MIP → focus stack |
|---|---|---|
| DCN_IC_inj 2025.03.06 s1 cell_000 video_000 (40x, 32 frames) | 15.0 → 7.9 | 20.5 → 19.5 |
| DCN_IC_inj 2025.02.06 s1 cell_001 video_000 (40x, 43 frames) | 29.2 → 15.4 | 26.6 → 31.0 |
| Thalamocortical Rig2 2022.12.15 cell_001 video_005 (40x, 60 frames, dim and noisy) | 41.3 → 38.5 | 1.8 → 2.8 |

The halo is about halved where there is a bright soma; noise is similar, and higher on the dim,
noisy stack. The soma itself is 25–30% dimmer than in the MIP, which includes the soma's own
flare.

## Other methods

`--video-projection mip` gives the MIP of the smoothed frames; `best` gives the single sharpest
frame (by the variance of the Laplacian relative to the brightness), close to the previous
behaviour of `acq4_mosaic` (which picked the sharpest of up to 40 frames).

## Display

Dye fills span a large range (a bright soma, dim processes). The projections are stored in their
original intensity units; the cell figures of `acq4_mosaic` show them on a square-root scale.

## References

P. J. Burt and E. H. Adelson, "Merging images through pattern decomposition," Proc. SPIE 575
(1985), doi:10.1117/12.966501.

P. J. Burt and R. J. Kolczynski, "Enhanced image capture through fusion," Proc. 4th
International Conference on Computer Vision (1993), 173–182, doi:10.1109/ICCV.1993.378222.

B. Forster, D. Van De Ville, J. Berent, D. Sage and M. Unser, "Complex wavelets for extended
depth-of-field: a new method for the fusion of multichannel microscopy images," Microscopy
Research and Technique 65 (2004), 33–42, doi:10.1002/jemt.20092.
