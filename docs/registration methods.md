# Registration methods: the constants in `acq4_mosaic`

This page explains the fixed constants in the two functions that prepare images for registration in
`src/acq4_mosaic.py`: `resampled()` and `features()`. None of these values was fitted to the data.
They are standard image-processing choices, and as far as the code shows, no other values have been
tested. The overall method (pairwise registration, least-squares stitching, group registration) is
described in the module docstring of `acq4_mosaic`.

## `resampled()`

`resampled()` returns a tile's image, or its registration features, at a given resolution (um per
pixel), with a mask. Results are cached in `_resample_cache`.

| Constant | Purpose |
|---|---|
| `round(res, 6)` in the cache key | The resolution is a float worked out in several places. Rounding it means 2.0000000001 and 2.0 find the same cached entry. `t.mag` is also in the key because the objective-correction step can change a tile's magnification and pixel size, which makes an older resampled copy wrong. |
| `zoom = t.pixel / res` | The scale factor. Below 1 means the tile is being shrunk. |
| `gaussian_filter(a, 0.5 / zoom)` | Smooths the image before shrinking so that fine detail does not turn into false patterns (anti-aliasing). The width is about half an output pixel, measured in source pixels; shrinking by 4× gives σ = 2 source pixels. It is a deliberately light filter, and it is skipped when enlarging. |
| `order=1` | Linear interpolation. Cubic interpolation would add overshoot ("ringing") at sharp edges, such as cell edges or a pipette tip, and that adds false edges for registration. |
| `border = 3` | Hides 3 pixels around each tile's edge in feature mode. The filters in `features()` produce false structure at the edge, and the hidden strip is about twice the 1–1.5 px filter widths. Without it, the tile outline itself could be matched as if it were an image feature. The hidden strip does not cover the edge effects of the wider 10 px background filter. Those change slowly and mostly just shift the local level, so they are left unmasked. |

## `features()`

`features()` turns an image into a band-passed, normalized feature image for normalized
cross-correlation.

All filter widths are in pixels **at the registration resolution**, not in microns. Registration runs
at 2× the native pixel size or coarser (in `register_group()`), so in microns the filters grow with
the resolution in use. In the objective check (`check_objectives()`), the image has already been
reduced 4×.

- **`np.percentile(a, 0.5)`**: the dark floor. The 0.5th percentile is used instead of the minimum so
  that a few dead or cold pixels do not set it.
- **`log1p(... / (median − lo))`**: dividing by the median brightness above the floor makes tiles
  taken with different exposure or gain comparable. The log then compresses bright regions (glare,
  bright fluorescent cells) so they do not dominate. `max(..., 1.0)` (in camera counts) prevents
  dividing by about zero for a blank or very dark tile.
- **DoG σ = 1 and 10** (mode `"dog"`): a band-pass filter. The 1 px blur removes pixel noise.
  Subtracting the 10 px blur removes slow trends such as vignetting, uneven illumination and focus
  changes across the field. What is left is mid-sized structure (cell bodies, vessels, slice edges),
  which is what should be aligned. This mode is used only when both sides are brightfield.
- **`gaussian_gradient_magnitude(a, 1.5)`** (mode `"grad"`): used when the channels differ
  (brightfield vs fluorescence). The same structure can look dark in one and bright in the other, and
  gradient magnitude ignores which way the contrast goes. The width is 1.5 rather than 1 because
  taking a derivative amplifies noise.
- **`g − gaussian_filter(g, 10)`**: removes the regional level of edge strength, so that a noisy or
  textured area does not outweigh a clean edge just by having more edge signal overall.

### The last two lines: a robust standardization

```python
mad = np.median(np.abs(f - med)) * 1.4826 + 1e-9
return np.clip((f - med) / mad, -6.0, 6.0)
```

- **Median and MAD (median absolute deviation) instead of mean and standard deviation.** The things
  that should not drive the match (saturated spots, dust, the electrode, a very bright cell) would
  also inflate the mean and SD. Median and MAD mostly ignore them.
- **`1.4826`** = 1/Φ⁻¹(0.75). For Gaussian noise, MAD × 1.4826 equals the standard deviation. The
  output is then in "sigma" units, so tiles end up on a common scale whatever their channel or
  exposure.
- **`1e-9`** prevents dividing by zero when MAD is 0, which happens with a flat or fully saturated
  tile. Such a tile then just gives zeros.
- **`clip(-6, 6)`** caps any pixel at ±6 robust sigma. Normalized cross-correlation (the matching
  score) adds up products of pixel values. A few huge pixels, like a saturated spot or the pipette
  tip, could otherwise decide the best offset by themselves. Gaussian noise almost never gets past 6σ
  (about 2 in 10⁹), so the cap only affects real outliers and leaves ordinary structure alone. The
  choice of 6 rather than 4 or 8 is a judgment call.

## Sensitivity

If the sensitivity to these choices is examined, the two most likely to matter are the 10 px
background width (it sets the largest structure that counts) and the clip level. The clip level
matters most for sections with strong fluorescence hot spots.
