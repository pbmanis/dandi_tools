# Slice mosaics (acq4_mosaic)

*Experimental.* `acq4_mosaic` assembles the images and videos of one slice into an aligned
mosaic, finds the recorded cells on it, and writes figures, an acq4 Mosaic Editor file and a
report that carries each image's metadata forward. It is meant for checking the placement of
slices and cells, and as input for aligning slices and cells to an atlas (e.g. with the
cochlear-nucleus alignment tool in ephys).

## Running it

```
uv run acq4_mosaic <slice_dir> [<slice_dir> ...] -o <output_dir>
uv run acq4_mosaic --random 4 --root <acq4 data root> -o <output_dir>
```

| Option | Meaning |
|---|---|
| `-o`, `--output` | Output directory (default `mosaic_test`) |
| `--random N --root R` | N cells chosen at random under root(s) R (repeat `--root` for several) |
| `--seed` | Random seed (printed when not given) |
| `--max-video-mb` | Leave out larger videos (default 500 MB; usually recordings that were not stopped) |
| `--video-projection` | How a video becomes one image: `focus` (default), `mip` or `best` ([](video_stacks.md)) |
| `--cell-region-um` | Width of the region around each cell in the cell figures (default 1500 µm) |
| `--dataset NAME` | Take settings from that dataset's entry in the configuration file: `mosaic_cell_region_um` (e.g. 1200-1500 for cochlear nucleus, ~2000 for cortex) |
| `--config` | The configuration file for `--dataset` (default `config/dandi_configuration.cfg`) |
| `--decisions FILE` | An exceptions workbook from an earlier run, with decisions filled in (see below) |

The raw data are only read.

## What it does

1. **Find the images.** `image_NNN.tif` and `video_NNN.ma` in the slice directory and its
   subdirectories (not protocol directories), with position, pixel size, objective and the
   manual correction (`userTransform`) from the `.index` files. Folders named garbage, bad, etc.
   are left out.
2. **Channels.** From folder names: FITC, GFP, YFP → green; TRITC, mCherry, tdTomato, DiI →
   red; Lucifer yellow, Alexa → yellow. Everything else is treated as brightfield. (Dye-fill
   videos in cell folders are therefore labelled brightfield.)
3. **Videos.** Each video (a z-stack) is projected to one image ([](video_stacks.md)).
4. **Flat field.** The illumination pattern of each set of 5 or more images of the same channel,
   magnification and size is divided out.
5. **Objective check.** An image that matches images in its folder labelled with another
   objective (and fewer with its own) is given that objective and placed from its best match.
6. **Base.** The slice-level folder with the most brightfield images at the lowest
   magnification is stitched: every pair of neighbouring tiles is registered (normalized
   cross-correlation of band-passed images, within 200 µm), each pair weighted, direction by
   direction, by how sharply its correlation peak is defined, and all translations are solved
   for together by least squares, dropping the most inconsistent pair until the rest agree.
7. **Groups.** The other images, one folder, channel and magnification at a time (lowest
   magnification first), are registered against everything placed so far: first as a group in a
   large window, then each image in a small one. Fluorescence is compared with brightfield by
   gradient magnitude. If a group does not match at its recorded magnification, the rig's other
   objectives are tried. Placements that are not clearly better than the alternatives are not
   used, and the images stay at their stage positions; the images of such a group are then
   stitched to each other (placement `stage+stitched`), so that, for example, the overlapping
   40x z-stacks of a filled cell are consistent with each other even though they cannot be
   registered to a 4x base.
8. **Joint solve.** All slice-level low-magnification images are stitched together.
9. **Cells.** From a Mosaic Editor file in the slice directory (`CellCanvasItem`, placed by
   hand), otherwise the centre of the cell's highest-magnification image (inferred).

Only translations are used (the stage only translates). A group move of more than 250 µm is
refused (the group stays at its stage position) and listed in the exceptions workbook: such
moves are usually the slice having been moved between images. The registration still searches
up to 1.5 mm, so that a clear far match is found and refused rather than a weak near one being
accepted. The constants of the registration are
explained in [](<../registration methods.md>).

## Output

For each slice, named `<set>~<day>~<slice>~...` in the output directory:

| File | Contents |
|---|---|
| `~mosaic.png` | The slice at the stage positions and as registered; higher-magnification images outlined; cells marked (magenta circles: placed by hand; blue squares: inferred) |
| `~qc.png` | For each moved group, its edges (magenta) on the base's edges (green), before and after; white where they coincide |
| `~cell_NNN.png` | 1500 µm (configurable) around each cell, every cell labelled, and the cell's images at higher magnification than the slice images, in file-number order (images, then videos), in 1-3 columns. Video projections on a square-root scale. Where the cell position is outside an image, the panel says how far away and in which direction. |
| `~aligned.mosaic` | The placements as an acq4 Mosaic Editor file (see below) |
| `~projections/` | The video projections, float32 TIFFs, with their placement in the TIFF metadata |
| `~report.json` | Everything: placements, corrections, problems, matches, cells, and each image's metadata and appearance measurements |
| `exceptions_<date>-<time>.xlsx` | One per run, for all of its slices: the exceptions to check by hand (see below) |

Every figure is stamped with the repository, its version and git commit at the top left
("+ uncommitted changes" when the code differs from the commit) and the processing date and
time at the bottom right.

### Metadata carried forward

Each entry of `tiles` in the report holds:

* `index`: the image's whole `.index` entry, as recorded (objective, binning, exposure, region,
  transform, deviceTransform, userTransform, time, ...);
* `corrected`: the placement made here, in the same terms:
  * `mag_index`, `mag_used`: the magnification in the `.index` and the one used (they differ when
    the objective was corrected);
  * `pixel_m`: the pixel size used (m);
  * `bounds_um`: the image's placed extent (xmin, xmax, ymin, ymax);
  * `transform`: the `.index` transform with the corrections applied (pos in m, scale in m per
    pixel with the original signs, angle). It places the image directly, with no
    userTransform;
  * `mosaic_userTransform`: the acq4 userTransform that moves the image from its `.index`
    transform to its placement here (it replaces the `.index` userTransform; see below);
* `correction_um` and `manual_correction_um`: this program's correction from the stage
  position, and the manual one in the `.index`;
* `projection` (videos): how the frames were projected.

### The .mosaic file

`~aligned.mosaic` is an acq4 Mosaic Editor save file (`MosaicEditor_save`, version 2.0). Each
image is an `ImageCanvasItem` whose `userTransform` places it as registered. In acq4 an item's
placement is its `.index` transform followed by the item's userTransform, and the item's
userTransform *replaces* the `.index` userTransform (the editor loads one or the other, and saves
the same value to both). A corrected magnification is a scale of the userTransform about the world
origin, so for those images `pos` can be large. Cells are `CellCanvasItem`s, with inferred
positions named `cell_NNN (inferred)`. Landmark (`MarkersCanvasItem`) and ruler items from the
slice's own Mosaic Editor files are copied unchanged. Videos are referenced as the original
`.ma` files (the editor shows their frames, not the projection).

`rootPath` is the slice directory and file names are relative to it. The acq4 Mosaic Editor
resolves file names relative to the directory that holds the `.mosaic` file, so to open the
file there, copy it into the slice directory. Readers that use `rootPath` (such as the mosaic
importer of the ephys cochlear-nucleus alignment tool) can read it where it is written.

## Exceptions and decisions

Each run writes `exceptions_<date>-<time>.xlsx` with a README sheet and one sheet per kind of
exception, each with a `decision` column (a drop-down list; empty means the default):

| Sheet | Lists | Decisions |
|---|---|---|
| image kind | Images at 10x or more whose appearance is ambiguous (or whose kind was set in the table), with a thumbnail and the measurements | `BF`, `fluorescence` (keeps a fluorescence folder's colour, otherwise yellow), `green`, `red`, `yellow`, `mixed` (fluorescence with brightfield for context), `exclude` |
| image kind (automatic) | Images whose channel was set from a confident classifier verdict, overriding the folder name, for checking | `BF`, `fluorescence` (keeps a fluorescence folder's colour, otherwise yellow), `green`, `red`, `yellow`, `mixed` (fluorescence with brightfield for context), `exclude` |
| refused matches | Group moves of more than 250 µm that were not used | `stage` (default), `accept` (use the move), `exclude` (leave the group out), `manual` (stage position, listed for a manual mosaic) |
| annotation mismatches | Cell positions placed by hand that are outside all of the cell's own higher-magnification images | `annotation` (default), `inferred` (centre of the cell's highest-magnification image), `exclude` |

Give the edited file to the next run with `--decisions`. The workbook written by that run lists
the same exceptions with the decisions filled in (including those the decisions resolved), so
one file can be kept up to date. Images can be refused for many reasons (the experimenter forgot
to turn the illumination on, the slice moved in the chamber, ...): use `exclude`.

### Brightfield or fluorescence: the classifier and the labelling page

`src/image_kind.py` measures each image's appearance (offset- and gain-independent) and a random
forest gives the probability that it is fluorescence. Between 0.3 and 0.8 the image is called
ambiguous and is left to you in the workbook. The measurements that carry the model are those
of relief: brightfield at 40-63x shows each edge as a bright/dark pair, so its band-passed image
goes as deep into the dark as into the bright (`dark_excursion` near 1, `structured_skew` near
0), while fluorescence features are bright only. Next come the bright tail relative to the local
background (`tail_ratio`) and the illumination gradient.

Fitted on 283 labelled images at 10x and above (234 fluorescence, 49 brightfield; 360 labelled
in two rounds and a review of the ambiguous cases, from 40 slices of the three datasets; also 42
mixed, 33 skip and 2 unsure, which are not used for training). Checked by slice: the slices are
divided into 10 groups and each image is predicted by a model fitted without any image of its
slice; mean of 5 divisions:

| Label | called fluorescence | called ambiguous | called BF |
|---|---|---|---|
| fluorescence (234) | 214.0 | 19.0 | 1.0 |
| BF (49) | 2.6 | 13.4 | 33.0 |

Mixed images (fluorescence taken with the room lights on, or both illumination paths on at
uncontrolled levels) are not a separate class: a three-class model could not recognise them (at
best 9 of 42). With the two-class model they come out fluorescence (21 of 42) or ambiguous
(17); as fluorescence they are registered by gradient magnitude, which suits mixed contrast.
Skip images (noise) mostly come out ambiguous (26 of 33), so they reach the workbook, where
`exclude` removes them. A logistic regression on the same measurements made more confident
errors; a first model fitted on round 1 only was tested on the 200 new images of round 2: of the
100 it called fluorescence none was brightfield.

**Default:** for images at 10x and above, a confident verdict sets the channel, overriding the
folder name: fluorescence (p >= 0.8) in a non-fluorescence folder becomes yellow, brightfield
(p < 0.3) in a fluorescence folder becomes BF. These are listed in the workbook sheet
`image kind (automatic)` for checking; ambiguous images keep the folder rule and are listed in
`image kind`. Decisions in either sheet override the classifier. `--no-auto-kind` turns this
off.

Labelling and refitting:

```
uv run python -m src.image_kind review <slice_dir> ... -o <dir> [--labels labels.csv] [--new 200 --low-p-fraction 0.5]
uv run python -m src.image_kind fit --labels labels.csv --measurements <dir>/measurements.csv
```

`review` writes `review.html` (thumbnails, click to switch linear and square-root display;
buttons BF, fluorescence, mixed, unsure and skip, for noise or otherwise unusable images, which
are shown with a red X and left out of the training; outline orange until labelled, then blue;
"Next unlabelled", "Download labels (CSV)", "Load labels (CSV)") and `measurements.csv`. Labels
are matched to images by path, so they survive rebuilding the page. With `--labels`, every
labelled image is on the page with its label; `--new N` adds N unlabelled images, of which the
fraction `--low-p-fraction` are those the classifier finds least likely to be fluorescence (to
find more brightfield examples) and the rest are random. `fit` refits the model (random forest by
default, `--method logistic` for the alternative), prints the by-slice validation and the images
misclassified in most divisions (candidates for checking the labels), and saves it to
`src/image_kind_model.json` and `src/image_kind_forest.pkl`.
