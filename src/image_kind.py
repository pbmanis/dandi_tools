"""Experimental: tell brightfield from fluorescence images by their appearance.

Folder names say which images are fluorescence only for some folders (FITC, TRITC, ...); images
and videos in the slice and cell folders can be either, or a mix (fluorescence with some
brightfield for context). This module measures the properties that tell them apart, following
the description of the images by the experimenter:

* Fluorescence (at 10x, 40x, 63x) has a few well-defined, local, bright, contiguous features on
  a darker background, darkest where the field includes areas that are not tissue. The features
  are filled somata (8-40 um, bright and featureless), filled processes (strands, from ~5 um
  near the soma to below 1 um), or retrogradely transported dye in lysosomes (spots of 1-3 um,
  clustered within a cell and usually not in the nucleus). There may be several filled cells.
* Brightfield often has an illumination gradient across the image (oblique illumination), and
  is granular all over.

The measurements (see :func:`measure`) do not depend on the camera offset or gain.
:func:`classify` combines four of them in a logistic model fitted to images labelled by the user
(:func:`fit`: a random forest, saved in ``image_kind_model.json`` and ``image_kind_forest.pkl``
next to this module), with a band of
probabilities in which the image is called ambiguous and left to the user. The measurements that
separate best are those of relief: brightfield at 40-63x shows each edge as a bright/dark pair,
so its band-passed image is about as deep in the dark as it is bright, while fluorescence
features are bright only.

Labelling: :func:`build_review` writes a page of thumbnails with buttons (BF, fluorescence,
mixed, unsure, skip); its labels are downloaded as CSV and given to :func:`fit`::

    python -m src.image_kind review <slice_dir> ... -o <dir> [--labels labels.csv]
    python -m src.image_kind fit --labels labels.csv --measurements <dir>/measurements.csv
"""

import json
from pathlib import Path

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np
import scipy.ndimage as ndi


@dataclass
class Appearance:
    """Measurements of one image that bear on brightfield vs fluorescence.

    Attributes
    ----------
    tail_ratio : float
        (p99.5 - p50) / (p50 - p0.5) of the image relative to its local background (the image
        divided by itself smoothed at 25 um, which removes vignetting and illumination
        gradients): how far the bright tail reaches, relative to the dark tail. Large for
        fluorescence (bright features, no dark ones); about 1 for brightfield (dark features as
        well as bright ones).
    raw_tail_ratio : float
        The same for the image itself (affected by vignetting and gradients).
    bright_area : float
        Fraction of the image in bright features (relative image, smoothed at 1 um, above the
        median by more than half of p99.9 - p50).
    n_bright : int
        Number of bright features of at least 1 um^2.
    largest_bright_um2 : float
        Area of the largest bright feature, um^2.
    dark_area : float
        Fraction of the image within 10% of the intensity range of the darkest level (dark
        background, e.g. outside the tissue in fluorescence).
    gradient : float
        Range of a plane fitted to the image smoothed at 20 um, relative to p99 - p1: a large
        smooth gradient (oblique brightfield illumination, but fluorescence images are often
        vignetted too, so it is recorded and not used by the rule).
    texture_spread : float
        Fraction of the relative image whose fine texture (band 0.5-2 um, averaged over 10 um) is at
        least half of the 95th percentile of the texture: near 1 when texture is everywhere
        (granular brightfield), small when it is concentrated in a few features (fluorescence).
    structured_skew : float
        Robust skewness, (p99.5 + p0.5 - 2 p50) / (p99.5 - p0.5), of the image band-passed at
        about 1 um (0.5-2 um), in the 20% of the image with the most structure (band energy
        averaged over 5 um): positive when features are bright only (fluorescence), near 0 when
        edges are bright/dark pairs (brightfield relief).
    dark_excursion : float
        Depth of the dark excursions relative to the bright ones, -p0.5 / p99.5, of the image
        band-passed at about 2 um (1-4 um), in the same structured areas: about 1 for
        brightfield, less for fluorescence.
    """

    tail_ratio: float
    raw_tail_ratio: float
    bright_area: float
    n_bright: int
    largest_bright_um2: float
    dark_area: float
    gradient: float
    texture_spread: float
    structured_skew: float = 0.0
    dark_excursion: float = 1.0

    def as_dict(self) -> dict:
        """The measurements as a dictionary (rounded), for reports.

        Returns
        -------
        dict
            Field name -> value.
        """
        return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in asdict(self).items()}


def measure(image: np.ndarray, pixel_um: float) -> Appearance:
    """Measure the appearance of an image (see :class:`Appearance`).

    The image is first reduced to no finer than 0.5 um per pixel, so the measurements mean the
    same at every magnification.

    Parameters
    ----------
    image : numpy.ndarray
        2D image (raw intensities, any offset and gain).
    pixel_um : float
        Pixel size, um.

    Returns
    -------
    Appearance
        The measurements.
    """
    a = np.asarray(image, dtype=np.float32)
    step = max(int(np.floor(0.5 / pixel_um)), 1)
    if step > 1:  # bin to about 0.5 um
        h, w = (a.shape[0] // step) * step, (a.shape[1] // step) * step
        a = a[:h, :w].reshape(h // step, step, w // step, step).mean(axis=(1, 3))
    px = pixel_um * step

    def um(v):
        """A length in um as a filter width in (binned) pixels.

        Parameters
        ----------
        v : float
            Length, um.

        Returns
        -------
        float
            Pixels (at least 0.5).
        """
        return max(v / px, 0.5)

    p0, p1, p50, p99, p995, p999 = np.percentile(a, [0.5, 1, 50, 99, 99.5, 99.9])
    spread = max(p99 - p1, 1e-6)
    raw_tail_ratio = float((p995 - p50) / max(p50 - p0, 1e-6 * spread))
    dark_area = float((a < p0 + 0.1 * (p995 - p0)).mean())

    # relative to the local background (above the darkest level): removes vignetting and gradients
    floor = p0 - 0.05 * spread
    background = ndi.gaussian_filter(a, um(25.0)) - floor
    r = (a - floor) / np.maximum(background, 0.05 * spread) - 1.0
    q0, q50, q995, q999 = np.percentile(r, [0.5, 50, 99.5, 99.9])
    tail_ratio = float((q995 - q50) / max(q50 - q0, 1e-6))

    smooth = ndi.gaussian_filter(r, um(1.0))
    p50, p999 = q50, q999
    level = p50 + 0.5 * (p999 - p50)
    bright = smooth > level
    labels, n = ndi.label(bright)
    sizes = ndi.sum(bright, labels, index=np.arange(1, n + 1)) * px * px if n else np.zeros(0)
    keep = sizes >= 1.0
    bright_area = float(bright.mean())

    coarse = ndi.gaussian_filter(a, um(20.0))
    yy, xx = np.mgrid[0 : a.shape[0], 0 : a.shape[1]]
    A = np.column_stack([np.ones(a.size), yy.ravel(), xx.ravel()])
    coef, *_ = np.linalg.lstsq(A, coarse.ravel(), rcond=None)
    plane = A @ coef
    gradient = float((plane.max() - plane.min()) / spread)

    band = ndi.gaussian_filter(r, um(0.5)) - ndi.gaussian_filter(r, um(2.0))
    energy = ndi.gaussian_filter(band * band, um(10.0))
    texture_spread = float((energy >= 0.5 * np.percentile(energy, 95)).mean())

    def structured_band(scale):
        """Band-pass of the image at about scale um, in its most structured 20%.

        Parameters
        ----------
        scale : float
            Centre of the band, um (the band is scale/2 to 2 x scale).

        Returns
        -------
        numpy.ndarray
            The band-passed values where the band energy (averaged over 5 um) is in its top 20%
            (all of them if the image has no structure).
        """
        b = ndi.gaussian_filter(a, um(scale / 2)) - ndi.gaussian_filter(a, um(scale * 2))
        e = ndi.gaussian_filter(b * b, um(5.0))
        top = b[e >= np.percentile(e, 80)]
        return top if top.size else b.ravel()  # a featureless image: use all of it

    b1 = structured_band(1.0)
    r1, r50, r99 = np.percentile(b1, [0.5, 50, 99.5])
    structured_skew = float((r99 + r1 - 2 * r50) / max(r99 - r1, 1e-9))
    b2 = structured_band(2.0)
    r1, r99 = np.percentile(b2, [0.5, 99.5])
    dark_excursion = float(-r1 / max(r99, 1e-9 * max(-r1, 1.0)))

    return Appearance(
        tail_ratio=tail_ratio,
        raw_tail_ratio=raw_tail_ratio,
        bright_area=bright_area,
        n_bright=int(keep.sum()),
        largest_bright_um2=float(sizes.max()) if len(sizes) else 0.0,
        dark_area=dark_area,
        gradient=gradient,
        texture_spread=texture_spread,
        structured_skew=structured_skew,
        dark_excursion=dark_excursion,
    )


model_file = Path(__file__).with_name("image_kind_model.json")
forest_file = Path(__file__).with_name("image_kind_forest.pkl")  # the fitted random forest (scikit-learn)
_model = None
_forest = None


def load_model() -> Optional[dict]:
    """The fitted model (``image_kind_model.json``), or None if there is none.

    Returns
    -------
    dict or None
        method ("forest" or "logistic"), features, log (which features are log-transformed),
        mean, scale, zone (the ambiguous band of probabilities), coef and intercept (logistic),
        and how it was trained and validated. A forest is loaded from ``image_kind_forest.pkl``.
    """
    global _model, _forest
    if _model is None and model_file.exists():
        _model = json.loads(model_file.read_text())
        if _model.get("method") == "forest":
            import pickle

            with open(forest_file, "rb") as f:
                _forest = pickle.load(f)
    return _model


def _vector(m: dict, model: dict) -> np.ndarray:
    """The model's features of one image, transformed and standardized.

    Parameters
    ----------
    m : dict
        Measurements (from :meth:`Appearance.as_dict` or a measurements.csv row).
    model : dict
        The model.

    Returns
    -------
    numpy.ndarray
        The standardized feature vector.
    """
    v = np.array([float(m[f]) for f in model["features"]])
    v = np.where(model["log"], np.log(np.maximum(v, 1e-3)), v)
    return (v - np.array(model["mean"])) / np.array(model["scale"])


def classify(m: Appearance) -> tuple:
    """Brightfield, fluorescence or ambiguous, by the fitted model.

    The model gives the probability that the image is fluorescence; below the lower limit of
    its ambiguous band the image is called BF, above the upper limit fluorescence, and in
    between ambiguous. Without a fitted model, a hand-set rule on the relief measurements is
    used (structured_skew > 0.25 and dark_excursion < 0.6: fluorescence; structured_skew < 0.1
    and dark_excursion > 0.9: BF).

    Parameters
    ----------
    m : Appearance
        The measurements.

    Returns
    -------
    kind : str
        "fluorescence", "BF" or "ambiguous".
    p_fluorescence : float
        Probability of fluorescence (NaN with the hand-set rule).
    """
    model = load_model()
    if model is None:
        if m.structured_skew > 0.25 and m.dark_excursion < 0.6:
            return "fluorescence", float("nan")
        if m.structured_skew < 0.1 and m.dark_excursion > 0.9:
            return "BF", float("nan")
        return "ambiguous", float("nan")
    v = _vector(m.as_dict(), model)
    if model.get("method") == "forest":
        p = float(_forest.predict_proba(v[None, :])[0, 1])
    else:
        z = float(v @ np.array(model["coef"]) + model["intercept"])
        p = 1.0 / (1.0 + np.exp(-z))
    lo, hi = model["zone"]
    kind = "fluorescence" if p >= hi else ("BF" if p < lo else "ambiguous")
    return kind, round(float(p), 3)


all_features = (
    "tail_ratio",
    "raw_tail_ratio",
    "bright_area",
    "n_bright",
    "largest_bright_um2",
    "dark_area",
    "gradient",
    "texture_spread",
    "structured_skew",
    "dark_excursion",
)
_logged = ("tail_ratio", "dark_excursion", "raw_tail_ratio", "bright_area", "largest_bright_um2")


def fit(
    labels_csv,
    measurements_csv,
    features: tuple = all_features,
    method: str = "forest",
    zone: tuple = (0.3, 0.8),
    splits: int = 5,
    save: bool = True,
) -> dict:
    """Fit the classifier to labelled images, check it by slice, and save it.

    Only images labelled BF or fluorescence are used (mixed, unsure and skip are counted but
    not used: skip marks images that are noise or otherwise unusable). Some features are
    log-transformed; all are standardized. The model is a random forest (default) or a logistic
    regression.

    Validation is by slice: the slices are divided into 10 groups, and each image is predicted
    by a model fitted without any image of its slice (images of one slice resemble each other,
    so leaving out one image at a time would flatter the model). This is repeated for
    ``splits`` random divisions of the slices and the counts are averaged. The images that are
    confidently misclassified in most divisions are listed (``disagreements``), as candidates
    for checking the labels.

    Parameters
    ----------
    labels_csv : str or pathlib.Path
        Labels downloaded from the review page (id, path, label).
    measurements_csv : str or pathlib.Path
        ``measurements.csv`` written with the review page (matched by path).
    features : tuple of str, optional
        Measurements used (default: all of them).
    method : {"forest", "logistic"}, optional
        Model (default "forest": 400 trees, at least 3 images per leaf).
    zone : tuple of float, optional
        Probabilities of fluorescence between which an image is called ambiguous (default
        (0.3, 0.8)).
    splits : int, optional
        Number of random divisions of the slices for validation (default 5).
    save : bool, optional
        Write the model to ``image_kind_model.json`` (and the forest to
        ``image_kind_forest.pkl``) (default True).

    Returns
    -------
    dict
        The model summary, with "trained_on" (counts), "validation" (mean counts by slice) and
        "disagreements" (confidently misclassified images).
    """
    import csv
    import pickle
    import time

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, cross_val_predict

    def make():
        """A new, unfitted model of the chosen kind."""
        if method == "forest":
            return RandomForestClassifier(n_estimators=400, min_samples_leaf=3, random_state=0)
        return LogisticRegression(max_iter=1000)

    labels = read_labels(labels_csv)
    rows = [r for r in csv.DictReader(open(measurements_csv)) if r["path"] in labels]
    use = [r for r in rows if labels[r["path"]] in ("BF", "fluorescence")]
    logged = [f in _logged for f in features]
    X = np.array([[float(r[f]) for f in features] for r in use])
    X = np.where(logged, np.log(np.maximum(X, 1e-3)), X)
    y = np.array([labels[r["path"]] == "fluorescence" for r in use]).astype(int)
    mean, scale = X.mean(axis=0), X.std(axis=0)
    Z = (X - mean) / scale
    slices = np.array(["/".join(r["path"].split("/")[:3]) for r in use])  # the slice title starts each path
    names = sorted(set(slices))
    lo, hi = zone
    counts, wrong = [], np.zeros(len(y))
    for seed in range(splits):
        order = np.random.RandomState(seed).permutation(len(names))
        groups = np.array([order[names.index(g)] for g in slices])
        p = cross_val_predict(make(), Z, y, cv=GroupKFold(n_splits=min(10, len(names))), groups=groups, method="predict_proba")[:, 1]
        f, b = p >= hi, p < lo
        counts.append([np.sum(f & (y == 1)), np.sum(~f & ~b & (y == 1)), np.sum(b & (y == 1)), np.sum(b & (y == 0)), np.sum(~f & ~b & (y == 0)), np.sum(f & (y == 0))])
        wrong += (b & (y == 1)) | (f & (y == 0))
    c = np.mean(counts, axis=0)
    final = make().fit(Z, y)
    model = {
        "method": method,
        "features": list(features),
        "log": logged,
        "mean": mean.tolist(),
        "scale": scale.tolist(),
        "zone": [lo, hi],
        "fitted": time.strftime("%Y-%m-%d %H:%M"),
        "labels": str(labels_csv),
        "trained_on": {k: int(sum(1 for r in rows if labels[r["path"]] == k)) for k in ("fluorescence", "BF", "mixed", "unsure", "skip")},
        "slices": len(names),
        "validation": {
            "method": f"by slice: 10 groups of slices, mean of {splits} random divisions",
            "fluorescence called fluorescence": round(float(c[0]), 1),
            "fluorescence called ambiguous": round(float(c[1]), 1),
            "fluorescence called BF": round(float(c[2]), 1),
            "BF called BF": round(float(c[3]), 1),
            "BF called ambiguous": round(float(c[4]), 1),
            "BF called fluorescence": round(float(c[5]), 1),
        },
        "disagreements": [
            {"path": r["path"], "label": labels[r["path"]], "times_wrong": int(w)} for r, w in zip(use, wrong) if w > splits / 2
        ],
    }
    if method == "logistic":
        model["coef"] = final.coef_[0].tolist()
        model["intercept"] = float(final.intercept_[0])
    if save:
        global _model, _forest
        model_file.write_text(json.dumps(model, indent=1))
        if method == "forest":
            with open(forest_file, "wb") as f:
                pickle.dump(final, f)
            _forest = final
        _model = model
    return model


def thumbnail_jpeg(image: np.ndarray, width: int = 320, gamma: float = 1.0) -> bytes:
    """A small JPEG of an image, stretched between its 0.5 and 99.8 percentiles.

    Parameters
    ----------
    image : numpy.ndarray
        2D image.
    width : int, optional
        Width of the thumbnail, pixels (default 320).
    gamma : float, optional
        Display gamma (default 1: linear; 0.5 shows dim structure next to bright).

    Returns
    -------
    bytes
        The JPEG file contents.
    """
    import io

    import matplotlib.image

    a = np.asarray(image, dtype=np.float32)
    step = max(a.shape[1] // width, 1)
    a = a[::step, ::step]
    lo, hi = np.percentile(a, [0.5, 99.8])
    g = np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) ** gamma
    buf = io.BytesIO()
    matplotlib.image.imsave(buf, g, cmap="gray", vmin=0, vmax=1, format="jpg", pil_kwargs={"quality": 85})
    return buf.getvalue()


_page = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>Image kind review</title>
<style>
body{font-family:-apple-system,Helvetica,sans-serif;margin:16px;background:#f4f4f4;color:#222}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:12px}
.card{background:#fff;border:4px solid #f08c00;border-radius:6px;padding:8px;position:relative}
.pic{position:relative}
.pic svg{position:absolute;inset:0;width:100%;height:100%;pointer-events:none;display:none}
.card.skipped .pic svg{display:block}
.card.skipped img{opacity:0.6}
.card.done{border-color:#1c64d8}
.card img{width:100%;cursor:pointer;background:#000}
.meta{font-size:11px;line-height:1.35;word-break:break-all}
.guess{font-size:11px;color:#555;margin-bottom:4px}
.labels{display:flex;flex-wrap:wrap;gap:4px}
.labels label{display:inline-flex;align-items:center;gap:6px;font-size:15px;padding:6px 8px;border-radius:5px;cursor:pointer;background:#eee}
.labels label:hover{background:#ddd}
.labels input{width:22px;height:22px;margin:0;cursor:pointer}
#bar{position:sticky;top:0;background:#f4f4f4;padding:8px 0;z-index:1;display:flex;flex-wrap:wrap;gap:10px;align-items:center}
.key{font-size:13px}.key b.o{color:#f08c00}.key b.b{color:#1c64d8}
button,.filebtn{font-size:14px;padding:5px 10px}
</style></head><body>
<div id="bar"><b>Image kind review</b>
<span class="key">Outline: <b class="o">orange</b> = not labelled yet, <b class="b">blue</b> = labelled.
"skip" (red X): noise or otherwise unusable, left out of the training.
Click a picture to switch between linear and square-root display.</span>
<span id="count"></span>
<button onclick="nextUnlabelled()">Next unlabelled</button>
<button onclick="download()">Download labels (CSV)</button>
<label class="filebtn">Load labels (CSV) <input type="file" accept=".csv" onchange="load(this.files[0])"></label>
</div>
<div class="grid" id="grid"></div>
<script>
const items = __ITEMS__;
const KEY = "image-kind-review-labels-by-path";
let labels = __LABELS__;  // labels given when the page was made, by path
try { Object.assign(labels, JSON.parse(localStorage.getItem(KEY) || "{}")); } catch (e) {}
function save(){ try { localStorage.setItem(KEY, JSON.stringify(labels)); } catch (e) {} count(); }
function count(){ const n = items.filter(it => labels[it.path]).length; document.getElementById("count").textContent = n + " of " + items.length + " labelled"; }
const grid = document.getElementById("grid");
function render(){
  grid.innerHTML = "";
  for (const it of items) {
    const card = document.createElement("div"); card.className = "card" + (labels[it.path] ? " done" : ""); card.id = "c" + it.id;
    if (labels[it.path] === "skip") card.classList.add("skipped");
    const pic = document.createElement("div"); pic.className = "pic";
    const img = document.createElement("img"); img.src = it.linear; img.dataset.alt = it.sqrt; img.loading = "lazy";
    img.onclick = () => { const s = img.src; img.src = img.dataset.alt; img.dataset.alt = s; };
    pic.appendChild(img);
    pic.insertAdjacentHTML("beforeend", '<svg viewBox="0 0 100 100" preserveAspectRatio="none"><line x1="3" y1="3" x2="97" y2="97" stroke="#e00" stroke-width="5"/><line x1="97" y1="3" x2="3" y2="97" stroke="#e00" stroke-width="5"/></svg>');
    card.appendChild(pic);
    const meta = document.createElement("div"); meta.className = "meta";
    meta.textContent = "#" + it.id + "  " + it.path + "  (" + it.mag + "x, " + it.kind + ")";
    card.appendChild(meta);
    const guess = document.createElement("div"); guess.className = "guess";
    guess.textContent = "folder says: " + it.folder + "   classifier says: " + it.guess;
    card.appendChild(guess);
    const lab = document.createElement("div"); lab.className = "labels";
    for (const k of ["BF", "fluorescence", "mixed", "unsure", "skip"]) {
      const l = document.createElement("label");
      const r = document.createElement("input"); r.type = "radio"; r.name = "r" + it.id; r.value = k;
      r.checked = labels[it.path] === k;
      r.onchange = () => { labels[it.path] = k; card.classList.add("done"); card.classList.toggle("skipped", k === "skip"); save(); };
      l.appendChild(r); l.appendChild(document.createTextNode(k)); lab.appendChild(l);
    }
    card.appendChild(lab); grid.appendChild(card);
  }
  count();
}
render();
function nextUnlabelled(){
  const it = items.find(it => !labels[it.path]);
  if (it) document.getElementById("c" + it.id).scrollIntoView({behavior: "smooth", block: "center"});
}
function parseCSV(text){
  const out = {};
  for (const line of text.split(/\r?\n/).slice(1)) {
    const m = line.match(/^[^,]*,"?([^"]*)"?,(.*)$/);
    if (m && m[2].trim()) out[m[1]] = m[2].trim();
  }
  return out;
}
function load(file){
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => { Object.assign(labels, parseCSV(reader.result)); save(); render(); };
  reader.readAsText(file);
}
function download(){
  const rows = ["id,path,label"];
  for (const it of items) rows.push(it.id + ',"' + it.path + '",' + (labels[it.path] || ""));
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([rows.join("\n") + "\n"], {type: "text/csv"}));
  a.download = "image_kind_labels.csv"; a.click();
}
</script></body></html>
"""


def read_labels(labels_csv) -> dict:
    """Labels from a CSV downloaded from the review page.

    Parameters
    ----------
    labels_csv : str or pathlib.Path
        The file (id, path, label).

    Returns
    -------
    dict
        path -> label, for the rows that have a label.
    """
    import csv

    return {r["path"]: r["label"] for r in csv.DictReader(open(labels_csv)) if r.get("label")}


def _round_robin(candidates: list, n: int, max_video_fraction: float = 1.0) -> list:
    """Take up to n candidates, in turn from each folder type, with a cap on videos.

    Parameters
    ----------
    candidates : list of tuple
        (folder type, slice directory, Tile), already shuffled.
    n : int
        How many to take.
    max_video_fraction : float, optional
        Largest fraction of videos among those taken (default 1: no cap).

    Returns
    -------
    list of tuple
        The candidates taken.
    """
    by_folder = {}
    for c in candidates:
        by_folder.setdefault(c[0], []).append(c)
    chosen, videos = [], 0
    while len(chosen) < n and any(by_folder.values()):
        for k in list(by_folder):
            while by_folder[k] and len(chosen) < n:
                c = by_folder[k].pop()
                if c[2].kind == "video" and videos >= max_video_fraction * n:
                    continue
                videos += c[2].kind == "video"
                chosen.append(c)
                break
    return chosen


def build_review(
    slice_dirs: list,
    output,
    min_mag: float = 10.0,
    max_items: int = 160,
    seed: int = 0,
    labels_csv=None,
    n_new: Optional[int] = None,
    low_p_fraction: float = 0.0,
    pool_factor: int = 3,
) -> list:
    """Measure the images of some slices and write a page for labelling them by hand.

    Images and videos (as focus projections) at ``min_mag`` or more are candidates. Every
    candidate already labelled in ``labels_csv`` is included, with its label shown. New images
    are added: a pool of ``pool_factor`` x ``n_new`` is drawn evenly from the folder types (cell
    folders, fluorescence folders, slice folders; at most a quarter videos) and measured; then
    ``low_p_fraction`` of the new images are those the classifier finds least likely to be
    fluorescence (to find more brightfield examples), and the rest are drawn at random from the
    pool. The images are shown in a random order.

    Writes ``review.html`` and ``measurements.csv``. The page shows thumbnails (click to switch
    linear and square-root display) with buttons BF, fluorescence, mixed, unsure and skip
    (noise or otherwise unusable: shown with a red X); cards are outlined orange until
    labelled, then blue. Labels are kept in the browser, can be loaded from and downloaded as
    CSV, and are matched to images by path, so they survive rebuilding the page.

    Parameters
    ----------
    slice_dirs : list of pathlib.Path
        Slice directories (those of the labelled images must be among them).
    output : pathlib.Path
        Output directory.
    min_mag : float, optional
        Smallest magnification included (default 10).
    max_items : int, optional
        Number of images on the page when ``n_new`` is not given (default 160).
    seed : int, optional
        Random seed (default 0).
    labels_csv : str or pathlib.Path, optional
        Labels from an earlier page: those images are included, with their labels.
    n_new : int, optional
        Number of new (unlabelled) images (default: max_items less the labelled ones).
    low_p_fraction : float, optional
        Fraction of the new images chosen as the least likely to be fluorescence (default 0).
    pool_factor : int, optional
        Size of the pool of new images measured, as a multiple of n_new (default 3).

    Returns
    -------
    list of dict
        One row per image on the page: id, path, magnification, kind, folder channel, the
        classifier's verdict and probability, and the measurements.
    """
    import base64
    import csv
    import random

    from src.acq4_mosaic import load_tile, scan_slice, slice_title

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    given = read_labels(labels_csv) if labels_csv is not None else {}
    candidates = []
    for slice_dir in slice_dirs:
        tiles, _ = scan_slice(Path(slice_dir))
        for t in tiles:
            if not t.excluded and np.isfinite(t.mag) and t.mag >= min_mag:
                folder = "fluor-folder" if t.channel_from == "folder" and t.channel != "BF" else ("cell" if t.cell else "slice")
                candidates.append((folder, Path(slice_dir), t))
    rng = random.Random(seed)
    rng.shuffle(candidates)

    def path_of(c):
        """Path of a candidate as shown on the page (slice title / image name).

        Parameters
        ----------
        c : tuple
            (folder type, slice directory, Tile).

        Returns
        -------
        str
            The path.
        """
        return slice_title(c[1]) + "/" + c[2].name

    labelled = [c for c in candidates if path_of(c) in given]
    missing = set(given) - {path_of(c) for c in labelled}
    if missing:
        print(f"   {len(missing):d} labelled images were not found among the slices given (e.g. {sorted(missing)[0]})")
    if n_new is None:
        n_new = max(max_items - len(labelled), 0)
    pool = _round_robin([c for c in candidates if path_of(c) not in given], n_new * pool_factor, max_video_fraction=0.25)

    measured = {}
    for c in labelled + pool:
        t = c[2]
        try:
            load_tile(t)
        except Exception as e:  # unreadable files are skipped
            print(f"   {t.path}: {e!r}")
            continue
        m = measure(t.data, t.pixel)
        kind, p = classify(m)
        measured[path_of(c)] = (c, m, kind, p, thumbnail_jpeg(t.data), thumbnail_jpeg(t.data, gamma=0.5))
        t.data = None
        print(f"   {len(measured):4d} {path_of(c)}: {kind} ({p})")
    new = [measured[path_of(c)] for c in pool if path_of(c) in measured]
    n_low = int(round(low_p_fraction * n_new))
    by_p = sorted(new, key=lambda x: x[3] if np.isfinite(x[3]) else 0.5)
    chosen_new = by_p[:n_low]
    rest = by_p[n_low:]
    rng.shuffle(rest)
    chosen_new += rest[: max(n_new - n_low, 0)]
    shown = [measured[path_of(c)] for c in labelled if path_of(c) in measured] + chosen_new
    rng.shuffle(shown)

    rows, items = [], []
    for i, (c, m, kind, p, lin, sq) in enumerate(shown):
        folder, slice_dir, t = c
        path = path_of(c)
        row = {"id": i, "path": path, "mag": t.mag, "kind": t.kind, "folder": f"{folder} ({t.channel})", "guess": kind, "p_fluorescence": p}
        row.update(m.as_dict())
        rows.append(row)
        items.append(
            {
                "id": i,
                "path": path,
                "mag": f"{t.mag:g}",
                "kind": "video (focus projection)" if t.kind == "video" else "image",
                "folder": row["folder"],
                "guess": kind + ("" if not np.isfinite(p) else f" (p fluorescence {p:.2f})"),
                "linear": "data:image/jpeg;base64," + base64.b64encode(lin).decode(),
                "sqrt": "data:image/jpeg;base64," + base64.b64encode(sq).decode(),
            }
        )
    page = _page.replace("__ITEMS__", json.dumps(items)).replace("__LABELS__", json.dumps(given))
    (output / "review.html").write_text(page)
    with open(output / "measurements.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"   {len(items):d} images on the page: {len(items) - len(chosen_new):d} labelled before, {len(chosen_new):d} new")
    return rows


def main():
    """Command line: ``python -m src.image_kind review ...`` or ``python -m src.image_kind fit ...``.

    ``review <slice_dir> ... -o <dir> [--labels labels.csv]`` writes the labelling page and the
    measurements; ``fit --labels labels.csv --measurements measurements.csv`` fits the
    classifier, prints its leave-one-out validation and saves it.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Experimental: brightfield vs fluorescence, labelling and fitting")
    sub = parser.add_subparsers(dest="command", required=True)
    rv = sub.add_parser("review", help="measure images and write a page for labelling them")
    rv.add_argument("slices", nargs="+", type=Path, help="slice directories")
    rv.add_argument("-o", "--output", type=Path, default=Path("mosaic_test/image_kind"), help="output directory")
    rv.add_argument("--labels", type=Path, default=None, help="labels (CSV from the page) to show already set")
    rv.add_argument("--max-items", type=int, default=160)
    rv.add_argument("--new", type=int, default=None, help="number of new (unlabelled) images")
    rv.add_argument(
        "--low-p-fraction", type=float, default=0.0, help="fraction of the new images chosen as least likely fluorescence"
    )
    rv.add_argument("--min-mag", type=float, default=10.0)
    ft = sub.add_parser("fit", help="fit the classifier to labelled images")
    ft.add_argument("--labels", type=Path, required=True, help="labels (CSV from the page)")
    ft.add_argument("--measurements", type=Path, required=True, help="measurements.csv written with the page")
    ft.add_argument("--method", choices=["forest", "logistic"], default="forest")
    args = parser.parse_args()
    if args.command == "review":
        build_review(
            args.slices,
            args.output,
            min_mag=args.min_mag,
            max_items=args.max_items,
            labels_csv=args.labels,
            n_new=args.new,
            low_p_fraction=args.low_p_fraction,
        )
    else:
        model = fit(args.labels, args.measurements, method=args.method)
        print(json.dumps({k: model[k] for k in ("trained_on", "slices", "validation", "disagreements")}, indent=1))
        print(f"saved {model_file}")


if __name__ == "__main__":
    main()
