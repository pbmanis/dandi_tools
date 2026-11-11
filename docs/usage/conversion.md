# From acq4 to the DANDI archive

The workflow has four steps: summarize the acq4 data, describe the dataset, convert, then check
and upload.

```
acq4 data ──dataSummary (ephys)──▶ datasummary database (pandas pickle)
                                          │
config/dandi_configuration.cfg ───────────┤
ephys experiments file (exclusions) ──────┤
                                          ▼
                             acq4_to_nwb <dataset> -c   ──▶  <output_dir>/<dataset_id>/*.nwb
                             acq4_to_nwb <dataset> -v   (nwbinspector, DANDI configuration)
                             acq4_to_nwb <dataset> -u   (download dandiset.yaml, organize, validate)
                             dandi upload               (by hand)
```

## 1. The acq4 data

acq4 stores an experiment as a directory tree: the day (`2021.04.05_000`), the slice
(`slice_000`), the cell (`cell_000`) and the protocol (`CCIV_long_000`), with one subdirectory
per sweep (`000`, `001`, ...). Each sweep holds a MetaArray file for each device (e.g.
`MultiClamp1.ma`); `.index` files in every directory hold the metadata (subject, solutions,
amplifier settings, stimulus parameters, image positions).

## 2. The datasummary database

The conversion is driven by the datasummary database made by the `dataSummary` tool in
[ephys](https://github.com/pbmanis/ephys). It scans the acq4 data and writes a pickled pandas
table with one row per cell; the `data_complete` column lists the protocols that were completed.
The ephys experiments file gives the database's location (`databasepath`, `directory`,
`datasummaryFilename`) and the protocols to leave out (`excludeIVs`, by day/slice or
day/slice/cell; `"all"` leaves out every protocol).

## 3. The dataset configuration

Each dataset is an entry in `config/dandi_configuration.cfg` (pyqtgraph configfile format;
see `config/dandi_example_configuration.cfg`):

```
DCN_IC_inj:
    dataset_id: "001423"
    keywords: [...]
    experimenters: [...]
    grants: "..."
    experiment_description: "..."
    iacuc_protocols: "..."
    device: "MultiClamp1.ma"      # acq4 amplifier file
    laser_power_W: 0.070          # nominal laser power (W)
    led_power_W: 0.010            # nominal LED power (W); only for LED datasets
    configfile: ".../experiments.cfg"   # the ephys experiments file
    experiment: "DCN_IC_inj"      # its experiment name, if different
    output_dir: "NWBFiles"        # the dataset id is added
```

Leave `led_power_W` out for datasets without LED protocols (the power is then recorded as NaN,
with a warning); do not set it to 0.0, which would record the LED pulses as 0 W.

## 4. Convert

```
uv run acq4_to_nwb <dataset> -c
```

Each cell becomes one NWB file, `<output_dir>/<dataset_id>/ExperimentName~Rig#~yyyy.mm.dd~SnnCmm.nwb`.
Only protocols of known types are converted (current-clamp and voltage-clamp IV series, LED
stimulation, laser scanning photostimulation maps); the recording mode of each protocol is read
from the amplifier settings. Within a file:

* each sweep is a stimulus and a response series, and a row of the intracellular and
  simultaneous recordings tables;
* the sweeps of a protocol are a sequential recording; repeats of a protocol
  (`CCIV_long_000`, `CCIV_long_001`) are a repetition;
* optical stimulation is stored in OptogeneticStimulusSite and OptogeneticSeries structures
  (one site per spot for laser scanning, one per protocol for the LED).

Options: `-t` reads the data and metadata without writing files; `-o` overwrites existing
files (otherwise cells with a file are skipped).

## 5. Check

```
uv run acq4_to_nwb <dataset> -v           # nwbinspector on every file, DANDI configuration
uv run acq4_to_nwb <dataset> -s <file>    # one file
```

or by hand, `nwbinspector <folder> --config dandi`. Errors marked CRITICAL should be fixed in the
acq4 data (the `.index` files) and the files converted again. Then `dandi validate <folder>`.

## 6. Upload

```
uv run acq4_to_nwb <dataset> -u
```

downloads the dandiset metadata (`dandiset.yaml`), organizes the converted files into the
dandiset directory (dry run first) and validates, stopping if a step fails. The upload itself,
`dandi upload` (or `dandi upload -i dandi-staging`), is run by hand.

## Images

The NWB files do not yet hold the slice and cell images. `acq4_mosaic` assembles them into
aligned mosaics with their metadata (see [](mosaics.md)); storing them in the archive is planned.
