# Installation

The project uses [uv](https://docs.astral.sh/uv/) and `pyproject.toml`. To create (or update) the
environment in `.venv`:

```
uv sync
```

This installs the runtime dependencies and the `dev` group (black, pytest, mypy, sphinx,
myst-parser, ...). Use `uv sync --no-dev` for the runtime dependencies only. Run the tools with
`uv run`, for example `uv run acq4_to_nwb --help`, or activate the environment with
`source .venv/bin/activate`.

The acq4 reader (`src/acq4_reader.py`) is included. `src/display_acq4.py` needs the
[ephys](https://github.com/pbmanis/ephys) package, which is not installed by default.

## Command line tools

| Command | What it does |
|---|---|
| `acq4_to_nwb` | Convert a dataset to NWB, check the files, prepare the DANDI upload ([](conversion.md)) |
| `display_nwb` | Plot each protocol of one NWB file ([](viewing.md)) |
| `browse_nwb` | Browse the NWB files of a directory by day and cell ([](viewing.md)) |
| `acq4_mosaic` | Assemble and align the images of a slice ([](mosaics.md)) |

## Building this documentation

```
uv run sphinx-build -b html docs docs/_build/html
```

then open `docs/_build/html/index.html`.
