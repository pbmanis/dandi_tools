# Viewing NWB files

## display_nwb

```
uv run display_nwb <file.nwb>
```

plots each protocol of one file in its own window, in recording order. The layout:

| Panel | Contents |
|---|---|
| Upper left | All sweeps, stacked (responses: V in current clamp, I in voltage clamp) |
| Middle left | The commands (current or voltage), superimposed |
| Lower left | The light-source commands (LED or laser), if any |
| Upper right | Map of the laser stimulus positions, coloured by the response in the 20 ms after the light onset (laser scanning only). Click a spot to show its trace. |
| Middle right | One sweep (the first, or the one clicked in the map) |
| Lower right | The mean response |

Units are volts and amperes. Protocols without optical stimulation leave the map and the
light-command panels empty.

## browse_nwb

```
uv run browse_nwb [directory]
```

opens a window listing the NWB files under a directory (searched recursively, so a dandiset
directory with one subdirectory per subject works), grouped by recording day (from the file
name), with the session and subject information of the selected file and a plot of one
protocol (the display_nwb layout).

| Key | Action |
|---|---|
| Ctrl+Right / Ctrl+Left | Next / previous protocol (continues into the next / previous file) |
| Ctrl+Down / Ctrl+Up | Next / previous file |
| PgDown / PgUp | Next / previous day |

## display_acq4

`src/display_acq4.py` plots one acq4 protocol directory directly, before conversion (for
comparison). It needs the ephys package.
