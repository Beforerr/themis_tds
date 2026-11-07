# THEMIS TDS event dataset

Local plasma and field conditions at THEMIS during wave bursts containing Time Domain Structures (Lang Qin's event list, `data/raw/AllData`), with solar-wind context to study TDS global distribution and infer their origin.

Methods, product choices, caveats and column definitions: `report/report.typ`.

## Using the dataset

The built tables are attached to each GitHub release: `tds_tha.parquet` (also `tds_tha.csv`, `columns.csv`), one row per
list interval; `tds_tha_pulses.parquet`, one row per TDS pulse; and `report.pdf`.

```python
import json, pandas as pd, pyarrow.parquet as pq
df = pd.read_parquet("tds_tha.parquet")
meta = pq.read_schema("tds_tha.parquet").metadata
cols = json.loads(meta[b"columns"])          # {column: {"unit", "source" ("Lang" or "ours"), "description"}}
print(meta[b"pipeline"].decode(), cols["efw_Epar"])
pulses = pd.read_parquet("tds_tha_pulses.parquet")   # same metadata keys; join on id
```

With CSV, parse the time columns: `pd.read_csv("tds_tha.csv", parse_dates=["start_time", "end_time", "tds_start", "tds_end", "efw_t_Ez", "efw_t_Epar"])`.

## Pipeline

```sh
julia --project=. scripts/dates.jl               # → data/dates.txt: event days and the neighbouring days their windows reach
scripts/download.sh l2/gmom; scripts/download.sh l2/fgm; scripts/download.sh l1/efw   # those days' CDFs → data/cdf/ (~220 GB)
julia --project=. -t auto scripts/build.jl       # → data/tds_tha.csv (+ typed .parquet with column dictionary), data/columns.csv (~1 min)
julia --project=. -t 8 detect/run.jl             # TDS detector → detect/out/ (needs region from tds_tha.csv)
julia --project=. -t auto scripts/build.jl       # again: tds_* columns, efw_* stats over the TDS windows, data/tds_tha_pulses.parquet
julia --project=. -t auto scripts/build.jl --sample[=24]  # rebuild days spread over the mission, print differences from data/tds_tha.csv (~15 s)
julia --project=. -t auto scripts/bgwindow.jl    # → data/bgwindow.csv (background-window test)
julia --project=. -t auto scripts/efwstages.jl   # → data/efw_stages.csv, data/efw_example.csv (waveform-cleaning diagnostics)
julia --project=. scripts/qa.jl                  # → report/figs/*.png, report/summary.json (numbers quoted in report.typ)

typst compile --root . report/report.typ
typst compile --root . detect/report.typ          # automatic TDS detection experiment
```

| file | role |
|---|---|
| `scripts/events.jl` | parse event list (`*-QLP.txt` replaces same-day plain file) |
| `scripts/dayfile.jl` | event-day CDF lookup (`dayfile`) |
| `scripts/daily.jl` | per-(probe, day) products (themis, efw, efw_stages, detect): all stale products of a day run in one pass over its files, results cached in `data/cache/<product>/` |
| `scripts/themis.jl` | FGM + GMOM event/background averages |
| `scripts/efw.jl` | L1 EFW calibration and per-event E stats |
| `scripts/omni.jl` | OMNI 1-min + Kp |
| `scripts/tdswindow.jl` | assigns detector TDS pulses to intervals: `det_status`, `tds_*`, `tds_tha_pulses.parquet` |
| `scripts/lang_check.jl` | reproduce Lang's window, B, V, `max_Ez`, `rms_Ew` from his code's definitions on random days → report/lang_check.json |
| `idl/` | Reference selection code (`findthmbrst_v2.pro`) and an overview script (FGM, ESA, FFF spectra) |
| `scripts/build.jl` | join, derived quantities, quality masks, column dictionary (`COLUMNS`) |
| `scripts/bgwindow.jl` | sensitivity of b̂ and V∥ to the background-window length |
| `scripts/efwstages.jl` | max-\|E∥\| spin phase without band-pass/mask, and one burst before/after cleaning |
| `scripts/table.jl` | load `data/tds_tha.parquet` |
| `scripts/qa.jl` | dataset-quality figures and `report/summary.json` |
