# THEMIS TDS event dataset

Local plasma and field conditions at THEMIS during wave bursts containing Time Domain Structures (Lang Qin's event list, `data/raw/AllData`), with solar-wind context to study TDS global distribution and infer their origin.

Methods, product choices, caveats and column definitions: `report/report.typ`.

## Using the dataset

The built table is attached to each GitHub release (`tds_tha.parquet`, `tds_tha.csv`, `columns.csv`, `report.pdf`); one row per
list interval.

```python
import json, pandas as pd, pyarrow.parquet as pq
df = pd.read_parquet("tds_tha.parquet")
meta = pq.read_schema("tds_tha.parquet").metadata
cols = json.loads(meta[b"columns"])          # {column: {"unit", "source" ("Lang" or "ours"), "description"}}
print(meta[b"pipeline"].decode(), cols["efw_Epar"])
```

With CSV, parse the time columns: `pd.read_csv("tds_tha.csv", parse_dates=["start_time", "end_time", "tds_start", "tds_end", "efw_t_Ez", "efw_t_Epar"])`.

## Pipeline

```sh
scripts/download.sh l2/gmom; scripts/download.sh l2/fgm; scripts/download.sh l1/efw   # event-day CDFs → data/cdf/ (~220 GB)
julia --project=. -t auto scripts/build.jl       # → data/tds_tha.csv (+ typed .parquet with column dictionary), data/columns.csv (~1 min)
julia --project=. -t 8 detect/run.jl             # TDS detector → detect/out/ (needs region from tds_tha.csv)
julia --project=. -t auto scripts/build.jl       # again: tds_* columns, and efw_* stats over the TDS windows
julia --project=. -t auto scripts/build.jl --sample[=24]  # rebuild days spread over the mission, print differences from data/tds_tha.csv (~15 s)
julia --project=. -t auto scripts/bgwindow.jl    # → data/bgwindow.csv (background-window test)
julia --project=. scripts/qa.jl                  # → report/figs/*.png, report/summary.json (numbers quoted in report.typ)

typst compile --root . report/report.typ
typst compile --root . detect/report.typ          # automatic TDS detection experiment
```

| file | role |
|---|---|
| `scripts/events.jl` | parse event list (`*-QLP.txt` replaces same-day plain file) |
| `scripts/themis.jl` | FGM + GMOM event/background averages |
| `scripts/efw.jl` | L1 EFW calibration and per-event E stats |
| `scripts/omni.jl` | OMNI 1-min + Kp |
| `scripts/tdswindow.jl` | detector TDS trains inside each interval (`det_status`, `tds_*`) |
| `scripts/lang_check.jl` | reproduce Lang's window, B, V, `max_Ez`, `rms_Ew` from his code's definitions on random days → report/lang_check.json |
| `idl/` | Reference selection code (`findthmbrst_v2.pro`) and an overview script (FGM, ESA, FFF spectra) |
| `scripts/build.jl` | join, derived quantities, quality masks, column dictionary (`COLUMNS`) |
| `scripts/bgwindow.jl` | sensitivity of b̂ and V∥ to the background-window length |
| `scripts/table.jl` | load `data/tds_tha.parquet` |
| `scripts/qa.jl` | dataset-quality figures and `report/summary.json` |
