# Load the built event table with its column dictionary and pipeline version, plus plotting helpers shared by
# scripts/qa.jl and downstream analyses.
using DataFrames, Parquet2, JSON, Statistics, CairoMakie, AlgebraOfGraphics

const TABLE = joinpath(@__DIR__, "..", "data", "tds_tha.parquet")
const REGIONS = ["sheath/SW", "dayside/inner", "PS", "PSBL", "lobe"]
const CPAL = (; palette = [[r => c for (r, c) in zip(REGIONS, Makie.wong_colors())]; "unknown" => :gray])
const CSCALE = scales(Color = CPAL)

"(table, column dictionary, pipeline version) from data/tds_tha.parquet."
function load_table(path = TABLE)
    ds = Parquet2.Dataset(path)
    m = Parquet2.metadata(ds)
    DataFrame(ds; copycols = true), JSON.parse(m["columns"]), get(m, "pipeline", "unknown")
end

fin(x) = !ismissing(x) && x isa Real && isfinite(x)
good(v) = filter(fin, v)
posgood(v) = filter(x -> fin(x) && x > 0, v)
med(v) = isempty(v) ? NaN : median(v)
frac(f, v) = isempty(v) ? NaN : mean(f, v)

jsonsafe(x::AbstractFloat) = isfinite(x) ? round(x; sigdigits = 4) : nothing
jsonsafe(x::AbstractDict) = Dict(string(k) => jsonsafe(v) for (k, v) in x)
jsonsafe(x) = x
write_summary(path, S) = write(path, JSON.json(jsonsafe(S); pretty = true))
