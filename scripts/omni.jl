using CDFDatasets, DataFrames, Dates, Statistics, TimeseriesUtilities

const OMNIDIR = joinpath(@__DIR__, "..", "data", "cdf", "omni")
const OMNI_VARS = ["BX_GSE", "BY_GSM", "BZ_GSM", "flow_speed", "proton_density", "Pressure", "AE_INDEX", "AL_INDEX", "SYM_H"]

"Minute-indexed OMNI table spanning all downloaded months (gaps = NaN)."
function load_omni()
    files = sort(filter(endswith(".cdf"), readdir(OMNIDIR; join = true)))
    parts = map(files) do f
        d = CDFDataset(f)
        DataFrame([:time => Array(d["Epoch"]); [Symbol(v) => Float64.(Array(d[v])) for v in OMNI_VARS]])
    end
    o = reduce(vcat, parts)
    sort!(o, :time)
end

"Hourly Kp from OMNI2 low-res ASCII (word 39 = 10·Kp, repeated over each 3-h interval)."
function load_kp(dir = joinpath(OMNIDIR, "..", "omni2"))
    rows = [split(l) for f in sort(readdir(dir; join = true)) for l in eachline(f)]
    DataFrame(time = [DateTime(parse(Int, r[1])) + Day(parse(Int, r[2]) - 1) + Hour(parse(Int, r[3])) for r in rows],
        Kp = [(k = parse(Int, r[39]); k == 99 ? NaN : k / 10) for r in rows])  # 99 = missing
end

function _mean_in(t, x, a, b)
    r = searchsortedfirst(t, a):searchsortedlast(t, b)
    tmean(view(x, r))
end
function _ext_in(f, t, x, a, b)
    r = searchsortedfirst(t, a):searchsortedlast(t, b)
    v = filter(isfinite, view(x, r))
    isempty(v) ? NaN : f(v)
end

"""Solar-wind and geomagnetic context per event. OMNI is already propagated to the bow-shock nose.
At-event values average the minutes overlapping the event; `_1h` columns summarize the hour before it."""
function omni_params(ev, o = load_omni(), kp = load_kp())
    t = o.time
    t0, t1 = floor.(ev.start_time, Minute), floor.(ev.end_time, Minute)
    df = DataFrame(id = ev.id)
    for v in OMNI_VARS
        df[!, Symbol(v)] = _mean_in.(Ref(t), Ref(o[!, v]), t0, t1)
    end
    df.BZ_GSM_1h = _mean_in.(Ref(t), Ref(o.BZ_GSM), t0 .- Minute(60), t0 .- Minute(1))
    df.AL_min_1h = _ext_in.(minimum, Ref(t), Ref(o.AL_INDEX), t0 .- Minute(60), t0 .- Minute(1))
    df.AE_max_1h = _ext_in.(maximum, Ref(t), Ref(o.AE_INDEX), t0 .- Minute(60), t0 .- Minute(1))
    i = searchsortedlast.(Ref(kp.time), ev.start_time)
    df.Kp = kp.Kp[i]
    df
end
