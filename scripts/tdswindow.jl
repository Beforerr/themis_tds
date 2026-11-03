# TDS windows inside each Lang interval, from the automatic detector (detect/; run after a first build.jl).
using CSV, DataFrames, Dates

const DETECT_OUT = joinpath(@__DIR__, "..", "detect", "out")

"""Per-event `det_status` and `tds_*` columns. Lang's [start, end] is a burst interval with 1 s floors, so it is
extended by 1 s, as in detect/detect.jl. `nothing` if the detector has not been run."""
function tds_windows(ev)
    fb, ft = joinpath(DETECT_OUT, "bursts.csv"), joinpath(DETECT_OUT, "detections.csv")
    isfile(fb) && isfile(ft) || (@warn "detector output missing; tds_* columns omitted (run detect/run.jl, then build.jl again)"; return nothing)
    B = CSV.read(fb, DataFrame; select = [:date, :burst, :t0, :t1, :hasb])
    T = CSV.read(ft, DataFrame; select = [:date, :burst, :start_time, :end_time, :n_tds, :peak_tds, :peak_tds_Epar_ratio, :clip, :class],
        types = Dict(:start_time => DateTime, :end_time => DateTime))
    filter!(r -> startswith(r.class, "tds"), T)
    Bd, Td = groupby(B, :date), groupby(T, :date)
    rows = map(eachrow(ev)) do e
        a, b = e.start_time, e.end_time + Second(1)
        d = (date = Date(a),)
        bs = haskey(Bd, d) ? filter(r -> r.t0 <= datetime2unix(b) && r.t1 >= datetime2unix(a), Bd[d]) : B[1:0, :]
        tr = haskey(Td, d) ? filter(r -> r.start_time <= b && r.end_time >= a, Td[d]) : T[1:0, :]
        status = !isempty(tr) ? "tds" : isempty(bs) ? "no_burst" : all(bs.hasb) ? "no_tds" : "no_b"
        isempty(tr) && return (id = e.id, det_status = status, tds_start = missing, tds_end = missing, tds_n_trains = 0,
            tds_n_pulses = 0, tds_class = missing, tds_peak = NaN, tds_Epar_ratio = NaN, tds_clip = NaN)
        k = argmax(abs.(tr.peak_tds))
        (id = e.id, det_status = status, tds_start = minimum(tr.start_time), tds_end = maximum(tr.end_time),
         tds_n_trains = nrow(tr), tds_n_pulses = sum(tr.n_tds), tds_class = String(tr.class[argmax(tr.n_tds)]),
         tds_peak = tr.peak_tds[k], tds_Epar_ratio = tr.peak_tds_Epar_ratio[k], tds_clip = maximum(tr.clip))
    end
    DataFrame(rows)
end
