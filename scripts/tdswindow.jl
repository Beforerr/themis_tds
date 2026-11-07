# TDS pulses, windows and summaries per Lang interval, from the automatic detector (detect/; run after a first build.jl).
using CSV, DataFrames, Dates

const DETECT_OUT = joinpath(@__DIR__, "..", "detect", "out")
const TDS_PAD = 0.01   # s; TDS window = first to last owned pulse ± this (the detector's ±isow), so edge pulses lie inside

"""
Each event's status and its TDS pulses (TDS-like pulses of TDS trains), each pulse owned by exactly one event. Lang's
[start, end] is a burst interval with 1 s floors, so it is extended by 1 s, as in detect/detect.jl. Adjacent intervals
are often consecutive pieces of one burst segment and a train can overlap several: a pulse goes to the event whose
extended interval contains it, else the nearest one whose interval its train overlaps. `det_status` is tds when the
event owns a pulse. `nothing` if the detector has not been run.
"""
function event_pulses(ev)
    fb, ft, fp = joinpath.(DETECT_OUT, ["bursts.csv", "detections.csv", "pulses.csv"])
    all(isfile, (fb, ft, fp)) || (@warn "detector output missing; tds_* columns omitted (run detect/run.jl, then build.jl again)"; return nothing)
    B = CSV.read(fb, DataFrame; select = [:date, :burst, :t0, :t1, :hasb])
    T = CSV.read(ft, DataFrame; select = [:date, :burst, :train, :start_time, :end_time, :clip, :class],
        types = Dict(:start_time => DateTime, :end_time => DateTime))
    filter!(r -> startswith(r.class, "tds"), T)
    iv = Dict(e.id => (datetime2unix(e.start_time), datetime2unix(e.end_time) + 1) for e in eachrow(ev))
    Bd, Td = groupby(B, :date), groupby(T, :date)
    cand = DataFrame(id = Int[], date = Date[], burst = Int[], train = Int[], train_class = String[], train_clip = Float64[])
    hasburst = Dict{Int, Bool}(); allb = Dict{Int, Bool}()
    for e in eachrow(ev)
        a, b = iv[e.id]; d = (date = Date(e.start_time),)
        bs = haskey(Bd, d) ? filter(r -> r.t0 <= b && r.t1 >= a, Bd[d]) : B[1:0, :]
        hasburst[e.id] = !isempty(bs); allb[e.id] = all(bs.hasb)
        haskey(Td, d) || continue
        for r in eachrow(Td[d])
            datetime2unix(r.start_time) <= b && datetime2unix(r.end_time) >= a &&
                push!(cand, (e.id, r.date, r.burst, r.train, String(r.class), r.clip))
        end
    end
    P = CSV.read(fp, DataFrame; select = [:date, :burst, :train, :t, :chan, :a, :w_ms, :fit_s_ms, :shape, :pol, :phi,
        :snr, :iso, :ctr, :kurt, :fit_R, :tdslike])
    P = innerjoin(cand, select!(filter!(:tdslike => identity, P), Not(:tdslike)); on = [:date, :burst, :train])
    P.dist = [max(iv[i][1] - t, t - iv[i][2], 0.0) for (i, t) in zip(P.id, P.t)]
    sort!(P, [:t, :chan, :dist, :id])
    unique!(P, [:date, :burst, :t, :chan])
    owned = Set(P.id)
    status = DataFrame(id = ev.id, det_status = [i in owned ? "tds" : !hasburst[i] ? "no_burst" : allb[i] ? "no_tds" : "no_b" for i in ev.id])
    (; status, pulses = sort!(P, [:id, :t]))
end

"Per-event `det_status` and `tds_*` columns from the event's own pulses."
function tds_windows(ev, ep = event_pulses(ev))
    isnothing(ep) && return nothing
    S = combine(groupby(ep.pulses, :id)) do g
        k = argmax(abs.(g.a))
        n = combine(groupby(g, [:burst, :train]), nrow => :n, :train_class => first => :class)
        (tds_start = unix2datetime(minimum(g.t) - TDS_PAD), tds_end = unix2datetime(maximum(g.t) + TDS_PAD),
         tds_n_trains = nrow(n), tds_n_pulses = nrow(g), tds_class = n.class[argmax(n.n)],
         tds_peak = g.a[k], tds_Epar_ratio = g.pol[k], tds_clip = maximum(g.train_clip))
    end
    S = leftjoin(ep.status, S; on = :id)
    S.tds_n_trains = coalesce.(S.tds_n_trains, 0); S.tds_n_pulses = coalesce.(S.tds_n_pulses, 0)
    for c in (:tds_peak, :tds_Epar_ratio, :tds_clip); S[!, c] = coalesce.(S[!, c], NaN); end
    sort!(S, :id)
end

"The pulse table (data/tds_tha_pulses.parquet)."
function tds_pulses(ep)
    isnothing(ep) && return nothing
    P = copy(ep.pulses)
    P.time = unix2datetime.(P.t)
    P.chan = String.(P.chan); P.shape = String.(P.shape)
    select!(P, :id, :time, :t => :t_unix, :burst, :train, :train_class, :chan, :a => :amplitude, :w_ms => :width_ms,
        :fit_s_ms, :shape => :pulse_shape, :pol => :Epar_ratio, :phi, :snr, :iso, :ctr, :kurt, :fit_R)
end
