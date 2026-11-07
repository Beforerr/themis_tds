# Waveform-cleaning diagnostics (after build.jl) → data/efw_stages.csv, data/efw_example.csv, plotted by qa.jl.
# Per event, the spin phase and size of max |E∥| over the efw window with neither band-pass nor phase mask ("raw"), with
# the band-pass only ("band-pass"); the final stage is the table's efw_phi_Epar. Plus one burst before/after cleaning.
using CSV, DataFrames, Dates, Statistics, LinearAlgebra
include("events.jl"); include("themis.jl"); include("efw.jl"); include("tdswindow.jl")
include("table.jl")

const DATA = joinpath(@__DIR__, "..", "data")

const STAGES = DailyProduct("efw_stages", ["scripts/efw.jl", "scripts/efwstages.jl"], [:efw, :fgm], efw_day(function (c)
    raw = efw_stats!(c; band = nothing, mask = ())   # leaves c.E calibrated-only for the band-pass stage
    (raw = select(raw, :phi_Epar, :Epar), bp = select(efw_stats!(c; mask = ()), :phi_Epar, :Epar))
end))

function stages(ev, tw)
    _, items = job = efw_job(STAGES, efw_windows(ev, tw))
    res = only(daily([job]; files = dayfile))
    stage(name, k) = reduce(vcat, map(items, res) do (_, _, a), r
        DataFrame(id = a.id, stage = name, phi = isnothing(r) ? missing : r[k].phi_Epar, Epar = isnothing(r) ? missing : r[k].Epar)
    end)
    vcat(stage("raw", :raw), stage("band-pass", :bp))
end

"""E∥ of one burst over [a, b] (unix s), calibrated only (`raw`) and band-passed (`bp`; qa.jl applies the phase mask), and
band-passed E⊥ along ê⊥1 = b̂ × ẑ (DSL spin axis), ê⊥2 = b̂ × ê⊥1."""
function example(day, a, b; probe = "tha")
    fg = CDFDataset(dayfile("fgm", day)); ew = CDFDataset(dayfile("efw", day))
    tl, Bl, Sl = Array(fg["$(probe)_fgl_time"]), Array(fg["$(probe)_fgl_dsl"]), Array(fg["$(probe)_fgl_ssl"])
    φ, Tspin = spinphase(tl, Sl, Bl)
    t = Array(ew["$(probe)_efw_time"])
    r = only(filter(s -> t[first(s)] <= a && t[last(s)] >= b, find_continuous_indexranges(t, 0.5)))
    E, _, ϕ = efi_dsl(t[r], Array(ew["$(probe)_efw"][:, r]), Array(ew["$(probe)_efw_hed_time"]), Array(ew["$(probe)_efw_hed_ac"]), φ, Tspin)
    F = bandpass!(copy(E), t[r])
    deg = mod.(rad2deg.(ϕ), 360)
    bh = normalize(vec(median(Bl[:, searchsortedfirst(tl, a):searchsortedlast(tl, b)]; dims = 2)))
    k = findall(x -> a <= x <= b, t[r])
    p1 = normalize(cross(bh, [0, 0, 1])); p2 = cross(bh, p1)
    DataFrame(t = t[r][k] .- a, phase = deg[k], raw = (bh' * E[:, k])', bp = (bh' * F[:, k])',
        perp1 = (p1' * F[:, k])', perp2 = (p2' * F[:, k])')
end

"""Example burst: a post-2016 TDS event whose raw max sits on the 45° shadow and whose final max is a pulse the detector
calls TDS-like, ≥ 15° from any shadow phase (so the demonstration doesn't hinge on the mask edges)."""
function write_example(S = CSV.read(joinpath(DATA, "efw_stages.csv"), DataFrame))
    df, _, _ = load_table()
    P = CSV.read(joinpath(DETECT_OUT, "pulses.csv"), DataFrame; select = [:t, :chan, :tdslike, :shape])
    P = sort!(P[P.chan .== "Epar", :], :t)
    tdslike(t) = (k = searchsortedfirst(P.t, t - 2e-3); k <= nrow(P) && P.t[k] <= t + 2e-3 && P.tdslike[k] && P.shape[k] == "bipolar")
    raw = select(S[S.stage .== "raw", :], :id, :phi => :phi_raw, :Epar => :Epar_raw)
    c = dropmissing(innerjoin(df, raw; on = :id), [:efw_t_Epar, :efw_phi_Epar, :phi_raw, :efw_clip])
    c = c[(c.det_status .== "tds") .& (c.start_time .>= DateTime(2016, 7)) .& (abs.(c.phi_raw .- 45) .< 12) .& (c.efw_clip .== 0) .&
          (abs.(c.efw_Epar) .> 5) .& (abs.(mod.(c.efw_phi_Epar, 90) .- 45) .> 15), :]
    c = c[tdslike.(datetime2unix.(c.efw_t_Epar)), :]
    e = c[argmax(abs.(c.Epar_raw) ./ abs.(c.efw_Epar)), :]
    @info "example" e.id e.start_time e.tds_start e.tds_end e.phi_raw e.Epar_raw e.efw_Epar e.efw_phi_Epar
    X = example(Date(e.tds_start), datetime2unix(e.tds_start), datetime2unix(e.tds_end))
    X.id .= e.id
    CSV.write(joinpath(DATA, "efw_example.csv"), X)
end

if abspath(PROGRAM_FILE) == @__FILE__()
    ev = read_events()
    S = stages(ev, tds_windows(ev))
    CSV.write(joinpath(DATA, "efw_stages.csv"), S)
    write_example(S)
end
