using CDFDatasets, DataFrames, Dates, Statistics, LinearAlgebra, TimeseriesUtilities

const CDFDIR = joinpath(@__DIR__, "..", "data", "cdf")
const BG_HALF = 300.0   # s; background-field window is event center ± 5 min
const SPIN_HALF = 1.5   # s; a spin-resolution moment sample represents its center ± half a spin
const CTX_HALF = 60.0   # s; before/after span for the dipolarization ΔBz
const LOCAL_HALF = 30.0 # s; sheath fields rotate within a minute, so V∥ there needs a short window (scripts/bgwindow.jl)

function dayfile(prod, date)
    dir = joinpath(CDFDIR, prod)
    pre = "tha_l$(prod == "efw" ? 1 : 2)_$(prod)_$(Dates.format(date, "yyyymmdd"))_v"
    i = findlast(f -> startswith(f, pre) && endswith(f, ".cdf"), readdir(dir))
    isnothing(i) ? nothing : joinpath(dir, readdir(dir)[i])
end

"Index range of sorted `t` within [a, b]."
window(t, a, b) = searchsortedfirst(t, a):searchsortedlast(t, b)

struct Series
    t::Vector{Float64}
    x::Matrix{Float64}   # components × time
end
Series() = Series(Float64[], zeros(0, 0))

readseries(ds, tvar, var) = haskey(ds, var) ? Series(Array(ds[tvar]), ds[var]) : Series()

# Series of a variable may share `t`, which is never mutated.
function Series(t::Vector{Float64}, var::AbstractArray)
    x = convert(Array{Float64}, Array(var))
    x = x isa AbstractVector ? reshape(x, 1, :) : x
    isempty(t) && return Series()
    data, ts = tsort(Timed(x, t; dim = 2))
    Series(ts, data)
end

function wmean(s::Series, a, b; stat = tmean)
    isempty(s.t) && return (fill(NaN, max(size(s.x, 1), 1)), 0)
    r = window(s.t, a, b)
    n = count(j -> isfinite(s.x[1, j]), r)
    n == 0 ? (fill(NaN, size(s.x, 1)), 0) : (stat(view(s.x, :, r); dim = 2), n)
end
# L2 fgl/fgs contain isolated multi-kilonT spikes (>10% of spin fits on some days); medians ignore them.
wmed(s::Series, a, b) = wmean(s, a, b; stat = tmedian)

"Signed value of component `i` with the largest magnitude in [a, b], after a 3-sample running median
(single-spin velocity spikes of several 1000 km/s occur in reduced-mode moments)."
function peak(s::Series, i, a, b)
    isempty(s.t) && return NaN
    v = filter(isfinite, view(s.x, i, window(s.t, a, b)))
    length(v) >= 3 && (v = [median(view(v, j-1:j+1)) for j in 2:length(v)-1])
    isempty(v) ? NaN : v[argmax(abs.(v))]
end

wor(s::Series, a, b) = (r = window(s.t, a, b); isempty(r) ? -1 : reduce(|, Int.(s.x[1, r])))

"""Burst moments are used when present and their density is within ×3 of reduced mode (when available):
some bursts carry SST-only sums (n ~ 1e-4 cm⁻³, T ~ 25 keV)."""
function burst_ok(bb, rf, a, b)
    (nb, k) = wmean(bb.n, a, b)
    k == 0 && return false
    (nr, kr) = wmean(rf.n, a, b)
    kr == 0 || 1 / 3 < nb[1] / nr[1] < 3
end

const MOM_VARS = (n = "density", T = "avgtemp", V = "velocity_gsm", T3 = "magt3", flag = "data_quality")

function moments(ds, kind)
    names = map(v -> "tha_$(kind)_$v", MOM_VARS)
    t = any(n -> haskey(ds, n), names) ? Array(ds["tha_$(kind)_time"]) : Float64[]
    map(n -> haskey(ds, n) ? Series(t, ds[n]) : Series(), names)
end

"Per-event THEMIS-A plasma and field parameters for all events on one day."
function day_params(date, evs)
    rows = [Dict{Symbol, Any}(:id => id) for id in evs.id]
    fgmf, gmf = dayfile("fgm", date), dayfile("gmom", date)
    (isnothing(fgmf) || isnothing(gmf)) && return rows
    fgm = CDFDataset(fgmf)
    fgl = readseries(fgm, "tha_fgl_time", "tha_fgl_gsm")
    fgs = readseries(fgm, "tha_fgs_time", "tha_fgs_gsm")
    gm = CDFDataset(gmf)
    # Burst (bb) moments cover ~98% of events at spin cadence with full angular resolution; reduced (rf) mode
    # is continuous but has too few electron look directions for T⊥/T∥ (≈1 always). Full (ff) mode at event
    # times is sparse and interpolated. Context windows (±5 min) use rf.
    ion, ele = moments(gm, "ptirf"), moments(gm, "pterf")
    ionb, eleb = moments(gm, "ptibb"), moments(gm, "ptebb")

    for (row, t0, t1) in zip(rows, datetime2unix.(evs.start_time), datetime2unix.(evs.end_time))
        tc = (t0 + t1) / 2
        # Event-window B: fgl when present, else spin-res fgs
        B, nB = wmed(fgl, t0, t1)
        nB == 0 && ((B, nB) = wmed(fgs, t0 - SPIN_HALF, t1 + SPIN_HALF))
        Bbg, nbg = wmed(fgs, tc - BG_HALF, tc + BG_HALF)
        nbg == 0 && ((Bbg, nbg) = wmed(fgl, tc - BG_HALF, tc + BG_HALF))
        a, b = t0 - SPIN_HALF, t1 + SPIN_HALF
        isrc = burst_ok(ionb, ion, a, b) ? ionb : ion
        esrc = burst_ok(eleb, ele, a, b) ? eleb : ele
        ni, nni = wmean(isrc.n, a, b)
        Ti, = wmean(isrc.T, a, b)
        V, = wmean(isrc.V, a, b)
        Ti3, = wmean(isrc.T3, a, b)
        ne, nne = wmean(esrc.n, a, b)
        Te, = wmean(esrc.T, a, b)
        Te3, = wmean(esrc.T3, a, b)
        Bz_pre, = wmed(fgs, tc - CTX_HALF, tc)
        Bz_post, = wmed(fgs, tc, tc + CTX_HALF)
        Vx_peak = peak(ion.V, 1, tc - BG_HALF, tc + BG_HALF)
        bhat = Bbg / norm(Bbg)
        Vpar = dot(V, bhat)
        Vperp = V - Vpar * bhat
        Bloc, = wmed(fgs, tc - LOCAL_HALF, tc + LOCAL_HALF)
        bloc = Bloc / norm(Bloc)
        Bpre, = wmed(fgs, tc - BG_HALF, tc); Bpost, = wmed(fgs, tc, tc + BG_HALF)
        merge!(row, Dict(
            :B_x => B[1], :B_y => B[2], :B_z => B[3], :B_t => norm(B), :n_fgm => nB,
            :Bbg_x => Bbg[1], :Bbg_y => Bbg[2], :Bbg_z => Bbg[3], :Bbg_t => norm(Bbg), :n_fgm_bg => nbg,
            :Ni => ni[1], :Ti => Ti[1], :Ti_perp => (Ti3[1] + Ti3[2]) / 2, :Ti_par => Ti3[3],
            :Vi_x => V[1], :Vi_y => V[2], :Vi_z => V[3], :n_ion => nni,
            :Ne => ne[1], :Te => Te[1], :Te_perp => (Te3[1] + Te3[2]) / 2, :Te_par => Te3[3], :n_ele => nne,
            :V_par_5min => Vpar, :V_perp_5min => norm(Vperp),
            :V_par_30s => dot(V, bloc), :V_perp_30s => norm(V - dot(V, bloc) * bloc),
            :bg_rot => acosd(clamp(dot(Bpre, Bpost) / (norm(Bpre) * norm(Bpost)), -1, 1)),
            :Vperp_5min_x => Vperp[1], :Vperp_5min_y => Vperp[2], :Vperp_5min_z => Vperp[3],
            :flag_ion => wor(isrc.flag, a, b), :flag_ele => wor(esrc.flag, a, b),
            :mom_burst => (isrc === ionb) + 2(esrc === eleb),
            :dBz_1min => Bz_post[3] - Bz_pre[3], :Vx_peak_5min => Vx_peak,
        ))
    end
    rows
end

function themis_params(ev)
    days = collect(groupby(transform(ev, :start_time => ByRow(Date) => :date), :date))
    out = Vector{Vector{Dict{Symbol, Any}}}(undef, length(days))
    Threads.@threads for i in eachindex(days)
        g = days[i]
        out[i] = try
            day_params(g.date[1], g)
        catch e
            @warn "day failed" g.date[1] exception = e
            [Dict{Symbol, Any}(:id => id) for id in g.id]
        end
    end
    rows = reduce(vcat, out)
    cols = union(keys.(rows)...)
    concrete(v) = all(x -> x isa Integer, v) ? Int.(v) : float.(v)   # counts/flags are Int unless a failed day left NaN
    DataFrame([c => concrete([get(r, c, c == :id ? 0 : NaN) for r in rows]) for c in cols])
end
