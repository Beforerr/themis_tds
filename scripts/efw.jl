# THEMIS-A EFI waveform (L1 efw/eff/efp) calibration to DSL and per-event TDS statistics.
# Calibration follows SPEDAS thm_cal_efi.pro / thm_efi_despin.pro, with deviations marked below.
# It matches L2 EFI `eff_dot0_dsl` to $10^(-6)$ mV/m on a test day.
using CDFDatasets, Statistics, LinearAlgebra, DataFrames, Dates, TimeseriesUtilities, DSP
@isdefined(ProbeDay) || include("daily.jl")

# tha row of http://themis.ssl.berkeley.edu/data/themis/tha/l1/eff/0000/tha_efi_calib_params.txt
# (single row, valid from 2002-01-01). Gains in V/count, boom lengths in m, DSC offset in mV/m (DSL).
const THA_EFI = (edc_gain = 0.000457771, eac_gain = 7.7516e-5, edc_offset = (70.0, 113.0, 0.0),
    boom = (49.6, 40.4, 5.63), shorting = (0.714286, 0.714286, 0.5), dsc_offset = (2.5, 0.0, 0.0))

# `searchsortedlast(v, x)` for queries that mostly land just past the previous result `i`.
function searchlast_from(v, x, i)
    1 <= i <= lastindex(v) && !isless(x, v[i]) || return searchsortedlast(v, x)
    for _ in 1:8
        i < lastindex(v) && !isless(x, v[i+1]) || return i
        i += 1
    end
    searchsortedlast(v, x)
end

"""
    spinphase(tfgm, Bssl, Bdsl) -> (f, Tspin)

Spin phase φ (rad, `B_dsl = Rz(φ) B_ssl`, the ssl2dsl convention) at times `t` as `f(t)`, from L2 FGM fgl/fgh, which
carry both frames: their ratio is the SPEDAS spin model (with corrections) evaluated at FGM time tags, so no
L1 state file is needed. Linear interpolation inside FGM gaps ≤ `maxgap` s, else extrapolation at
the median rate (`Tspin`) up to `maxext` s from the nearest FGM sample, else NaN.
"""
function spinphase(tf, s, d; maxgap = 1.0, maxext = 5.0)
    # L2 time tags step back at some packet boundaries and repeat on some days
    wrapped, ts = tsort(Timed([Float64(atan(d[2, i], d[1, i]) - atan(s[2, i], s[1, i])) for i in axes(d, 2)], tf))
    ψ = unwrap_periodic(wrapped, 2π)
    n = length(ts)
    Tspin = 2π / median!([(ψ[i+1] - ψ[i]) / (ts[i+1] - ts[i]) for i in 1:n-1 if 0 < ts[i+1] - ts[i] <= maxgap])
    function f(t)
        ϕ = parent(tinterp(Timed(ψ, ts), t; max_gap = maxgap))
        for k in findall(isnan, ϕ)
            x = t[k]
            i = searchsortedlast(ts, x)
            1 <= i < n && ts[i+1] - ts[i] <= maxgap && continue
            j = i < 1 || (i < n && ts[i+1] - x < x - ts[i]) ? i + 1 : i
            abs(x - ts[j]) <= maxext && (ϕ[k] = ψ[j] + 2π / Tspin * (x - ts[j]))
        end
        ϕ
    end
    f, Tspin
end

# IDL smooth(y, w, /nan, /edge_truncate): centered boxcar skipping NaNs, edge values replicated.
# `out`, `c` and `m` are scratch reused across segments.
function boxcar_edge!(out, c, m, y, w)
    n = length(y); h = w ÷ 2
    resize!(out, n); resize!(c, n + 1); resize!(m, n + 1)
    c[1] = 0.0; m[1] = 0
    for i in 1:n
        nan = isnan(y[i])
        c[i+1] = c[i] + (nan ? 0.0 : Float64(y[i])); m[i+1] = m[i] + !nan
    end
    y1, yn = y[1], y[n]
    for i in 1:n
        a, b = i - h, i + h
        lo, hi = clamp(a - 1, 0, n) + 1, clamp(b, 0, n) + 1
        pa, pb = isnan(y1) ? 0 : max(1 - a, 0), isnan(yn) ? 0 : max(b - n, 0)
        out[i] = (c[hi] - c[lo] + (pa > 0 ? pa * y1 : 0.0) + (pb > 0 ? pb * yn : 0.0)) / (m[hi] - m[lo] + pa + pb)
    end
    out
end

"""
    efi_dsl(t, raw, hed_t, hed_ac, φ, Tspin; cal = THA_EFI, nspins = 20) -> (E, ac, ϕ)

`E`: 3×n Matrix{Float64} (mV/m, DSL); `ac`: per-sample AC-coupling flag; `ϕ`: per-sample spin phase (rad).

Raw L1 EFI E12/E34/E56 counts → DSL. Coupling per sample from the nearest header's `hed_ac` (EAC vs
EDC gain). Per segment (split at 0.5 s gaps and coupling changes), DC-coupled data follow
thm_cal_efi: E12/E34 offset = boxcar over `nspins` spins, truncated to whole spins inside short
segments, NaN if < 1 spin; E56 offset from the cal file; DSL `dsc_offset` subtracted.
AC-coupled data deviate from SPEDAS: the EAC high-pass already removes DC and spin tone, so the
offset of all three channels is the segment mean (any length) and `dsc_offset` is not applied
(SPEDAS keeps the E56 ADC offset and adds a spurious -2.5 mV/m to AC Ex).
"""
function efi_dsl(t, raw, hed_t, hed_ac, φ, Tspin; cal = THA_EFI, nspins = 20)
    n = length(t); E = fill(NaN, 3, n)
    n == 0 && return E, Bool[], Float64[]
    ac = Vector{Bool}(undef, n); h = 0
    for (k, x) in enumerate(t)
        h = searchlast_from(hed_t, x, h)
        j = clamp(h, 1, length(hed_t))
        j < length(hed_t) && hed_t[j+1] - x < x - hed_t[j] && (j += 1)
        ac[k] = hed_ac[j] != 0
    end
    ϕ = φ(t)::Vector{Float64}
    b = [0; findall(i -> t[i+1] - t[i] > 0.5 || ac[i+1] != ac[i], 1:n-1); n]
    leff = cal.boom .* cal.shorting
    y, smooth, cum, cnt = Float64[], Float64[], Float64[], Int[]
    for k in 1:length(b)-1
        r = b[k]+1:b[k+1]; isac = ac[r[1]]
        g = isac ? cal.eac_gain : cal.edc_gain
        nps = length(r) > 1 ? abs(Tspin) * (length(r) - 1) / (t[r[end]] - t[r[1]]) : Inf  # samples per spin
        nsp = length(r) > nspins * nps ? nspins : floor(Int, (length(r) - 1) / nps)
        resize!(y, length(r))
        for c in 1:3
            y .= view(raw, c, r)
            off = isac ? mean(Iterators.filter(!isnan, y)) : c == 3 ? cal.edc_offset[3] :
                  nsp >= 1 ? boxcar_edge!(smooth, cum, cnt, y, ceil(Int, nsp * nps)) : NaN
            @views E[c, r] .= (-1000g / leff[c]) .* (y .- off)
        end
        dsc = isac ? (0.0, 0.0, 0.0) : cal.dsc_offset
        for i in r
            e1, e2, e3 = E[1, i], E[2, i], E[3, i]
            sx = (e2 - e1) / √2; sy = -(e1 + e2) / √2   # sensor (spg) → SSL, from thm_efi_despin
            sφ, cφ = sincos(ϕ[i])
            E[1, i] = sx * cφ - sy * sφ - dsc[1]
            E[2, i] = sx * sφ + sy * cφ - dsc[2]
            E[3, i] = e3 - dsc[3]
        end
    end
    E, ac, ϕ
end

"""
Zero-phase Butterworth band-pass (`band` Hz, −6 dB at the edges) of each row of `E`, per contiguous chunk of `t`
(bursts carry ~30 ms packet gaps). NaN (rail/fill) samples are bridged for filtering and stay NaN; samples within
`edge` s of a chunk end, where the filter rings, become NaN.
"""
function bandpass!(E, t; band = (100.0, 3000.0), order = 4, edge = 0.03)
    for r in find_continuous_indexranges(t, 1.5 * median(diff(t)))
        length(r) < 64 && (E[:, r] .= NaN; continue)
        fs = (length(r) - 1) / (t[r[end]] - t[r[1]])
        f = digitalfilter(Bandpass(band[1], min(band[2], 0.45fs)), Butterworth(order); fs)
        e = round(Int, edge * fs)
        for c in axes(E, 1)
            x = E[c, r]; bad = isnan.(x)
            all(bad) && continue
            any(bad) && bridge!(x, bad)
            y = filtfilt(f, x); y[bad] .= NaN
            y[1:min(e, end)] .= NaN; y[max(1, end - e + 1):end] .= NaN
            E[c, r] = y
        end
    end
    E
end

"Linear interpolation over `bad` samples of `x` from the nearest good ones (held flat at the ends)."
function bridge!(x, bad)
    g = findall(!, bad)
    for k in findall(bad)
        j = searchsortedfirst(g, k)
        x[k] = j > length(g) ? x[g[end]] : j == 1 ? x[g[1]] :
            (a = g[j-1]; b = g[j]; x[a] + (x[b] - x[a]) * (k - a) / (b - a))
    end
    x
end

_unix(t::DateTime) = datetime2unix(t)
_unix(t::Real) = Float64(t)

# Spin phases (deg) whose samples are excluded. Near φ ≡ 45° mod 90° a spin-plane boom points along the sun line;
# E56 and the spin-plane pairs spike there (up to ~80 mV/m, broadband, so the band-pass keeps them) and would otherwise
# set the event maxima. From ~2016 Q3 the 45° and 225° spikes split to c ± 4.5–6° (also found by detect/), and the 45°
# one leaks to c + 10° (excess of max-|E∥| phases at 52° and 55°, post-2016 only), hence the wider, asymmetric intervals.
# scripts/efwstages.jl shows the max-|E∥| phases without band-pass and mask.
const EFW_MASK = ((38.0, 57.0), (131.0, 139.0), (218.0, 232.0), (311.0, 319.0))

"""
    efw_events(d::ProbeDay, windows; mask = EFW_MASK, band = (100.0, 3000.0)) -> DataFrame

Per (t0, t1) window (DateTime or unix s), stats of calibrated efw E (mV/m, DSL), band-passed by `bandpass!`, over
samples in [t0, t1]. The band (100–3000 Hz, as Lang's `max_Ez`) holds ms-scale structures and removes the spin tone and
slow wake/photoelectron excursions; its fixed upper edge makes 8192 and 16384 S/s bursts comparable. Background b̂ = median FGM fgl_dsl (fgs_dsl if no fgl) over the window. Columns:
`n` (efw samples; 0 → NaN stats), `ac` (fraction AC-coupled), `Bx, By, Bz` (median DSL B, nT);
`t_Ez, Ez, Eperp_Ez, ratio_Ez, phi_Ez`: time (unix s) and signed Ez at max |Ez|, |E⊥| there, |Ez|/|E⊥|, spin phase
(deg, SSL→DSL rotation angle of `spinphase`) there; `t_Epar, Epar, Eperp_Epar, ratio_Epar, phi_Epar`: same for E∥ = E·b̂;
`rms_Ez`, `rms_Exy` (√mean(Ex²+Ey²)); `clip`: fraction of samples with a raw channel at the 16-bit ADC
rail (−32768 is also FILLVAL, so it reads as NaN; +32767 does not); maxima of events with clip > 0 are lower bounds. E∥, E⊥ and spin-plane stats need FGM spin phase (else NaN).

`mask`: spin-phase intervals [lo, hi) (deg) excluded from all stats, see `EFW_MASK`; `band = nothing` skips the band-pass.
`n` counts all samples.
"""
efw_events(d::ProbeDay, windows; kw...) = efw_stats!(efw_calibrated(d, windows); kw...)

"Calibrated E of the bursts touching `windows`, with each window's sample range `j` into it and its b̂ field `B`."
function efw_calibrated(d::ProbeDay, windows)
    probe = d.probe
    rd(kind, k) = dayvar(d, kind, k)
    tl, Bl, Sl = rd(:fgm, "$(probe)_fgl_time"), rd(:fgm, "$(probe)_fgl_dsl"), rd(:fgm, "$(probe)_fgl_ssl")
    ts, Bs = rd(:fgm, "$(probe)_fgs_time"), rd(:fgm, "$(probe)_fgs_dsl")
    φ, Tspin = tl === nothing || Sl === nothing ? (t -> fill(NaN, length(t)), 3.0) : spinphase(tl, Sl, Bl)
    t::Vector{Float64} = rd(:efw, "$(probe)_efw_time")
    win = [(_unix(a), _unix(b)) for (a, b) in windows]
    idx = [searchsortedfirst(t, a):searchsortedlast(t, b) for (a, b) in win]
    # calibrate only the bursts touching an event; offsets need the whole 0.5 s-gap segment.
    # Reading just those records of the (lazy) waveform variable avoids loading the whole day.
    segs = filter(r -> any(i -> !isempty(intersect(i, r)), idx), find_continuous_indexranges(t, 0.5))
    before = cumsum([0; length.(segs)])
    rawv = dayset(d, :efw)["$(probe)_efw"]
    tk = Vector{Float64}(undef, before[end]); raw = Matrix{eltype(rawv)}(undef, 3, before[end])
    for (k, r) in enumerate(segs)
        o = before[k]+1:before[k+1]
        tk[o] = view(t, r)
        copyto!(view(raw, :, o), view(rawv, :, r))
    end
    # A window lies wholly inside kept segments, so its samples stay contiguous in `tk`.
    pos(x) = (k = searchsortedlast(segs, x; by = first); before[k] + x - first(segs[k]) + 1)
    j = [isempty(i) ? (1:0) : (pos(first(i)):pos(last(i))) for i in idx]
    E, ac, ϕ = efi_dsl(tk, raw, rd(:efw, "$(probe)_efw_hed_time"), rd(:efw, "$(probe)_efw_hed_ac"), φ, Tspin)::Tuple{Matrix{Float64}, Vector{Bool}, Vector{Float64}}
    medB(tb, B, a, b) = (tb === nothing ? nothing :
        (r = searchsortedfirst(tb, a):searchsortedlast(tb, b); isempty(r) ? nothing : vec(median(B[:, r], dims = 2))))  # median: FGM L2 has kilo-nT spikes
    B = [something(medB(tl, Bl, a, b), medB(ts, Bs, a, b), fill(NaN, 3)) for (a, b) in win]
    (; tk, raw, E, ac, ϕ, j, B)
end

"`efw_events` columns from `efw_calibrated` output `c`. Band-pass and mask overwrite `c.E`; with neither, it is left as is."
function efw_stats!(c; mask = EFW_MASK, band = (100.0, 3000.0))
    (; tk, raw, E, ac, ϕ) = c
    isnothing(band) || bandpass!(E, tk; band)
    if !isempty(mask)
        for k in eachindex(ϕ)
            x = mod(rad2deg(ϕ[k]), 360)
            any(((lo, hi),) -> lo <= x < hi, mask) && (E[:, k] .= NaN)
        end
    end
    nanrms(x) = (m = count(!isnan, x); m == 0 ? NaN : sqrt(sum(v -> isnan(v) ? 0.0 : v^2, x) / m))
    rows = map(c.j, c.B) do j, B
        n = length(j)
        bh = B ./ norm(B)
        ez = view(E, 3, j)
        par = [bh[1] * E[1, k] + bh[2] * E[2, k] + bh[3] * E[3, k] for k in j]
        perp = [sqrt(max(E[1, k]^2 + E[2, k]^2 + E[3, k]^2 - p^2, 0)) for (k, p) in zip(j, par)]
        function pick(x)
            m = isempty(x) ? (-Inf, 0) : findmax(v -> isnan(v) ? -Inf : abs(v), x)
            m[1] == -Inf && return (NaN, NaN, NaN, NaN, NaN)
            k = m[2]; (tk[j[k]], x[k], perp[k], m[1] / perp[k], mod(rad2deg(ϕ[j[k]]), 360))
        end
        tz, z, pz, rz, fz = pick(ez); tp, p, pp, rp, fp = pick(par)
        (n = n, ac = n == 0 ? NaN : mean(view(ac, j)), Bx = B[1], By = B[2], Bz = B[3],
         t_Ez = tz, Ez = z, Eperp_Ez = pz, ratio_Ez = rz, phi_Ez = fz,
         t_Epar = tp, Epar = p, Eperp_Epar = pp, ratio_Epar = rp, phi_Epar = fp,
         rms_Ez = nanrms(ez), rms_Exy = n == 0 ? NaN : hypot(nanrms(view(E, 1, j)), nanrms(view(E, 2, j))),
         clip = n == 0 ? NaN : count(k -> any(v -> !(abs(v) < 32767), view(raw, :, k)), j) / n)
    end
    DataFrame(rows)
end

"""Per-event efw windows `w0`, `w1` and their kind `window`: the detector's TDS window (`tw`, from tdswindow.jl) where it
found one, else Lang's [start, end + 1 s] (his times are 1 s floors)."""
efw_windows(ev, tw = nothing) = isnothing(tw) ?
    transform(ev, [:start_time, :end_time] => ByRow((a, b) -> (a, b + Second(1), "list")) => [:w0, :w1, :window]) :
    transform(leftjoin(ev, select(tw, :id, :tds_start, :tds_end); on = :id), [:start_time, :end_time, :tds_start, :tds_end] =>
        ByRow((a, b, c, d) -> ismissing(c) ? (a, b + Second(1), "list") : (c, d, "tds")) => [:w0, :w1, :window])

"Day function of efw products: `f(efw_calibrated(d, windows))`, `nothing` without efw or FGM files."
efw_day(f) = (d, w) -> isnothing(input(d, :efw)) || isnothing(input(d, :fgm)) ? nothing : f(efw_calibrated(d, collect(zip(w.w0, w.w1))))

const EFW = DailyProduct("efw", ["scripts/efw.jl"], [:efw, :fgm], efw_day(efw_stats!))

efw_job(p, W) = perday(p, select(W, :id, :w0, :w1), Date.(W.start_time))

"Per-event `efw_*` columns from `daily` results `res` of `efw_job(EFW, W)`."
function efw_table(W, (_, items), res)
    kind = Dict(zip(W.id, W.window))
    out = map(items, res) do (_, _, a), e
        df = DataFrame(id = a.id, window = [kind[i] for i in a.id])
        isnothing(e) ? df : hcat(df, select(e, Not([:Bx, :By, :Bz])))
    end
    df = reduce((a, b) -> vcat(a, b; cols = :union), out)
    rename!(c -> c == "id" ? c : "efw_" * c, df)
    for c in (:efw_t_Ez, :efw_t_Epar)
        df[!, c] = [ismissing(x) || isnan(x) ? missing : unix2datetime(x) for x in df[!, c]]
    end
    df
end
