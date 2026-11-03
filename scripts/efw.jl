# THEMIS-A EFI waveform (L1 efw/eff/efp) calibration to DSL and per-event TDS statistics.
# Calibration follows SPEDAS thm_cal_efi.pro / thm_efi_despin.pro, with deviations marked below.
# It matches L2 EFI `eff_dot0_dsl` to $10^(-6)$ mV/m on a test day.
using CDFDatasets, Statistics, LinearAlgebra, DataFrames, Dates, TimeseriesUtilities, DSP

# tha row of http://themis.ssl.berkeley.edu/data/themis/tha/l1/eff/0000/tha_efi_calib_params.txt
# (single row, valid from 2002-01-01). Gains in V/count, boom lengths in m, DSC offset in mV/m (DSL).
const THA_EFI = (edc_gain = 0.000457771, eac_gain = 7.7516e-5, edc_offset = (70.0, 113.0, 0.0),
    boom = (49.6, 40.4, 5.63), shorting = (0.714286, 0.714286, 0.5), dsc_offset = (2.5, 0.0, 0.0))

# `searchsortedlast(v, x)` for queries that mostly land just past the previous result `i`."
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

"""
    efw_events(efw_path, fgm_path, windows; probe = "tha", shadow, wake, wake_halfwidth) -> DataFrame

Per (t0, t1) window (DateTime or unix s), stats of calibrated efw E (mV/m, DSL), band-passed by `bandpass!`, over
samples in [t0, t1]. The band (100–3000 Hz, as Lang's `max_Ez`) holds ms-scale structures and removes the spin tone and
slow wake/photoelectron excursions; its fixed upper edge makes 8192 and 16384 S/s bursts comparable. Background b̂ = median FGM fgl_dsl (fgs_dsl if no fgl) over the window. Columns:
`n` (efw samples; 0 → NaN stats), `ac` (fraction AC-coupled), `Bx, By, Bz` (mean DSL B, nT);
`t_Ez, Ez, Eperp_Ez, ratio_Ez, phi_Ez`: time (unix s) and signed Ez at max |Ez|, |E⊥| there, |Ez|/|E⊥|, spin phase
(deg, SSL→DSL rotation angle of `spinphase`) there; `t_Epar, Epar, Eperp_Epar, ratio_Epar, phi_Epar`: same for E∥ = E·b̂;
`rms_Ez`, `rms_Exy` (√mean(Ex²+Ey²)); `clip`: fraction of samples with a raw channel at the 16-bit ADC
rail (−32768 is also FILLVAL, so it reads as NaN; +32767 does not); maxima of events with clip > 0 are lower bounds. E∥, E⊥ and spin-plane stats need FGM spin phase (else NaN).

`shadow` (phase => half-width, deg): samples near the phases where a spin-plane boom points along the sun line
(φ ≡ 45° mod 90°) are excluded from all stats; E56 and the spin-plane pairs show spikes there
(up to ~30 mV/m) that otherwise set the event maxima. From ~2016 Q3 the 45° and 225° spikes split to c ± 4.5–6°
(also found by detect/), hence the wider mask there. `wake` phases (deg, ± `wake_halfwidth`) are a
second artifact pair found empirically from the clustering of max-|E∥| phases over all events (~180° apart).
`n` counts all samples. Residual clustering is flagged in build.jl (`efw_phase_flag`).
"""
function efw_events(efw_path, fgm_path, windows; probe = "tha", shadow = (45 => 7.0, 135 => 4.0, 225 => 7.0, 315 => 4.0),
        wake = (165.0, 346.0), wake_halfwidth = 4.0)
    rd(ds, k) = haskey(ds, k) && length(ds[k]) > 0 ? Array(ds[k]) : nothing
    fg = CDFDataset(fgm_path)
    tl, Bl, Sl = rd(fg, "$(probe)_fgl_time"), rd(fg, "$(probe)_fgl_dsl"), rd(fg, "$(probe)_fgl_ssl")
    ts, Bs = rd(fg, "$(probe)_fgs_time"), rd(fg, "$(probe)_fgs_dsl")
    φ, Tspin = tl === nothing || Sl === nothing ? (t -> fill(NaN, length(t)), 3.0) : spinphase(tl, Sl, Bl)
    ew = CDFDataset(efw_path)
    t::Vector{Float64} = rd(ew, "$(probe)_efw_time")
    win = [(_unix(a), _unix(b)) for (a, b) in windows]
    idx = [searchsortedfirst(t, a):searchsortedlast(t, b) for (a, b) in win]
    # calibrate only the bursts touching an event; offsets need the whole 0.5 s-gap segment.
    # Reading just those records of the (lazy) waveform variable avoids loading the whole day.
    segs = filter(r -> any(i -> !isempty(intersect(i, r)), idx), find_continuous_indexranges(t, 0.5))
    before = cumsum([0; length.(segs)])
    rawv = ew["$(probe)_efw"]
    tk = Vector{Float64}(undef, before[end]); raw = Matrix{eltype(rawv)}(undef, 3, before[end])
    for (k, r) in enumerate(segs)
        o = before[k]+1:before[k+1]
        tk[o] = view(t, r)
        copyto!(view(raw, :, o), view(rawv, :, r))
    end
    # A window lies wholly inside kept segments, so its samples stay contiguous in `tk`.
    pos(x) = (k = searchsortedlast(segs, x; by = first); before[k] + x - first(segs[k]) + 1)
    E, ac, ϕ = efi_dsl(tk, raw, rd(ew, "$(probe)_efw_hed_time"), rd(ew, "$(probe)_efw_hed_ac"), φ, Tspin)::Tuple{Matrix{Float64}, Vector{Bool}, Vector{Float64}}
    bandpass!(E, tk)
    masked(x) = any(((c, h),) -> abs(mod(x - c + 180, 360) - 180) < h, shadow) || any(w -> abs(mod(x - w + 180, 360) - 180) < wake_halfwidth, wake)
    for k in eachindex(ϕ)
        masked(rad2deg(ϕ[k])) && (E[:, k] .= NaN)
    end
    medB(tb, B, a, b) = (tb === nothing ? nothing :
        (r = searchsortedfirst(tb, a):searchsortedlast(tb, b); isempty(r) ? nothing : vec(median(B[:, r], dims = 2))))  # median: FGM L2 has kilo-nT spikes
    nanrms(x) = (m = count(!isnan, x); m == 0 ? NaN : sqrt(sum(v -> isnan(v) ? 0.0 : v^2, x) / m))
    rows = map(zip(win, idx)) do ((a, b), i)
        n = length(i); j = n == 0 ? (1:0) : (pos(first(i)):pos(last(i)))
        B = something(medB(tl, Bl, a, b), medB(ts, Bs, a, b), fill(NaN, 3))
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

"""Per-event `efw_*` columns for all events (needs themis.jl for `dayfile`). Window: the detector's TDS window
(`tw`, from tdswindow.jl) where it found one, else Lang's [start, end + 1 s] (his times are 1 s floors)."""
function efw_params(ev, tw = nothing)
    ev = isnothing(tw) ? transform(ev, [:start_time, :end_time] => ByRow((a, b) -> (a, b + Second(1), "list")) => [:w0, :w1, :window]) :
        transform(leftjoin(ev, select(tw, :id, :tds_start, :tds_end); on = :id), [:start_time, :end_time, :tds_start, :tds_end] =>
            ByRow((a, b, c, d) -> ismissing(c) ? (a, b + Second(1), "list") : (c, d, "tds")) => [:w0, :w1, :window])
    days = collect(groupby(transform(ev, :start_time => ByRow(Date) => :date), :date))
    out = Vector{DataFrame}(undef, length(days))
    Threads.@threads for i in eachindex(days)
        g = days[i]
        ef, fg = dayfile("efw", g.date[1]), dayfile("fgm", g.date[1])
        out[i] = DataFrame(id = g.id, window = g.window)
        (isnothing(ef) || isnothing(fg)) && continue
        try
            e = efw_events(ef, fg, collect(zip(g.w0, g.w1)))
            out[i] = hcat(out[i], select(e, Not([:Bx, :By, :Bz])))
        catch err
            @warn "efw failed" g.date[1] exception = err
        end
    end
    df = reduce((a, b) -> vcat(a, b; cols = :union), out)
    rename!(c -> c == "id" ? c : "efw_" * c, df)
    for c in (:efw_t_Ez, :efw_t_Epar)
        df[!, c] = [ismissing(x) || isnan(x) ? missing : unix2datetime(x) for x in df[!, c]]
    end
    df
end
