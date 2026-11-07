# Automatic TDS detection in THEMIS-A EFW burst waveforms. See detect/report.typ.
include(joinpath(@__DIR__, "..", "scripts", "dayfile.jl"))
include(joinpath(@__DIR__, "..", "scripts", "efw.jl"))
include(joinpath(@__DIR__, "..", "scripts", "events.jl"))
using CSV, StatsBase, TimeseriesUtilities

struct Burst
    t::Vector{Float64}      # unix s
    E::Matrix{Float64}      # 3×n DSL mV/m
    bh::Matrix{Float64}     # 3×n DSL unit B (NaN if no FGM)
    φ::Vector{Float64}      # spin phase, deg in [0, 360)
    bt::Vector{Float64}     # |B| nT
    fs::Float64
    rail::BitVector         # any raw channel at the ADC rail, as efw_events' `clip`
end

"All calibrated efw bursts (segments split at 0.5 s gaps) of one probe-day."
load_day(date::Date; probe = "tha") = load_day(ProbeDay(probe, date, dayfile))
function load_day(d::ProbeDay)
    probe = d.probe
    rd(kind, k) = dayvar(d, kind, k)
    (isnothing(input(d, :efw)) || isnothing(input(d, :fgm))) && return Burst[]
    tl, Bl, Sl = rd(:fgm, "$(probe)_fgl_time"), rd(:fgm, "$(probe)_fgl_dsl"), rd(:fgm, "$(probe)_fgl_ssl")
    (tl === nothing || Sl === nothing) && return Burst[]
    φf, Tspin = spinphase(tl, Sl, Bl)
    th, Bh = rd(:fgm, "$(probe)_fgh_time"), rd(:fgm, "$(probe)_fgh_dsl")
    th === nothing && ((th, Bh) = (tl, Bl))
    t, raw = rd(:efw, "$(probe)_efw_time"), rd(:efw, "$(probe)_efw")
    t === nothing && return Burst[]
    E, _, ϕ = efi_dsl(t, raw, rd(:efw, "$(probe)_efw_hed_time"), rd(:efw, "$(probe)_efw_hed_ac"), φf, Tspin)
    bursts = Burst[]
    for r in find_continuous_indexranges(t, 0.5)
        length(r) < 4096 && continue
        tr = t[r]
        B = parent(tinterp(tsort(Timed(Float64.(Bh), th; dim = 2)), tr; max_gap = 1.0))
        fs = round((length(r) - 1) / (tr[end] - tr[1]) / 1024) * 1024
        bt = vec(sqrt.(sum(abs2, B, dims = 1)))
        rail = BitVector([any(v -> !(abs(v) < 32767), view(raw, :, k)) for k in r])
        push!(bursts, Burst(tr, E[:, r], B ./ bt', mod.(rad2deg.(view(ϕ, r)), 360), bt, fs, rail))
    end
    bursts
end

"Centered boxcar mean of odd width `w` (edges: shrinking window)."
function boxmean(x, w)
    n = length(x); h = w ÷ 2; c = [0.0; cumsum(x)]
    [(c[min(i + h, n)+1] - c[max(i - h, 1)]) / (min(i + h, n) - max(i - h, 1) + 1) for i in 1:n]
end

"High-pass: `x` minus a twice-applied boxcar of `τ` s (triangular kernel, cutoff ≈ 0.6/τ Hz)."
highpass(x, fs, τ) = (w = 2round(Int, τ * fs / 2) + 1; x .- boxmean(boxmean(x, w), w))

"""
`highpass` per contiguous chunk of `t` (bursts carry ~30 ms packet gaps). Returns the filtered signal
and a mask of samples within `τ` of a chunk edge, where the filter rings.
"""
function highpass_chunks(x, t, fs, τ)
    y = similar(x, Float64); edge = falses(length(x)); e = round(Int, τ * fs)
    for r in find_continuous_indexranges(t, 1.5 / fs)
        y[r] .= highpass(view(x, r), fs, τ)
        edge[r[1]:min(r[end], r[1] + e)] .= true; edge[max(r[1], r[end] - e):r[end]] .= true
    end
    y, edge
end


# Raw count -32768 is both the CDF fill value and negative ADC full scale, so CDFDatasets turns
# clipped samples into NaN (+32767 stays finite; `Burst.rail` catches both).
"Linearly interpolate NaN runs in place; return the NaN mask."
function fillnan!(x)
    bad = isnan.(x); all(bad) && (fill!(x, 0.0); return bad)
    g = findall(!, bad)
    for i in findall(bad)
        k = searchsortedfirst(g, i)
        x[i] = k == 1 ? x[g[1]] : k > length(g) ? x[g[end]] :
               x[g[k-1]] + (x[g[k]] - x[g[k-1]]) * (i - g[k-1]) / (g[k] - g[k-1])
    end
    bad
end

const PARAMS = (τ = 0.004,          # s, high-pass time scale (cutoff ≈ 150 Hz)
    qblk = 0.005, qq = 0.2,         # quiet floor σq: qq-quantile of qblk-s block rms
    snr = 6.0,                      # candidate peak ≥ snr σq
    amin = 1.0,                     # mV/m floor on the filtered peak
    wmin = 2, wmax = 2.5e-3,        # half-max lobe width: ≥ 2 samples, ≤ 2.5 ms
    iso = 3.0, isow = 0.010,        # peak / rms in ±isow s excluding the pulse
    nms = 0.5e-3,                   # s, non-maximum suppression half-window
    satw = 0.02,                    # s, spikes this close to a clipped sample are flagged `sat`
    gap = 0.25,                     # s, spikes closer than this form one train
    fitw = 0.003,                   # s, minimum shape-fit half-window
    kw = 0.025,                     # s, half-window for the pulse's local kurtosis
    # TDS-like pulse: artifact-free, |E∥|/|E⊥| ≥ pol and ≥ votes of: local kurtosis ≥ kurt (intermittent),
    # contrast ≥ ctr (isolated), shape-fit residual ≤ rfit (coherent single structure)
    pol = 1.0, kurt = 8.0, ctr = 6.0, rfit = 0.35, votes = 2,
    art = 0.5,                      # artifact train: fraction of spin-recurrent / shadow pulses ≥
    acmax = 0.5,                    # coherent wave train: median periodicity ac ≥ acmax
    shapefrac = 2 / 3,              # tds_bipolar / tds_unipolar if ≥ this fraction of TDS-like pulses
    nmin = 3)                       # TDS-like pulses per TDS train

"Quiet-gap noise floor: `p.qq` quantile of the rms in `p.qblk` blocks free of chunk edges."
function quiet_sigma(x, edge, fs; p = PARAMS)
    m = max(8, round(Int, p.qblk * fs))
    r = [sqrt(sum(abs2, view(x, k:k+m-1)) / m) for k in 1:m:length(x)-m+1 if !any(view(edge, k:k+m-1))]
    isempty(r) ? NaN : quantile(r, p.qq)
end

"""
Candidate spikes in filtered channel `x` with noise floor `σ`. NamedTuples with sample index `i`, signed
peak `a`, `snr`, half-max lobe width `w` (samples), `opp` (largest opposite-sign excursion within
±(2w + nms) over |a|; ≳ 0.4 ⇒ bipolar), `iso`, `ctr` (local contrast: |a| over the median-based σ in ±isow
excluding the pulse; a Gaussian fluctuation's extremes reach only ≈ 3–4), `ac`/`aclag`: max normalized autocorrelation in ±isow at lags
2w … isow/2 and its lag (samples), ≈ 1 for a quasi-monochromatic packet. Samples flagged in `edge` are skipped.
"""
function spikes(x, fs, σ, edge; p = PARAMS)
    n = length(x); ax = abs.(x)
    thr = max(p.snr * σ, p.amin)
    m = max(1, round(Int, p.nms * fs)); W = round(Int, p.isow * fs)
    out = NamedTuple{(:i, :a, :snr, :w, :opp, :iso, :ctr, :ac, :aclag), Tuple{Int, Float64, Float64, Int, Float64, Float64, Float64, Float64, Int}}[]
    buf = Float64[]
    for i in 1:n
        (ax[i] >= thr && !edge[i]) || continue
        lo, hi = max(1, i - m), min(n, i + m)
        any(j -> ax[j] > ax[i] || (ax[j] == ax[i] && j < i), lo:hi) && continue
        s = sign(x[i]); h = ax[i] / 2
        L = i; while L > 1 && s * x[L-1] >= h; L -= 1; end
        R = i; while R < n && s * x[R+1] >= h; R += 1; end
        w = R - L + 1
        (w >= p.wmin && w / fs <= p.wmax) || continue
        k = 2w + m
        opp = maximum(j -> -s * x[j], max(1, i - k):min(n, i + k)) / ax[i]
        ss = 0.0; empty!(buf)
        for j in max(1, i - W):min(n, i + W)
            abs(j - i) <= k && continue
            ss += x[j]^2; push!(buf, ax[j])
        end
        iso = isempty(buf) ? Inf : ax[i] / sqrt(ss / length(buf))
        ctr = isempty(buf) ? Inf : ax[i] / (1.4826 * median!(buf))
        lo, hi = max(1, i - W), min(n, i + W); ac = 0.0; lag = 0
        for L in 2w:(hi - lo) ÷ 4
            sxy = sxx = syy = 0.0
            @inbounds @simd for j in lo:hi-L
                sxy += x[j] * x[j+L]; sxx += x[j]^2; syy += x[j+L]^2
            end
            r = sxy / sqrt(sxx * syy)
            r > ac && ((ac, lag) = (r, L))
        end
        push!(out, (i = i, a = x[i], snr = ax[i] / σ, w = w, opp = opp, iso = iso, ctr = ctr, ac = ac, aclag = lag))
    end
    out
end

const SQE = sqrt(exp(1.0))   # normalizes u·exp(-u²/2) to unit peak

"""
Least squares x ≈ c + a G + b D on samples lo:hi, with G = exp(-u²/2) (unipolar), D = √e u exp(-u²/2)
(bipolar, derivative of Gaussian), u = (j - t0)/s, both unit peak. Returns (a, b, residual ss, total ss
about the mean).
"""
function pulsefit(x, t0, s, lo, hi)
    n = hi - lo + 1; sx = sG = sD = 0.0
    @inbounds for j in lo:hi
        u = (j - t0) / s; e = exp(-u^2 / 2)
        sx += x[j]; sG += e; sD += SQE * u * e
    end
    mx, mG, mD = sx / n, sG / n, sD / n
    GG = GD = DD = Gx = Dx = xx = 0.0
    @inbounds for j in lo:hi
        u = (j - t0) / s; e = exp(-u^2 / 2)
        G = e - mG; D = SQE * u * e - mD; y = x[j] - mx
        GG += G * G; GD += G * D; DD += D * D; Gx += G * y; Dx += D * y; xx += y * y
    end
    d = GG * DD - GD^2
    d <= 0 && return (0.0, 0.0, xx, xx)
    a = (DD * Gx - GD * Dx) / d; b = (GG * Dx - GD * Gx) / d
    a, b, max(xx - (a * Gx + b * Dx), 0.0), xx
end

"""
Best Gaussian/derivative-of-Gaussian fit to the pulse at sample `i` with half-max lobe width `w` (samples):
grid over width s and centre t0, window t0 ± max(10s, `hmin` samples) so the flanks (quiet for an isolated pulse,
oscillating for a wave) count. Returns
(t0, s, lo, hi, a, b, R) with R = residual / total ss (≈ 0 for a clean isolated pulse; ragged or
oscillating waveforms leave R ≳ 0.5).
"""
function shapefit(x, i, w, hmin)
    best = (0.0, 0.0, 0, 0, 0.0, 0.0, Inf)
    for f in (0.3, 0.42, 0.55, 0.7, 0.9), δ in (-1.0, -0.5, 0.0, 0.5, 1.0)
        s = max(f * w, 0.5); t0 = i + δ * s; h = max(10s, hmin)
        lo, hi = max(1, floor(Int, t0 - h)), min(length(x), ceil(Int, t0 + h))
        hi - lo < 6 && continue
        a, b, r, tot = pulsefit(x, t0, s, lo, hi)
        R = tot > 0 ? r / tot : Inf
        R < best[7] && (best = (t0, s, lo, hi, a, b, R))
    end
    best
end

"""
Band-limited polarization at a pulse: the fitted template (t0, s) is regressed on each filtered DSL
component, so E∥ and E⊥ are compared in the pulse's own band. Returns `pol` = |E∥|/|E⊥| (3-D) and
`pair` = largest share of spin-plane template energy on one sensor boom pair (E12 or E34), ≈ 1 for sun-shadow
pulses but also for field-aligned pulses when b̂ lies along a pair; NaN without spin-plane signal. Boom response: a structure of scale L ≲ boom length is
underestimated by the 50 m spin-plane booms but not the 5.6 m axial ones; L = v × width needs the unknown
structure speed, so no correction is applied (for v ≳ 300 km/s and widths ≳ 0.3 ms, L ≳ 100 m and the
spin-plane response is ≳ 0.8).
"""
function polarization(hp, bh, φdeg, t0, s, lo, hi)
    A = [(r = pulsefit(hp[c], t0, s, lo, hi); (r[1], r[2])) for c in 1:3]
    par = (sum(bh[c] * A[c][1] for c in 1:3), sum(bh[c] * A[c][2] for c in 1:3))
    tot = sum(sum(abs2, a) for a in A)
    pol = all(isfinite, bh) ? hypot(par...) / sqrt(max(tot - sum(abs2, par), 1e-30)) : NaN
    # DSL → SSL (inverse of efi_dsl's rotation), then SSL → sensor pairs (inverse of thm_efi_despin)
    c, sn = cosd(φdeg), sind(φdeg)
    e12 = e34 = 0.0
    for k in 1:2
        x, y = A[1][k], A[2][k]
        sx, sy = x * c + y * sn, -x * sn + y * c
        e12 += ((sx + sy) / √2)^2; e34 += ((sx - sy) / √2)^2
    end
    # Spin-plane channels clipped across the fit are flat after `fillnan!`, leaving only rounding noise here.
    pol, e12 + e34 > eps() * tot ? max(e12, e34) / (e12 + e34) : NaN
end

"""
Spin-synchronous recurrence of spikes at times `t`, phases `φ` (deg) within one burst of spin period `T`
spanning `dur` s. Flagged when a spike recurs (±dφ) in at least half of the other spins the burst covers and
well above the chance count `λ` that a dense train would give.
"""
function recurrent(t, φ, T, dur; dφ = 1.0)
    isfinite(T) || return fill(false, length(t))
    nother = max(1, floor(Int, dur / T))
    map(eachindex(t)) do i
        spins = Set{Int}(); far = 0
        for j in eachindex(t)
            abs(t[j] - t[i]) > T / 2 || continue
            far += 1
            abs(mod(φ[j] - φ[i] + 180, 360) - 180) < dφ && push!(spins, round(Int, (t[j] - t[i]) / T))
        end
        rec = length(spins); λ = far * 2dφ / 360
        rec >= cld(nother, 2) && rec >= 3λ + 1
    end
end

# Spin phases (deg) where sun-shadow spikes pile up beyond what `recurrent` removes (excess bins of the
# phase histogram of accepted spikes, all event days): at c = 45° mod 90°, (c - 2, c + 1.5); for c = 45°, 225°
# the peak splits to c ± 4.5° from ~2016 Q3 on, so (c ± 3.5, c ± 5.5) is masked too, in all years.
const ARTIFACT_PHASES = [[(c - 2.0, c + 1.5) for c in (45, 135, 225, 315)];
                         [(c + 3.5, c + 5.5) for c in (45, 225)]; [(c - 5.5, c - 3.5) for c in (45, 225)]]
masked(φ) = any(((lo, hi),) -> lo <= mod(φ, 360) < hi, ARTIFACT_PHASES)

"""
Pulses and trains of one burst. Candidates on E∥ (if b̂ known) and Ez, noise floor from quiet blocks
(`quiet_sigma`), are deduplicated within `nms`; those with iso ≥ p.iso and |a| ≥ p.amin are pulses. Each
gets a shape fit on its detection channel, local kurtosis (±kw), band-limited polarization, artifact flags
(`sat` clipping nearby, spin-recurrent `rec`, phase mask → `art`) and the `tdslike` verdict.
Pulses ≤ p.gap apart form trains, with features on the filtered E∥ (Ez without b̂) over the train ± isow:
`kurt`, `crest`, `duty` (fraction of samples ≥ 3σq), `clip` (fraction at the ADC rail: amplitudes are lower bounds when > 0), TDS-like pulse
count, shape fractions, width spread, `fce` (Hz) and `class` (`label`).
"""
function burst_pulses(b::Burst; p = PARAMS)
    E = copy(b.E); sat = copy(b.rail)
    for c in 1:3; sat .|= fillnan!(view(E, c, :)); end
    hpe = [highpass_chunks(E[c, :], b.t, b.fs, p.τ) for c in 1:3]
    hp = first.(hpe); edge = last(hpe[1])
    hasb = !any(isnan, b.bh)
    chans = Pair{Symbol, Vector{Float64}}[:Ez => hp[3]]
    hasb && pushfirst!(chans, :Epar => hp[1] .* b.bh[1, :] .+ hp[2] .* b.bh[2, :] .+ hp[3] .* b.bh[3, :])
    σq = Dict(ch => quiet_sigma(x, edge, b.fs; p) for (ch, x) in chans)
    m = max(1, round(Int, p.nms * b.fs)); ms = round(Int, p.satw * b.fs)
    cand = [(; c..., chan = ch) for (ch, x) in chans for c in spikes(x, b.fs, σq[ch], edge; p)
            if c.iso >= p.iso && abs(c.a) >= p.amin]
    sort!(cand, by = c -> -c.snr)
    taken = Int[]; P = []; K = round(Int, p.kw * b.fs)
    for c in cand
        any(t -> abs(t - c.i) <= m, taken) && continue
        push!(taken, c.i)
        x = Dict(chans)[c.chan]
        t0, s, lo, hi, fa, fb, R = shapefit(x, c.i, c.w, p.fitw * b.fs)
        pol, pair = polarization(hp, view(b.bh, :, c.i), b.φ[c.i], t0, s, lo, hi)
        y = view(x, max(1, c.i - K):min(length(x), c.i + K)); m2 = mean(abs2, y)
        push!(P, (i = c.i, t = b.t[c.i], chan = c.chan, a = c.a, snr = c.snr, w_ms = 1000c.w / b.fs, opp = c.opp,
            iso = c.iso, ctr = c.ctr, ac = c.ac, f_ac = c.aclag > 0 ? b.fs / c.aclag : NaN,
            fit_R = R, fit_s_ms = 1000s / b.fs, fit_uni = fa, fit_bip = fb, pol = pol, pair = pair,
            kurt = mean(v -> v^4, y) / m2^2,
            phi = b.φ[c.i], sat = any(view(sat, max(1, c.i - ms):min(length(sat), c.i + ms)))))
    end
    sort!(P, by = r -> r.i)
    T = 360 * (b.t[end] - b.t[1]) / sum(x -> mod(x + 180, 360) - 180, diff(b.φ))
    rec = recurrent([r.t for r in P], [r.phi for r in P], T, b.t[end] - b.t[1])
    P = map(enumerate(P)) do (k, r)
        # No single-pair cut near the shadow phases: a field-aligned pulse also sits on one pair whenever b̂ lies along
        # it, and outside the mask such pulses are no more frequent near the shadow phases than far from them.
        art = rec[k] || !isfinite(r.phi) || masked(r.phi)
        v = (r.kurt >= p.kurt) + (r.ctr >= p.ctr) + (r.fit_R <= p.rfit)
        shape = abs(r.fit_bip) >= 2abs(r.fit_uni) ? :bipolar : abs(r.fit_uni) >= 2abs(r.fit_bip) ? :unipolar : :mixed
        (; r..., rec = rec[k], art = art, shape = shape,
           tdslike = !art && !r.sat && isfinite(r.pol) && r.pol >= p.pol && v >= p.votes)
    end
    # trains
    x = first(chans).second; σx = σq[first(chans).first]; W = round(Int, p.isow * b.fs)
    trains = []; tid = zeros(Int, length(P))
    for (k, rk) in enumerate(find_continuous_indexranges([q.t for q in P], p.gap))
        g = P[rk]; tid[rk] .= k
        r = max(1, g[1].i - W):min(length(x), g[end].i + W); y = view(x, r)
        m2 = mean(abs2, y)
        L = [q for q in g if q.tdslike]; ok = [!q.art for q in g]
        j = argmax([abs(q.a) for q in g]); jl = isempty(L) ? 0 : argmax([abs(q.a) for q in L])
        r = (start_time = unix2datetime(g[1].t), end_time = unix2datetime(g[end].t), n = length(g), n_tds = length(L),
            peak = g[j].a, peak_chan = g[j].chan, peak_phi = g[j].phi,
            peak_tds = jl == 0 ? NaN : L[jl].a, peak_tds_Epar_ratio = jl == 0 ? NaN : L[jl].pol,
            clip = mean(view(b.rail, max(1, g[1].i - W):min(length(sat), g[end].i + W))),
            art_frac = 1 - mean(ok), sigma_q = σx, kurt = mean(v -> v^4, y) / m2^2, crest = maximum(abs, y) / sqrt(m2),
            duty = mean(v -> abs(v) >= 3σx, y),
            f_bip = isempty(L) ? NaN : mean(q.shape == :bipolar for q in L),
            f_uni = isempty(L) ? NaN : mean(q.shape == :unipolar for q in L),
            width_ms = isempty(L) ? NaN : median(q.w_ms for q in L),
            width_spread = length(L) >= 3 ? iqr([q.fit_s_ms for q in L]) / median([q.fit_s_ms for q in L]) : NaN,
            pol = isempty(L) ? NaN : median(q.pol for q in L), ac = median(q.ac for q in g), opp = median(q.opp for q in g),
            f_ac = median(q.f_ac for q in g), fce = 28.0 * median(view(b.bt, max(1, g[1].i - W):min(length(x), g[end].i + W))))
        push!(trains, (; r..., class = label(r, hasb; p)))
    end
    [(; r..., train = tid[k]) for (k, r) in enumerate(P)], trains, T
end

"""
Train class: artifact (mostly spin-recurrent / shadow pulses); tds_* when ≥ p.nmin TDS-like pulses, by the
majority fitted shape (bipolar: electron-hole-like; unipolar: double-layer-like, net potential); otherwise
no_b (no b̂, E∥ unknown), wave_coherent (periodic, median ac ≥ p.acmax) or wave_broadband.
"""
function label(r, hasb; p = PARAMS)
    r.art_frac >= p.art && return "artifact"
    if r.n_tds >= p.nmin
        r.f_bip >= p.shapefrac && return "tds_bipolar"
        r.f_uni >= p.shapefrac && return "tds_unipolar"
        return "tds_mixed"
    end
    hasb || return "no_b"
    r.ac >= p.acmax ? "wave_coherent" : "wave_broadband"
end
istds(c) = startswith(c, "tds")

"""
Bursts, pulses and trains of one day. `evd`: Lang events of the day; a burst is Lang-positive if any
[start, end + 1 s] overlaps it (Lang's times are 1 s floors).
"""
function process_day(d::ProbeDay, evd; p = PARAMS)
    bs = load_day(d); date = d.date
    brows = []; prows = []; trows = []
    lw = [(datetime2unix(a), datetime2unix(b) + 1) for (a, b) in zip(evd.start_time, evd.end_time)]
    for (k, b) in enumerate(bs)
        ids = [evd.id[j] for (j, (a, c)) in enumerate(lw) if a <= b.t[end] && c >= b.t[1]]
        P, Tr, T = burst_pulses(b; p)
        for q in P; push!(prows, (date = date, burst = k, q...)); end
        for (j, r) in enumerate(Tr); push!(trows, (date = date, burst = k, train = j, r...)); end
        push!(brows, (date = date, burst = k, t0 = b.t[1], t1 = b.t[end], fs = b.fs, Tspin = T,
            clip = mean(b.rail), hasb = !any(isnan, b.bh),
            lang = !isempty(ids), lang_ids = join(ids, ";")))
    end
    DataFrame(brows), DataFrame(prows), DataFrame(trows)
end

const DETECT = DailyProduct("detect", ["detect/detect.jl", "scripts/efw.jl"], [:efw, :fgm], process_day)

"Bursts, pulses and trains of the days of events `ev` (`process_day` with `PARAMS`, cached per day)."
function run_days(ev)
    res = only(daily([perday(DETECT, select(ev, :id, :start_time, :end_time), Date.(ev.start_time))]; files = dayfile))
    Tuple(reduce(vcat, DataFrame[isnothing(r) ? DataFrame() : r[k] for r in res]; cols = :union) for k in 1:3)
end
