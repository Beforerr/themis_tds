# Reproduce Lang's list quantities (idl/findthmbrst_v2.pro) from our data on random days → report/lang_check.json.
# Needs data/tds_tha.csv; fetches FBK L2 for the sampled days into data/cdf/fbk. `julia --project=. -t auto scripts/lang_check.jl`
using CSV, DataFrames, Statistics, Dates, Random, LinearAlgebra, StatsBase, FFTA, CDFDatasets, TimeseriesUtilities, JSON
include("themis.jl"); include("efw.jl"); include("events.jl")

const NDAYS = 80
const EXAMPLE = DateTime(2011, 1, 2, 15, 34, 59)   # quoted in the report

function fbk_file(date)
    dir = mkpath(joinpath(CDFDIR, "fbk")); y = year(date); d = Dates.format(date, "yyyymmdd")
    f = findlast(x -> occursin("_$(d)_", x), readdir(dir))
    isnothing(f) || return joinpath(dir, readdir(dir)[f])
    base = "https://themis.ssl.berkeley.edu/data/themis/tha/l2/fbk/$y/"
    m = match(Regex("tha_l2_fbk_$(d)_v\\d+\\.cdf"), read(`curl -sf $base`, String))
    isnothing(m) && return nothing
    path = joinpath(dir, m.match)
    run(`curl -sf --retry 3 -o $path.part $base$(m.match)`); mv("$path.part", path)
    path
end

"Lang's rms_Ew series: √(3 Σ π²A²/8) over FBK bands 572, 144, 36 Hz of fb_edc12, else fb_edc34."
function fbk_ew(path)
    ds = CDFDataset(path)
    has(k) = haskey(ds, k) && length(ds[k*"_time"]) > 0
    k = has("tha_fb_edc12") ? "tha_fb_edc12" : has("tha_fb_edc34") ? "tha_fb_edc34" : return nothing
    y = Array(ds[k])
    (t = Array(ds[k*"_time"]), ew = vec(sqrt.(3 .* sum(y[2:4, :] .^ 2 .* π^2 ./ 8; dims = 1))), edc12 = k == "tha_fb_edc12")
end

"max |E∥| and |Ez| in [a, b] after an FFT band-pass f1–f2 per continuous segment; masked spin phases skipped."
function bandmax(t, E, ϕ, bh, segs, a, b; f1 = 100.0, f2 = 3000.0, masked)
    par = z = -Inf
    for r in segs
        (t[last(r)] < a || t[first(r)] > b || length(r) < 1024) && continue
        n = length(r); dt = median(diff(t[r]))
        y = [E[c, k] for k in r, c in 1:3]
        y[.!isfinite.(y)] .= 0; y .-= mean(y; dims = 1)
        F = mapslices(fft, complex.(y); dims = 1)
        F[[!(f1 <= min(i, n - i) / (n * dt) <= f2) for i in 0:n-1], :] .= 0
        yf = real.(mapslices(ifft, F; dims = 1))
        for (j, k) in enumerate(r)
            (a <= t[k] <= b && !masked(rad2deg(ϕ[k]))) || continue
            par = max(par, abs(dot(bh, view(yf, j, :)))); z = max(z, abs(yf[j, 3]))
        end
    end
    (par = par, z = z)
end

function day_rows(day, g)
    ef, fgf, gmf = dayfile("efw", day), dayfile("fgm", day), dayfile("gmom", day)
    (isnothing(ef) || isnothing(fgf)) && return NamedTuple[]
    fg = CDFDataset(fgf)
    fgs = readseries(fg, "tha_fgs_time", "tha_fgs_gsm")
    tl, Bl, Sl = Array(fg["tha_fgl_time"]), Array(fg["tha_fgl_dsl"]), Array(fg["tha_fgl_ssl"])
    φf, Tspin = spinphase(tl, Sl, Bl)
    ew = CDFDataset(ef); t = Array(ew["tha_efw_time"])
    win = [(datetime2unix(r.start_time), datetime2unix(r.end_time) + 1) for r in eachrow(g)]   # Lang's times are 1 s floors
    keep = filter(s -> any(w -> t[first(s)] <= w[2] && t[last(s)] >= w[1], win), find_continuous_indexranges(t, 0.5))
    idx = reduce(vcat, collect.(keep)); tk = t[idx]
    E, _, ϕ = efi_dsl(tk, Array(ew["tha_efw"])[:, idx], Array(ew["tha_efw_hed_time"]), Array(ew["tha_efw_hed_ac"]), φf, Tspin)
    segs = find_continuous_indexranges(tk, 0.5)
    ion = isnothing(gmf) ? nothing : moments(CDFDataset(gmf), "ptirf")
    fb = (f = fbk_file(day); isnothing(f) ? nothing : fbk_ew(f))
    masked(x) = abs(mod(x, 90) - 45) < 4 || any(w -> abs(mod(x - w + 180, 360) - 180) < 4, (165.0, 346.0))
    map(zip(eachrow(g), win)) do (r, (a, b))
        s = segs[argmin([abs(tk[first(s)] - a) for s in segs])]
        bs, be = tk[first(s)] - a, tk[last(s)] - b
        B30, = wmean(fgs, a - 30, a + 30)
        r1 = searchsortedfirst(tl, a):searchsortedlast(tl, b)
        bh = isempty(r1) ? fill(NaN, 3) : normalize(vec(median(Bl[:, r1]; dims = 2)))
        m = bandmax(tk, E, ϕ, bh, segs, a, b; masked)
        vx = isnothing(ion) ? NaN : (w = window(ion.V.t, a - 30, a + 30); isempty(w) ? NaN : (v = ion.V.x[1, w]; v[argmax(abs.(v))]))
        fw = isnothing(fb) ? Int[] : findall(x -> a - 30 <= x <= a + 30, fb.t)
        (id = r.id, start = r.start_time, inside = bs <= 1 && be >= -1.5, exact = abs(bs) <= 1.5 && abs(be) <= 1.5,
         lang_Bx = r.lang_Bx, lang_Bz = r.lang_Bz, Bx30 = B30[1], Bz30 = B30[3],
         lang_Vx = r.lang_Vx, Vx30 = vx, Vi_x = r.Vi_x,
         lang_Ez = r.lang_max_Ez, efw_Epar = abs(r.efw_Epar), bp_par = m.par, bp_z = m.z, flag = r.efw_phase_flag === true,
         lang_rms = r.lang_rms_Ew, rms30 = isempty(fw) ? NaN : maximum(fb.ew[fw]),
         rms_day = isnothing(fb) ? NaN : maximum(filter(isfinite, fb.ew)), edc12 = isnothing(fb) ? missing : fb.edc12)
    end
end

ok(v...) = [all(x -> !ismissing(x) && x isa Real && isfinite(x), t) for t in zip(v...)]
rho(x, y) = (m = ok(x, y); corspearman(Float64.(x[m]), Float64.(y[m])))

function lang_check()
    d = CSV.read(joinpath(@__DIR__, "..", "data", "tds_tha.csv"), DataFrame; types = Dict(:start_time => DateTime, :end_time => DateTime))
    ev = read_events()
    d = leftjoin(d, select(ev, :start_time, :end_time, :lang_Vx), on = [:start_time, :end_time])
    d.date = Date.(d.start_time)
    Random.seed!(1)
    days = [sample(setdiff(unique(d.date), [Date(EXAMPLE)]), NDAYS; replace = false); Date(EXAMPLE)]
    out = Vector{Vector{NamedTuple}}(undef, length(days))
    Threads.@threads for i in eachindex(days)
        out[i] = try
            day_rows(days[i], d[d.date .== days[i], :])
        catch e
            @warn "day failed" days[i] exception = e
            NamedTuple[]
        end
    end
    R = DataFrame(reduce(vcat, out))
    dB = max.(abs.(R.lang_Bx .- R.Bx30), abs.(R.lang_Bz .- R.Bz30))
    mv = ok(R.lang_Vx, R.Vx30) .& (R.lang_Vx .!= 0)
    me = ok(R.lang_Ez, R.bp_par, R.efw_Epar) .& .!R.flag
    f12 = coalesce.(R.edc12, false); f34 = coalesce.(.!R.edc12, false)
    near(x, y) = abs(x - y) <= 0.011    # Lang prints 2 decimals
    ex = only(eachrow(R[R.start .== EXAMPLE, :]))
    # day-constant rms_Ew from the list alone: all events of a ≥ 2-event day share one value
    g = transform(groupby(d, :date), nrow => :n, :lang_rms_Ew => (x -> length(unique(x)) == 1) => :one)
    g.const_day = (g.n .>= 2) .& g.one
    ev.plain = .!occursin.("QLP", ev.source)
    S = Dict(
        "n_days" => length(unique(Date.(R.start))), "n_events" => nrow(R),
        "window_inside" => mean(R.inside), "window_exact" => mean(R.exact),
        "B_within1" => mean(dB[ok(dB)] .<= 1), "B_maxdiff" => maximum(dB[ok(dB)]), "rho_B_x" => rho(R.lang_Bx, R.Bx30),
        "rho_V" => rho(R.lang_Vx[mv], R.Vx30[mv]), "V_med_absdiff" => median(abs.(R.lang_Vx[mv] .- R.Vx30[mv])),
        "V_ratio_q" => quantile(R.lang_Vx[mv] ./ R.Vx30[mv], [0.05, 0.5, 0.95]),
        "V_signflips" => count((sign.(R.lang_Vx[mv]) .!= sign.(R.Vx30[mv])) .& (abs.(R.Vx30[mv]) .> 50)), "n_V" => count(mv),
        "rho_V_eventmean" => rho(R.lang_Vx[mv], R.Vi_x[mv]),
        "example" => Dict("lang_Vx" => ex.lang_Vx, "Vx30" => ex.Vx30, "Vi_x" => ex.Vi_x),
        "frac_VxVyVz_plain" => mean(((ev.lang_Vx .== ev.lang_Vy) .& (ev.lang_Vy .== ev.lang_Vz))[ev.plain]),
        "rho_Ez_raw" => rho(R.lang_Ez[me], R.efw_Epar[me]), "rho_Ez_bp" => rho(R.lang_Ez[me], R.bp_par[me]),
        "Ez_ratio_bp" => median(R.bp_par[me] ./ R.lang_Ez[me]), "n_Ez" => count(me),
        "rms_edc12_match" => mean(near.(R.lang_rms[f12], R.rms30[f12])), "n_rms_edc12" => count(f12),
        "rms_edc34_daymax" => mean(near.(R.lang_rms[f34], R.rms_day[f34])), "n_rms_edc34" => count(f34),
        "rms_const_by_year" => Dict(string(k.yr) => mean(v.const_day) for (k, v) in pairs(groupby(transform(g, :date => ByRow(year) => :yr), :yr))),
        "rms_const_all" => mean(g.const_day),
    )
    write(joinpath(@__DIR__, "..", "report", "lang_check.json"), JSON.json(S; pretty = true))
    S
end

abspath(PROGRAM_FILE) == @__FILE__() && display(lang_check())
