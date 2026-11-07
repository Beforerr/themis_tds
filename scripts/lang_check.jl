# Reproduce Lang's list quantities (idl/findthmbrst_v2.pro) from our data on random days → report/lang_check.json.
# Needs data/tds_tha.csv; fetches FBK L2 for the sampled days into data/cdf/fbk. `julia --project=. -t auto scripts/lang_check.jl`
using CSV, DataFrames, Statistics, Dates, Random, StatsBase, CDFDatasets, TimeseriesUtilities, JSON
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

function day_rows(day, g)
    d = ProbeDay("tha", day, dayfile)
    (isnothing(input(d, :efw)) || isnothing(input(d, :fgm))) && return NamedTuple[]
    win = [(datetime2unix(r.start_time), datetime2unix(r.end_time) + 1) for r in eachrow(g)]   # Lang's times are 1 s floors
    span = (minimum(first.(win)) - 30, maximum(first.(win)) + 30)
    fgs = readseries(d, :fgm, "tha_fgs_time", "tha_fgs_gsm", span)
    t = dayvar(d, :efw, "tha_efw_time")
    keep = filter(s -> any(w -> t[first(s)] <= w[2] && t[last(s)] >= w[1], win), find_continuous_indexranges(t, 0.5))
    tk = t[reduce(vcat, collect.(keep))]
    segs = find_continuous_indexranges(tk, 0.5)
    L = efw_events(d, win)   # our efw stats over the list interval instead of the efw window
    ion = isnothing(input(d, :gmom)) ? nothing : moments(d, "ptirf", span)
    fb = (f = fbk_file(day); isnothing(f) ? nothing : fbk_ew(f))
    map(enumerate(zip(eachrow(g), win))) do (i, (r, (a, b)))
        s = segs[argmin([abs(tk[first(s)] - a) for s in segs])]
        bs, be = tk[first(s)] - a, tk[last(s)] - b
        B30, = wmean(fgs, a - 30, a + 30)
        vx = isnothing(ion) ? NaN : (w = window(ion.V.t, a - 30, a + 30); isempty(w) ? NaN : (v = ion.V.x[1, w]; v[argmax(abs.(v))]))
        fw = isnothing(fb) ? Int[] : findall(x -> a - 30 <= x <= a + 30, fb.t)
        (id = r.id, start = r.start_time, inside = bs <= 1 && be >= -1.5, exact = abs(bs) <= 1.5 && abs(be) <= 1.5,
         lang_Bx = r.lang_Bx, lang_Bz = r.lang_Bz, Bx30 = B30[1], Bz30 = B30[3],
         lang_Vx = r.lang_Vx, Vx30 = vx, Vi_x = r.Vi_x,
         lang_Ez = r.lang_max_Ez, efw_Epar = abs(r.efw_Epar), list_Epar = abs(L.Epar[i]),
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
    me = ok(R.lang_Ez, R.list_Epar, R.efw_Epar)
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
        "rho_Ez" => rho(R.lang_Ez[me], R.efw_Epar[me]), "Ez_ratio" => median(R.efw_Epar[me] ./ R.lang_Ez[me]),
        "rho_Ez_list" => rho(R.lang_Ez[me], R.list_Epar[me]), "Ez_ratio_list" => median(R.list_Epar[me] ./ R.lang_Ez[me]),
        "n_Ez" => count(me),
        "rms_edc12_match" => mean(near.(R.lang_rms[f12], R.rms30[f12])), "n_rms_edc12" => count(f12),
        "rms_edc34_daymax" => mean(near.(R.lang_rms[f34], R.rms_day[f34])), "n_rms_edc34" => count(f34),
        "rms_const_by_year" => Dict(string(k.yr) => mean(v.const_day) for (k, v) in pairs(groupby(transform(g, :date => ByRow(year) => :yr), :yr))),
        "rms_const_all" => mean(g.const_day),
    )
    write(joinpath(@__DIR__, "..", "report", "lang_check.json"), JSON.json(S; pretty = true))
    S
end

abspath(PROGRAM_FILE) == @__FILE__() && display(lang_check())
