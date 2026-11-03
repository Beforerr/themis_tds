# Full pipeline: pulses/trains/classes over all Lang event days, agreement with Lang and with the visual
# audit, report numbers (detect/out/summary.json) and figures (detect/figs/). `julia --project=. -t 8 detect/run.jl`
include(joinpath(@__DIR__, "detect.jl"))
include(joinpath(@__DIR__, "plots.jl"))
using JSON

const OUT = joinpath(@__DIR__, "out"); const FIGS = joinpath(@__DIR__, "figs")
f_b, f_p, f_t = joinpath.(OUT, ("bursts.csv", "pulses.csv", "detections.csv"))

ev = read_events()
tab = CSV.read(joinpath(@__DIR__, "..", "data", "tds_tha.csv"), DataFrame)
tab.date = Date.(tab.start_time)
dreg = Dict(r.date => r.region for r in eachrow(combine(groupby(tab, :date), :region => mode => :region)))
if !isfile(f_t)
    B, P, T = run_days(sort(unique(Date.(ev.start_time))), ev)
    B.region = [get(dreg, d, "unknown") for d in B.date]
    cls = combine(groupby(T, [:date, :burst]), :class => (c -> any(istds, c)) => :det,
        :class => (c -> count(istds, c)) => :ntds_trains, :n_tds => sum => :ntds_pulses)
    B = leftjoin(B, cls, on = [:date, :burst]); sort!(B, [:date, :burst])
    B.det = coalesce.(B.det, false); B.ntds_trains = coalesce.(B.ntds_trains, 0); B.ntds_pulses = coalesce.(B.ntds_pulses, 0)
    CSV.write(f_b, B); CSV.write(f_p, P); CSV.write(f_t, T)
end
B = CSV.read(f_b, DataFrame; types = Dict(:lang_ids => String)); B.lang_ids = coalesce.(B.lang_ids, "")
# explicit DateTime: CSV.jl on Julia 1.14 infers ISO times with ms as Dates.Timestamp{Nanosecond}
P = CSV.read(f_p, DataFrame); T = CSV.read(f_t, DataFrame; types = Dict(:start_time => DateTime, :end_time => DateTime))
T = leftjoin(T, select(B, :date, :burst, :lang, :region), on = [:date, :burst])

# ---- agreement with Lang (burst level) ----
cm(lang, det) = (tp = count(lang .& det), fp = count(.!lang .& det), fn = count(lang .& .!det), tn = count(.!lang .& .!det))
pr(c) = (; c..., precision = c.tp / (c.tp + c.fp), recall = c.tp / (c.tp + c.fn))
Ttds = T[istds.(T.class), :]
evhit = let tw = Dict{Date, Vector{Tuple{Float64, Float64}}}()
    for r in eachrow(Ttds); push!(get!(tw, r.date, Tuple{Float64, Float64}[]), (datetime2unix(r.start_time), datetime2unix(r.end_time))); end
    count(eachrow(ev)) do r
        a, b = datetime2unix(r.start_time), datetime2unix(r.end_time) + 1
        any(((x, y),) -> x <= b && y >= a, get(tw, Date(r.start_time), Tuple{Float64, Float64}[]))
    end
end
# re-threshold from the pulse table: TDS-like pulses with |a| ≥ amin, ≥ nmin per non-artifact train
sweep = let g = combine(groupby(P, [:date, :burst, :train]), :art => mean => :af,
        [:tdslike, :a] => ((l, a) -> (n1 = count(l .& (abs.(a) .>= 1)), n2 = count(l .& (abs.(a) .>= 2)), n5 = count(l .& (abs.(a) .>= 5)))) => AsTable)
    map(Iterators.product((1.0, 2.0, 5.0), (1, 3, 10))) do (amin, nmin)
        h = g[(g.af .< PARAMS.art) .& (g[!, Symbol("n", Int(amin))] .>= nmin), :]
        s = Set(zip(h.date, h.burst)); d = [(a, b) in s for (a, b) in zip(B.date, B.burst)]
        (; amin, nmin, pr(cm(B.lang, d))...)
    end |> vec |> DataFrame
end

# ---- visual audit ----
audit = CSV.read(joinpath(@__DIR__, "audit.csv"), DataFrame)
audit = leftjoin(audit, select(B, :date, :burst, :det, :lang), on = [:date, :burst])
grp(c) = c in ("tds", "weak") ? "TDS" : c in ("ambiguous", "none") ? "amb" : "neg"
audit.blind = grp.(audit.class)
conf(df, col) = Dict(string(g, "/", d ? "det" : "nodet") => count((df[!, col] .== g) .& (df.det .== d)) for g in ("TDS", "neg", "amb"), d in (true, false))
lconf(df, col) = Dict(string(g, "/", l ? "lang" : "nolang") => count((df[!, col] .== g) .& (df.lang .== l)) for g in ("TDS", "neg", "amb"), l in (true, false))
jconf(df, col) = Dict(g => Dict(k => count((df[!, col] .== g) .& (df.det .== d) .& (df.lang .== l))
    for (k, d, l) in (("both", true, true), ("det", true, false), ("lang", false, true), ("none", false, false))) for g in ("TDS", "neg", "amb"))
audit = leftjoin(audit, combine(groupby(T, [:date, :burst]), :class => (c -> join(["$(n)×$(x)" for (x, n) in
    sort(collect(countmap(string.(c))), by = last, rev = true)], " ")) => :classes), on = [:date, :burst])
audit.classes = coalesce.(audit.classes, "")
# adjudication shift: +1 when the relabel moves towards the detector's verdict (agree > amb > disagree)
score(g, d) = g == "amb" ? 0 : (g == "TDS") == d ? 1 : -1
shift = [sign(score(a, d) - score(b, d)) for (a, b, d) in zip(audit.adjudicated, audit.blind, audit.det)]

okp = P[P.tdslike, :]
summary = Dict(
    "n_days" => length(unique(B.date)), "n_bursts" => nrow(B), "n_lang_bursts" => count(B.lang),
    "burst_hours" => sum(B.t1 .- B.t0) / 3600, "n_events" => nrow(ev), "n_events_hit" => evhit,
    "n_pulses" => nrow(P), "n_tdslike" => nrow(okp), "n_trains" => nrow(T),
    "class_trains" => Dict(g.class[1] => nrow(g) for g in groupby(T, :class)),
    "class_trains_lang" => Dict(g.class[1] => mean(g.lang) for g in groupby(T, :class)),
    "n_det_bursts" => count(B.det), "nob_bursts" => count(.!B.hasb), "frac_masked_phase" => sum(((a, b),) -> b - a, ARTIFACT_PHASES) / 360,
    "tds_clip_frac" => mean(Ttds.clip .> 0),
    "train_dur_med" => median(Dates.value.(Ttds.end_time .- Ttds.start_time) ./ 1000), "train_n_med" => median(Ttds.n_tds),
    "pulse_w_med" => median(okp.w_ms), "pulse_bip" => mean(okp.shape .== "bipolar"), "pulse_uni" => mean(okp.shape .== "unipolar"),
    "pulse_amp_q" => quantile(abs.(okp.a), [0.5, 0.9, 0.99]),
    "all" => pr(cm(B.lang, B.det)),
    "region" => Dict(g.region[1] => pr(cm(g.lang, g.det)) for g in groupby(B, :region)),
    "sweep" => [Dict(pairs(r)) for r in eachrow(sweep)],
    "audit_blind" => Dict(s => conf(audit[audit.set .== s, :], :blind) for s in unique(audit.set)),
    "audit_adj" => Dict(s => conf(audit[audit.set .== s, :], :adjudicated) for s in unique(audit.set)),
    "audit_blind_all" => conf(audit, :blind), "audit_adj_all" => conf(audit, :adjudicated),
    "audit_lang_blind" => lconf(audit, :blind), "audit_lang_adj" => lconf(audit, :adjudicated),
    "audit_sheath_lang_adj" => lconf(audit[audit.set .== "sheath", :], :adjudicated),
    "audit_joint_blind" => jconf(audit, :blind), "audit_joint_adj" => jconf(audit, :adjudicated),
    "audit_adj_shift" => Dict("toward" => count(==(1), shift), "away" => count(==(-1), shift)),
    "audit_lang_fp_eye" => countmap(audit.class[(audit.adjudicated .== "neg") .& audit.lang]),
    "audit_det_errors" => [Dict(pairs(r)) for r in eachrow(select(audit[((audit.adjudicated .== "TDS") .& .!audit.det) .|
        ((audit.adjudicated .== "neg") .& audit.det), :], :set, :date, :burst, :class, :adjudicated, :det, :lang, :classes))],
    "class_region" => Dict(g.region[1] => Dict(h.class[1] => nrow(h) for h in groupby(g, :class)) for g in groupby(T, :region)),
    "class_clip" => Dict(g.class[1] => mean(g.clip .> 0) for g in groupby(T, :class)),
    "burst_clip" => Dict("all" => mean(B.clip .> 0), "lang" => mean(B.clip[B.lang] .> 0), "det" => mean(B.clip[B.det] .> 0)),
    "tds_clip_frac_1pct" => mean(Ttds.clip .> 0.01),
    "params" => Dict(pairs(PARAMS)),
    "cuts" => Dict("art" => count(P.art), "sat" => count(P.sat), "rec" => count(P.rec), "shadow1" => count(P.shadow1),
        "pol" => count(.!P.art .& .!(coalesce.(P.pol, NaN) .>= PARAMS.pol))),
)
open(io -> JSON.print(io, summary, 1), joinpath(OUT, "summary.json"), "w")

# ---- figures ----
burst_of(t::DateTime) = (bs = load_day(Date(t)); u = datetime2unix(t); (Date(t), findfirst(b -> b.t[1] - 1 <= u <= b.t[end] + 1, bs)))
function lab(d, k)
    r = B[(B.date .== d) .& (B.burst .== k), :][1, :]
    c = T[(T.date .== d) .& (T.burst .== k), :class]
    "$(r.region), Lang $(r.lang ? "yes" : "no"), detector $(r.det ? "TDS" : "no TDS"): " *
        join(["$(n)×$(x)" for (x, n) in sort(collect(countmap(c)), by = last, rev = true)], " ")
end
ex1 = [(burst_of(ev.start_time[182])..., "TDS train, Ez sun-shadow pulses rejected"),
       (Date(2015, 6, 10), 9, "Lang miss: E∥ bipolar TDS"),
       (Date(2020, 9, 8), 4, "Lang false positive: only sun-shadow pulses"),
       (Date(2018, 11, 30), 8, "sheath: isolated E∥ pulses in weak background")]
ex2 = [(Date(2015, 3, 29), 6, "sheath broadband bursts"),
       (Date(2021, 3, 15), 7, "coherent packet (in Lang's list)"),
       (Date(2018, 8, 18), 8, "sheath: TDS on broadband background"),
       (Date(2021, 6, 4), 4, "spin-plane shadow pulses (in Lang's list)")]
for (f, ex) in (("examples_a.png", ex1), ("examples_b.png", ex2))
    gallery([(d, k, "$t  [$(lab(d, k))]") for (d, k, t) in ex], P, ev, joinpath(FIGS, f))
end

let fig = Figure(size = (1000, 600), fontsize = 13)
    late = P.date .>= Date(2016, 7, 1); fin = isfinite.(P.phi)
    for (j, (sel, name)) in enumerate(((.!late, "before 2016-07"), (late, "from 2016-07")))
        ax = Axis(fig[j, 1], ylabel = "pulses / 0.5°", title = name, xlabel = j == 2 ? "spin phase φ (deg)" : "")
        for (m, c, l) in ((sel, :gray60, "all pulses"), (sel .& .!P.rec .& .!P.shadow1, :steelblue, "− recurrent − single-pair shadow"), (sel .& P.tdslike, :black, "TDS-like"))
            stairs!(ax, 0.25:0.5:359.75, fit(Histogram, mod.(P.phi[m .& fin], 360), 0:0.5:360).weights, color = c, label = l, step = :center)
        end
        for (a, b) in ARTIFACT_PHASES; vspan!(ax, a, b, color = (:red, 0.15)); end
        ax.xticks = 0:45:360; j == 1 && axislegend(ax, position = :rt, labelsize = 10)
    end
    save(joinpath(FIGS, "phase.png"), fig)
end

let fig = Figure(size = (1000, 380), fontsize = 13)
    ax = Axis(fig[1, 1], xlabel = "recall vs Lang", ylabel = "precision vs Lang", title = "burst level, threshold sweep")
    for (nmin, c) in zip((1, 3, 10), (:orange, :black, :steelblue))
        r = sweep[sweep.nmin .== nmin, :]
        scatterlines!(ax, r.recall, r.precision, color = c, label = "≥ $nmin TDS-like pulses/train")
        nmin == 3 && text!(ax, r.recall, r.precision, text = ["$(a) mV/m" for a in r.amin], fontsize = 10, offset = (4, 2))
    end
    axislegend(ax, position = :lb, labelsize = 10); xlims!(ax, nothing, 1.05)
    regs = ["PS", "PSBL", "lobe", "sheath/SW", "dayside/inner"]
    ax2 = Axis(fig[1, 2], xticks = (1:5, regs), xticklabelrotation = π / 6, ylabel = "fraction", title = "per region")
    v = [summary["region"][r] for r in regs]
    barplot!(ax2, [1:5; 1:5], [getfield.(v, :precision); getfield.(v, :recall)], dodge = [fill(1, 5); fill(2, 5)],
        color = [fill(:gray, 5); fill(:black, 5)])
    Legend(fig[1, 3], [PolyElement(color = :gray), PolyElement(color = :black)], ["precision", "recall"])
    save(joinpath(FIGS, "agreement.png"), fig)
end

# pulse features by class of the train they belong to (not artifact)
let fig = Figure(size = (1000, 600), fontsize = 13)
    Q = innerjoin(P[.!P.art, :], select(T, :date, :burst, :train, :class), on = [:date, :burst, :train])
    Q.grp = [istds(c) ? "TDS train" : c for c in Q.class]
    grps = ["TDS train", "wave_broadband", "wave_coherent"]
    specs = [(Q.kurt, "local kurtosis (±25 ms)", log10, 0.3:0.04:2.3), (Q.ctr, "contrast |a| / local σ", log10, 0.3:0.03:1.8),
             (Q.fit_R, "shape-fit residual fraction", identity, 0:0.02:1), (Q.pol, "|E∥|/|E⊥| (in band)", log10, -1.5:0.05:1.5)]
    for (i, (x, l, f, e)) in enumerate(specs)
        ax = Axis(fig[(i - 1) ÷ 2 + 1, (i - 1) % 2 + 1], xlabel = f === log10 ? "log10 " * l : l, ylabel = "pdf")
        for (g, c) in zip(grps, Makie.wong_colors())
            y = f.(x[(Q.grp .== g) .& isfinite.(x) .& (x .> 0)])
            stairs!(ax, e[1:end-1] .+ step(e) / 2, normalize(fit(Histogram, y, e), mode = :pdf).weights, color = c, label = g, step = :center)
        end
        th = (log10(PARAMS.kurt), log10(PARAMS.ctr), PARAMS.rfit, log10(PARAMS.pol))[i]
        vlines!(ax, th, color = :red, linestyle = :dash)
        i == 1 && axislegend(ax, labelsize = 10)
    end
    save(joinpath(FIGS, "features.png"), fig)
end

# TDS-like pulse properties by region
let fig = Figure(size = (1000, 330), fontsize = 13), A = leftjoin(okp, select(B, :date, :burst, :region), on = [:date, :burst])
    regs = ["PS", "PSBL", "sheath/SW", "dayside/inner"]
    specs = [(abs.(A.a), "log10 peak |E| (mV/m, filtered)", log10, 0:0.05:2.5),
             (A.w_ms, "half-max lobe width (ms)", identity, 0.0305:0.12207:2.5),   # bins = 8192 S/s quantum
             (A.pol, "log10 |E∥|/|E⊥|", log10, 0:0.05:2)]
    for (i, (x, l, f, e)) in enumerate(specs)
        ax = Axis(fig[1, i], xlabel = l, ylabel = "pdf")
        for (r, c) in zip(regs, Makie.wong_colors())
            y = f.(x[(A.region .== r) .& isfinite.(x)])
            stairs!(ax, e[1:end-1] .+ step(e) / 2, normalize(fit(Histogram, y, e), mode = :pdf).weights, color = c, label = r, step = :center)
        end
        i == 1 && axislegend(ax, labelsize = 9)
    end
    save(joinpath(FIGS, "pulses.png"), fig)
end
