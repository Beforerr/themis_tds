# Dataset quality figures and the numbers quoted in report/report.typ → report/figs/, report/summary.json.
# `julia --project=. scripts/qa.jl` (after build.jl and bgwindow.jl)
using CSV, DataFrames, DataFramesMeta, Dates, Statistics, StatsBase
include("table.jl"); include("efw.jl")

const ROOT = joinpath(@__DIR__, "..")
const FIG = mkpath(joinpath(ROOT, "report", "figs"))
save_fig(name, f) = (save(joinpath(FIG, name * ".png"), f; px_per_unit = 2); f)

df, _, version = load_table()
S = Dict{String, Any}("pipeline" => version)
S["n_events"] = nrow(df)
S["n_days"] = length(unique(Date.(df.start_time)))
S["frac_moments"] = mean(fin, df.Ni)
S["frac_fgm"] = mean(fin, df.B_t)
S["frac_efw"] = mean(x -> fin(x) && x > 0, coalesce.(df.efw_n, 0))
S["frac_burst_moments"] = mean(>=(2), coalesce.(df.mom_burst, 0))
S["mom_mask"] = Dict("ion_lowN" => mean(df.mom_mask .& 1 .> 0), "ele_lowN" => mean(df.mom_mask .& 2 .> 0),
    "ion_isotropic" => mean(df.mom_mask .& 4 .> 0))
S["frac_clip"] = mean(>(0), good(df.efw_clip))
S["frac_bg_unsteady"] = mean(>(30), good(df.bg_rot))
S["region"] = Dict(g.region[1] => Dict("N" => nrow(g), "clip" => frac(>(0), good(g.efw_clip)),
    "bg_unsteady" => frac(>(30), good(g.bg_rot))) for g in groupby(df, :region))
S["debye_sheath"] = med(posgood(df.debye_m[df.region .== "sheath/SW"]))

S["det_status"] = Dict(k => count(==(k), df.det_status) / nrow(df) for k in unique(df.det_status))
let tw = df[df.det_status .== "tds", :], m = fin.(tw.tds_peak) .& fin.(tw.efw_Epar)
    S["tds_span_med"] = median(datetime2unix.(tw.tds_end) .- datetime2unix.(tw.tds_start))
    S["tds_efw_corr"] = corspearman(abs.(tw.tds_peak[m]), abs.(tw.efw_Epar[m]))
    S["tds_pulses_q"] = quantile(tw.tds_n_pulses, [0.1, 0.5, 0.9])
    S["tds_multi_train"] = mean(tw.tds_n_trains .> 1)
    S["tds_efw_within20"] = mean(abs.(abs.(tw.efw_Epar[m]) ./ abs.(tw.tds_peak[m]) .- 1) .< 0.2)
    S["tds_frac_by_region"] = Dict(g.region[1] => mean(g.det_status .== "tds") for g in groupby(df, :region))
end

# waveform cleaning (scripts/efwstages.jl): one burst before/after, and the spin phase of max |E∥| at each stage;
# clusters at fixed phases are spin-synchronous artifacts
"One burst (scripts/efwstages.jl) in `g`: spin phase, E∥ calibrated and cleaned with their maxima, zooms on the cleaned max."
function fig_example!(g, X, ex; zoom = (0.1, 0.005))
    masked = [any(((lo, hi),) -> lo <= x < hi, EFW_MASK) for x in X.phase]
    m = findall(masked)
    runs = [(X.t[m[first(r)]], X.t[m[last(r)]]) for r in find_continuous_indexranges(m, 1)]
    X = transform(X, [:bp, :phase] => ByRow((v, p) -> any(((lo, hi),) -> lo <= p < hi, EFW_MASK) ? NaN : coalesce(v, NaN)) => :final)
    t0 = ex.tds_start
    tr = X.t[argmax(abs.(coalesce.(X.raw, 0)))]
    tf = Dates.value(ex.efw_t_Epar - t0) / 1000   # efw_t_Epar has ms precision: refine to the sample
    near = findall(x -> abs(x - tf) <= 2e-3, X.t)
    tf = X.t[near[argmax(abs.(replace(X.final[near], NaN => 0)))]]
    ax0 = Axis(g[1, 1]; ylabel = "spin phase\n(deg)", yticks = 0:180:360, xticklabelsvisible = false)
    hspan!(ax0, collect(first.(EFW_MASK)), collect(last.(EFW_MASK)); color = (:orange, 0.3))
    s = 1:40:nrow(X)
    scatter!(ax0, X.t[s], X.phase[s]; markersize = 1.5, color = :black)
    axs = [ax0]
    for (k, (c, lab, tm, col)) in enumerate([(:raw, "calibrated", tr, :red), (:final, "100–3000 Hz,\nmasked", tf, :blue)])
        ax = Axis(g[k+1, 1]; ylabel = "E∥ (mV/m)\n$lab", xlabel = k == 2 ? "s after $(Dates.format(t0, "yyyy-mm-dd HH:MM:SS.sss"))" : "",
            xticklabelsvisible = k == 2)
        vspan!(ax, first.(runs), last.(runs); color = (:orange, 0.3))
        lines!(ax, X.t, coalesce.(X[!, c], NaN); color = :black, linewidth = 0.4)
        i = argmin(abs.(X.t .- tm))
        scatter!(ax, [X.t[i]], [X[i, c]]; color = col, markersize = 8, marker = :circle, strokewidth = 0)
        push!(axs, ax)
    end
    for (tm, col) in ((tr, :red), (tf, :blue)), ax in axs
        vlines!(ax, tm; color = (col, 0.5), linewidth = 1, linestyle = :dash)
    end
    linkxaxes!(axs...)
    # A TDS is a train of ms pulses with E∥ ≫ E⊥; the calibrated trace checks the band-pass did not shape the pulse.
    rel(h) = (w = abs.(X.t .- tf) .<= h; (w, 1e3 .* (X.t[w] .- tf)))
    w, ms = rel(zoom[1])
    zg = GridLayout(g[1:3, 2])
    z1 = Axis(zg[1, 1]; title = "cleaned E∥ max: train", xlabel = "ms", ylabel = "mV/m")
    lines!(z1, ms, X.final[w]; color = :blue, linewidth = 0.6)
    w, ms = rel(zoom[2])
    z2 = Axis(zg[2, 1]; title = "single pulse", xlabel = "ms from the cleaned E∥ max", ylabel = "mV/m")
    lines!(z2, ms, coalesce.(X.raw[w], NaN) .- median(skipmissing(X.raw[w])); color = :gray, linewidth = 1, label = "E∥ calibrated")
    lines!(z2, ms, coalesce.(X.perp1[w], NaN); color = :orange, label = "E⊥1")
    lines!(z2, ms, coalesce.(X.perp2[w], NaN); color = :green, label = "E⊥2")
    lines!(z2, ms, X.final[w]; color = :blue, label = "E∥")
    axislegend(z2; framevisible = false, labelsize = 9, position = :rb, rowgap = 0)
    colsize!(g, 2, Relative(0.27))
    rowsize!(g, 1, Relative(0.16))
    g
end

let f = Figure(size = (1000, 800)), St = CSV.read(joinpath(ROOT, "data", "efw_stages.csv"), DataFrame), edges = 0:2:360
    X = CSV.read(joinpath(ROOT, "data", "efw_example.csv"), DataFrame)
    ex = only(eachrow(df[df.id .== X.id[1], :]))
    fig_example!(GridLayout(f[1, 1:3]), X, ex)
    i = argmax(abs.(coalesce.(X.raw, 0)))
    S["efw_example"] = Dict("date" => string(Date(ex.tds_start)), "raw_max" => abs(X.raw[i]), "raw_t" => X.t[i], "raw_phi" => X.phase[i],
        "max" => abs(ex.efw_Epar), "phi" => ex.efw_phi_Epar)
    late = Dict(r.id => r.start_time >= DateTime(2016, 7) for r in eachrow(df))
    final = DataFrame(id = df.id, stage = "final", phi = df.efw_phi_Epar)
    x = edges[1:end-1] .+ step(edges) / 2
    counts(v) = fit(Histogram, collect(good(v)), edges).weights
    for (k, (g, lab)) in enumerate([("raw", "calibrated"), ("band-pass", "+ 100–3000 Hz"), ("final", "+ phase mask (table)")])
        s = g == "final" ? final : St[St.stage .== g, :]
        L = [get(late, i, false) for i in s.id]
        ax = Axis(f[2, k]; title = lab, xlabel = "spin phase at max |E∥| (deg)", ylabel = k == 1 ? "events" : "", xticks = 0:90:360)
        vspan!(ax, collect(first.(EFW_MASK)), collect(last.(EFW_MASK)); color = (:orange, 0.3))
        early = counts(s.phi[.!L])
        barplot!(ax, x, early; gap = 0, color = :gray, label = "before 2016-07")
        barplot!(ax, x, counts(s.phi[L]); offset = early, gap = 0, color = :steelblue, label = "from 2016-07")
        k == 1 && axislegend(ax; framevisible = false, labelsize = 10, position = :rt)
    end
    for (r, l) in ((1, "a"), (2, "b"))
        Label(f[r, 1, TopLeft()], l; fontsize = 16, font = :bold, padding = (0, 10, 5, 0), halign = :right)
    end
    rowsize!(f.layout, 2, Relative(0.3))
    save_fig("efw_quality", f)
end

# consistency with Lang's list
let f = Figure(size = (1000, 330))
    for (k, (a, b, lab, lg)) in enumerate([(:lang_Bz, :B_z, "Bz GSM (nT)", false), (:lang_max_Ez, :efw_Epar, "max |E∥| (mV/m)", true)])
        sc = lg ? log10 : identity
        ax = Axis(f[1, k]; xlabel = "Lang: $lab", ylabel = "this work", xscale = sc, yscale = sc)
        x, y = Float64.(coalesce.(df[!, a], NaN)), Float64.(coalesce.(df[!, b], NaN))
        lg && (x = abs.(x); y = abs.(y))
        s = fin.(x) .& fin.(y) .& (lg ? (x .> 0) .& (y .> 0) : abs.(x) .< 200)
        scatter!(ax, x[s], y[s]; markersize = 2, color = (:black, 0.3))
        m = extrema([x[s]; y[s]]); lines!(ax, [m...], [m...]; color = :red)
        S["lang_corr_"*string(b)] = corspearman(x[s], y[s])
    end
    save_fig("lang_compare", f)
end

# background-window test (scripts/bgwindow.jl)
let w = CSV.read(joinpath(ROOT, "data", "bgwindow.csv"), DataFrame)
    w = leftjoin(w, select(@rsubset(w, :W == 300), :id, :V_par => :V_par_ref); on = :id)
    w.dVpar = abs.(w.V_par .- w.V_par_ref)
    m = @by(@rsubset(w, :region in REGIONS), [:region, :W], :ang_event = med(good(:ang_event)),
        :ang_halves = med(good(:ang_halves)), :dVpar = med(good(:dVpar)))
    L = stack(m, [:ang_event, :ang_halves, :dVpar]; variable_name = :metric)
    lab = ["ang_event" => "∠(b̂_W, interval B) (deg)", "ang_halves" => "∠(b̂ before, b̂ after) (deg)", "dVpar" => "|V∥(W) − V∥(5 min)| (km/s)"]
    p = data(L) * mapping(:W => (x -> x / 60) => "half-window W (min)", :value => "median";
        color = :region, col = :metric => renamer(lab)) * visual(ScatterLines)
    save_fig("bgwindow", draw(p, CSCALE; facet = (; linkxaxes = :none, linkyaxes = :none),
        axis = (; xscale = log10, xticks = [0.25, 0.5, 1, 2, 5, 10, 20]), figure = (; size = (1100, 330))))
    S["bgwindow"] = Dict(string(Int(k.W)) => Dict("ang_event" => med(good(g.ang_event)), "ang_halves" => med(good(g.ang_halves)),
        "halves_gt30" => mean(>(30), good(g.ang_halves)), "dVpar" => med(good(g.dVpar))) for (k, g) in pairs(groupby(w, :W)))
end

write_summary(joinpath(ROOT, "report", "summary.json"), S)
