# Dataset quality figures and the numbers quoted in report/report.typ → report/figs/, report/summary.json.
# `julia --project=. scripts/qa.jl` (after build.jl and bgwindow.jl)
using CSV, DataFrames, DataFramesMeta, Dates, Statistics, StatsBase
include("table.jl")

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
S["frac_phase_flag"] = mean(coalesce.(df.efw_phase_flag, false)[coalesce.(df.efw_n, 0) .> 0])
S["frac_bg_unsteady"] = mean(>(30), good(df.bg_rot))
S["region"] = Dict(g.region[1] => Dict("N" => nrow(g), "clip" => frac(>(0), good(g.efw_clip)),
    "bg_unsteady" => frac(>(30), good(g.bg_rot))) for g in groupby(df, :region))
S["debye_sheath"] = med(posgood(df.debye_m[df.region .== "sheath/SW"]))

S["det_status"] = Dict(k => count(==(k), df.det_status) / nrow(df) for k in unique(df.det_status))
let tw = df[df.det_status .== "tds", :], m = fin.(tw.tds_peak) .& fin.(tw.efw_Epar)
    S["tds_span_med"] = median(datetime2unix.(tw.tds_end) .- datetime2unix.(tw.tds_start))
    S["tds_efw_corr"] = corspearman(abs.(tw.tds_peak[m]), abs.(tw.efw_Epar[m]))
    S["tds_frac_by_region"] = Dict(g.region[1] => mean(g.det_status .== "tds") for g in groupby(df, :region))
end

# spin phase of the efw maxima: clusters are spin-synchronous artifacts (efw_phase_flag)
let f = Figure(size = (1000, 300))
    for (k, (c, flag, lab)) in enumerate([(:efw_phi_Epar, :efw_phase_flag, "max |E∥|"), (:efw_phi_Ez, :efw_phase_flag, "max |Ez|")])
        ax = Axis(f[1, k]; xlabel = "spin phase at $lab (deg)", ylabel = k == 1 ? "events" : "", xticks = 0:45:360)
        kept = .!coalesce.(df[!, flag], true)
        hist!(ax, good(df[!, c]); bins = 0:2:360, color = :gray, label = "all")
        hist!(ax, good(df[kept, c]); bins = 0:2:360, label = "not flagged")
        k == 1 && axislegend(ax; framevisible = false, labelsize = 10)
    end
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
