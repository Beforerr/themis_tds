# Build the event table
# Check a change on a few days against that table: scripts/build.jl --sample[=24]
using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Serialization, Parquet2, JSON
include("events.jl"); include("themis.jl"); include("omni.jl")
include("efw.jl"); include("compare.jl"); include("tdswindow.jl")

const μ0, ε0, qe, mp = 4π * 1e-7, 8.8541878128e-12, 1.602176634e-19, 1.67262192e-27

"Shue et al. (1998) magnetopause radius [RE] at angle θ from +X, for Bz [nT] and Pdyn [nPa]."
function shue98(θ, Bz, Pd)
    r0 = (10.22 + 1.29tanh(0.184(Bz + 8.14))) * Pd^(-1 / 6.6)
    α = (0.58 - 0.007Bz) * (1 + 0.024log(Pd))
    r0 * (2 / (1 + cos(θ)))^α
end

function derive!(df)
    # n < 1e-3 cm⁻³ comes from ESA-less (SST-only) moment sums, with T ~ 20–30 keV
    df.mom_mask = zeros(Int, nrow(df))
    for (bit, n, cols) in ((2, :Ne, (:Ne, :Te, :Te_par, :Te_perp)), (1, :Ni, (:Ni, :Ti, :Ti_par, :Ti_perp, :Vi_x, :Vi_y, :Vi_z, :V_par_5min, :V_perp_5min, :Vperp_5min_x, :Vperp_5min_y, :Vperp_5min_z, :V_par_30s, :V_perp_30s)))
        bad = df[!, n] .< 1e-3
        df.mom_mask[bad] .|= bit
        for c in cols
            df[bad, c] .= NaN
        end
    end
    # Reduced-mode ion moments in some (mostly sheath) intervals carry no angular information: |V| ~ 1e-3 km/s, T⊥/T∥ ≡ 1
    iso = @. sqrt(df.Vi_x^2 + df.Vi_y^2 + df.Vi_z^2) < 1
    df.mom_mask[iso] .|= 4
    for c in (:Ti_par, :Ti_perp, :Vi_x, :Vi_y, :Vi_z, :V_par_5min, :V_perp_5min, :Vperp_5min_x, :Vperp_5min_y, :Vperp_5min_z, :V_par_30s, :V_perp_30s)
        df[iso, c] .= NaN
    end
    pB = @. (df.B_t * 1e-9)^2 / 2μ0                     # magnetic pressure, Pa
    @. df.beta_i = df.Ni * 1e6 * df.Ti * qe / pB
    @. df.beta_e = df.Ne * 1e6 * df.Te * qe / pB
    @. df.beta = df.beta_i + df.beta_e
    @. df.debye_m = sqrt(ε0 * df.Te / (df.Ne * 1e6 * qe))
    @. df.fpe_fce = 8980sqrt(df.Ne) / (27.99df.B_t)
    @. df.VA_kms = df.B_t * 1e-9 / sqrt(μ0 * df.Ni * 1e6 * mp) / 1e3
    @. df.theta_B = atand(df.B_z, hypot(df.B_x, df.B_y))
    @. df.theta_Bbg = atand(df.Bbg_z, hypot(df.Bbg_x, df.Bbg_y))
    @. df.Bx_over_Bt = abs(df.B_x) / df.B_t
    @. df.Te_Ti = df.Te / df.Ti
    @. df.Ti_aniso = df.Ti_perp / df.Ti_par
    @. df.Te_aniso = df.Te_perp / df.Te_par
    @. df.Epar_norm = abs(df.efw_Epar) * 1e-3 * df.debye_m / df.Te
    @. df.Ez_norm = abs(df.efw_Ez) * 1e-3 * df.debye_m / df.Te
    @. df.Vi_t = sqrt(df.Vi_x^2 + df.Vi_y^2 + df.Vi_z^2)
    r = @. sqrt(df.gsm_x^2 + df.gsm_y^2 + df.gsm_z^2)
    @. df.r_over_rmp = r / shue98(acos(df.gsm_x / r), df.BZ_GSM, df.Pressure)
    df.region = region.(df.r_over_rmp, df.Ni, df.Ti, df.beta, df.gsm_x)
    df
end

"Heuristic plasma region; thresholds follow common usage (β_PS > 0.5, β_PSBL 0.05–0.5; e.g. Baumjohann et al. 1989)."
function region(rr, n, Ti, β, x)
    isfinite(n) && isfinite(β) || return "unknown"
    (rr > 1 || (n > 5 && Ti < 1000)) && return "sheath/SW"
    x > -5 && return "dayside/inner"
    β > 0.5 ? "PS" : β > 0.05 ? "PSBL" : "lobe"
end

# (columns, unit, description). One row per group in data/columns.csv. "Interval" is the list interval [start_time, end_time]
# (Lang's wave-burst interval, tc its center); the TDS window is [tds_start, tds_end]; the efw window is given by efw_window.
# Columns copied from Lang Qin's list (idl/findthmbrst_v2.pro); all others are ours.
const LANG_COLUMNS = Set([:start_time, :end_time, :duration, :radius, :gsm_x, :gsm_y, :gsm_z, :mlt, :mlat,
    :lang_max_Ez, :lang_rms_Ew, :lang_Bx, :lang_By, :lang_Bz])
const COLUMNS = [
    ([:id], "", "event (list interval) index, sorted by start time"),
    ([:start_time, :end_time], "UT", "wave-burst interval containing the TDS (not the TDS itself; see tds_start/tds_end), floored to 1 s"),
    ([:duration], "s", "end − start (unfloored)"),
    ([:radius], "RE", "geocentric distance, burst mean"),
    ([:gsm_x, :gsm_y, :gsm_z], "RE", "position GSM, burst mean"),
    ([:mlt, :mlat], "h, deg", "MLT and magnetic latitude: SM, mean of burst start and end"),
    ([:lang_max_Ez, :lang_rms_Ew], "mV/m", "Lang's max |E∥| (100–3000 Hz, unrecorded sub-window) and FBK rms E (start ±30 s; whole-day max on many 2016–2022 days); report.typ, Lang's list section"),
    ([:lang_Bx, :lang_By, :lang_Bz], "nT", "Lang's B GSM: mean fgs over start ±30 s"),
    ([:B_x, :B_y, :B_z, :B_t], "nT", "FGM fgl GSM B, component medians over the interval (fgs ±1.5 s if no fgl); B_t = |B|"),
    ([:Bbg_x, :Bbg_y, :Bbg_z, :Bbg_t], "nT", "background B GSM: fgs medians over tc ±5 min"),
    ([:theta_B, :theta_Bbg], "deg", "elevation atan(Bz/√(Bx²+By²)) of interval / background B"),
    ([:Bx_over_Bt], "", "|Bx|/|B|: ~1 in lobe/PSBL, ~0 at the neutral sheet"),
    ([:dBz_1min], "nT", "median Bz in [tc, tc+1 min] − [tc−1 min, tc]; dipolarization indicator"),
    ([:Ni, :Ti, :Ti_par, :Ti_perp], "cm⁻³, eV, eV, eV", "ion density, average/parallel/perpendicular temperature: mean of spin samples within the interval ±1.5 s; GMOM ESA+SST burst (ptibb), else reduced (ptirf)"),
    ([:Vi_x, :Vi_y, :Vi_z, :Vi_t], "km/s", "ion bulk velocity GSM and |V|, same samples"),
    ([:V_par_5min, :V_perp_5min], "km/s", "V·b̂ and |V⊥|, b̂ from the ±5-min background B (steady magnetosphere; window-dependent if bg_rot > 30°)"),
    ([:Vperp_5min_x, :Vperp_5min_y, :Vperp_5min_z], "km/s", "V⊥ GSM, same b̂"),
    ([:V_par_30s, :V_perp_30s], "km/s", "same with b̂ from fgs medians over tc ±30 s (sheath, where the field rotates within a minute)"),
    ([:bg_rot], "deg", "angle between fgs medians over [tc−5 min, tc] and [tc, tc+5 min]; large ⇒ background window straddles a boundary"),
    ([:Vx_peak_5min], "km/s", "ion Vx (reduced mode) with largest |Vx| within tc ±5 min, 3-spin running median (flow-burst indicator)"),
    ([:Ne, :Te, :Te_par, :Te_perp], "cm⁻³, eV, eV, eV", "electron moments, same samples; GMOM burst (ptebb), else reduced (pterf)"),
    ([:beta_i, :beta_e, :beta], "", "2μ0 n T / B² for ions, electrons, total"),
    ([:debye_m], "m", "electron Debye length √(ε0 Te / n_e e²)"),
    ([:fpe_fce], "", "f_pe/f_ce"),
    ([:VA_kms], "km/s", "Alfvén speed (protons)"),
    ([:Te_Ti, :Ti_aniso, :Te_aniso], "", "Te/Ti, Ti⊥/Ti∥, Te⊥/Te∥ (Te_aniso meaningful only when mom_burst & 2)"),
    ([:flag_ion, :flag_ele], "", "OR of GMOM quality flags over the moment samples (-1: no data)"),
    ([:mom_burst], "", "moment source bits: 1 = ion burst, 2 = electron burst; else reduced mode"),
    ([:mom_mask], "", "masked-moment bits: 1 = ion, 2 = electron n < 1e-3 cm⁻³ (SST-only sums, all moments NaN); 4 = ion V ≈ 0 and T⊥ ≡ T∥ (no angular information; V, T∥, T⊥ NaN)"),
    ([:n_fgm, :n_fgm_bg, :n_ion, :n_ele], "", "sample counts: FGM interval, FGM background, ion and electron moments"),
    ([:BX_GSE, :BY_GSM, :BZ_GSM], "nT", "OMNI 1-min IMF (bow-shock-nose time), mean of the minutes overlapping the interval"),
    ([:BZ_GSM_1h], "nT", "IMF Bz mean over the hour before the interval"),
    ([:flow_speed, :proton_density, :Pressure], "km/s, cm⁻³, nPa", "solar wind speed, density, dynamic pressure, same minutes"),
    ([:AE_INDEX, :AL_INDEX, :SYM_H], "nT", "AE, AL, SYM-H, same minutes"),
    ([:AL_min_1h, :AE_max_1h], "nT", "min AL, max AE in the hour before the interval"),
    ([:Kp], "", "Kp (OMNI2 hourly) at start_time"),
    ([:r_over_rmp], "", "r / Shue98 magnetopause radius at the same angle"),
    ([:region], "", "heuristic region: sheath/SW, dayside/inner (X>−5), PS, PSBL, lobe"),
    ([:efw_window], "", "window of all efw_* stats: tds = detector TDS window (tds_start–tds_end), list = Lang's [start, end + 1 s] where no TDS was detected"),
    ([:efw_n, :efw_ac], "", "# efw waveform samples in the window (8192 or 16384 S/s), fraction AC-coupled"),
    ([:efw_t_Ez, :efw_Ez, :efw_Eperp_Ez, :efw_ratio_Ez, :efw_phi_Ez], "UT, mV/m, mV/m, –, deg",
        "at max |Ez| (DSL spin axis; own L1 calibration, 100–3000 Hz band-pass, boom-shadow phases masked): time, signed Ez, in-band |E⊥B|, |Ez|/|E⊥|, spin phase"),
    ([:efw_t_Epar, :efw_Epar, :efw_Eperp_Epar, :efw_ratio_Epar, :efw_phi_Epar], "UT, mV/m, mV/m, –, deg",
        "same at max |E∥|, E∥ = E·b̂ with b̂ the window-median B"),
    ([:efw_rms_Ez, :efw_rms_Exy], "mV/m", "in-band rms Ez and spin-plane |E| over the window"),
    ([:efw_clip], "", "fraction of samples with a raw channel at the ADC rail (±32767, or −32768 = fill); > 0 ⇒ maxima are lower bounds (spin-plane AC saturates near 70 mV/m)"),
    ([:Epar_norm, :Ez_norm], "", "e|E∥|λ_D/T_e (≈ eφ/T_e for a Debye-scale structure), same with |Ez|"),
    ([:det_status], "", "automatic detector (detect/report.typ) on the bursts overlapping [start, end + 1 s]: tds (the event owns TDS pulses, see tds_tha_pulses), no_tds, no_b (no b̂, never TDS), no_burst (none ≥ 4096 samples)"),
    ([:tds_start, :tds_end], "UT", "first to last of the event's TDS pulses ± 10 ms: the TDS window, which Lang's list does not record"),
    ([:tds_n_trains, :tds_n_pulses, :tds_class], "", "number of TDS trains holding the event's TDS pulses, and of those pulses; class (bipolar/unipolar/mixed) of the train holding most of them"),
    ([:tds_peak, :tds_Epar_ratio, :tds_clip], "mV/m, –, –", "largest of the event's TDS pulses (signed, E∥ or Ez, 150 Hz high-pass) and its in-band |E∥|/|E⊥|; max rail-sample fraction near its trains (> 0 ⇒ lower bound)"),
]

# data/tds_tha_pulses.parquet: one row per TDS pulse, with its owning event (event_pulses); detector definitions in detect/detect.jl.
const PULSE_COLUMNS = [
    ([:id], "", "event (row of tds_tha) owning the pulse: the one whose [start, end + 1 s] contains it, else the nearest one its TDS train overlaps"),
    ([:time, :t_unix], "UT, s", "pulse peak time (ms precision), and as unix seconds (µs precision)"),
    ([:burst, :train, :train_class], "", "burst number within the day and train number within the burst (as in detect/out), class of the train"),
    ([:chan], "", "detection channel: Epar (E·b̂) or Ez (spin axis, without b̂)"),
    ([:amplitude], "mV/m", "signed peak on the 150 Hz high-passed detection channel"),
    ([:width_ms, :fit_s_ms], "ms", "half-maximum width of the peak lobe; σ of the best Gaussian / Gaussian-derivative fit"),
    ([:pulse_shape], "", "bipolar (Gaussian-derivative fit dominates), unipolar (Gaussian dominates) or mixed"),
    ([:Epar_ratio], "", "|E∥|/|E⊥| of the fitted template in the pulse's own band"),
    ([:phi], "deg", "spin phase at the peak"),
    ([:snr, :iso, :ctr], "", "peak over the burst's quiet-time σ; peak over the rms within ±10 ms excluding the pulse; peak over the median-based σ there"),
    ([:kurt, :fit_R], "", "kurtosis within ±25 ms; residual fraction of the template fit (≈ 0 for an isolated single structure)"),
]

const CACHE = joinpath(@__DIR__, "..", "data", "cache")

"""Stage output cached in data/cache, keyed by the content of `deps` (stage source files and the event list),
so iterating on derived quantities does not re-read the CDFs."""
function cached(f, name, deps)
    key = string(hash(read.(deps)); base = 16)
    path = joinpath(mkpath(CACHE), "$name-$key.jls")
    isfile(path) && return deserialize(path)
    x = f()
    foreach(p -> startswith(p, "$name-") && rm(joinpath(CACHE, p)), readdir(CACHE))
    serialize(path, x)
    x
end

"Events on `n` days spread evenly over the mission."
function sample_days(ev, n)
    days = sort!(unique(Date.(ev.start_time)))
    keep = Set(days[round.(Int, range(1, length(days), length = min(n, length(days))))])
    ev[in.(Date.(ev.start_time), Ref(keep)), :]
end

"Commit of this repo the table is built from (`-dirty`: uncommitted changes), so downstream results can name their input."
pipeline_version() = try
    readchomp(`git -C $(@__DIR__) describe --always --dirty --abbrev=10`)
catch
    "unknown"
end

const FULL_TABLE = joinpath(@__DIR__, "..", "data", "tds_tha.csv")
const PULSE_TABLE = joinpath(@__DIR__, "..", "data", "tds_tha_pulses.parquet")

function build(; out = FULL_TABLE, sample = nothing)
    ev = read_events()
    isnothing(sample) || (ev = sample_days(ev, sample))
    # OMNI is one whole-mission stage, so a sample bypasses its cache
    om = isnothing(sample) ? cached(() -> omni_params(ev), "omni", [joinpath(@__DIR__, "..", "data", "raw", "AllData.zip");
        joinpath.(@__DIR__, ["events.jl", "omni.jl"])]) : omni_params(ev)
    ep = event_pulses(ev)
    tw = tds_windows(ev, ep)
    W = efw_windows(ev, tw)
    jobs = [themis_job(ev), efw_job(EFW, W)]
    rth, ref = daily(jobs; files = dayfile)
    th, ef = themis_table(jobs[1], rth), efw_table(W, jobs[2], ref)
    df = leftjoin(leftjoin(ev, th; on = :id), om; on = :id)
    df = leftjoin(df, ef; on = :id)
    isnothing(tw) || (df = leftjoin(df, tw; on = :id))
    sort!(df, :id)
    derive!(df)
    named = reduce(vcat, first.(COLUMNS))
    cols = [named; setdiff(propertynames(df), named, [:source, :lang_Vx, :lang_Vy, :lang_Vz])]
    df = select(df, intersect(cols, propertynames(df)))
    CSV.write(out, df)
    src_of(c) = all(in(LANG_COLUMNS), c) ? "Lang" : "ours"
    if isnothing(sample)
        CSV.write(joinpath(dirname(out), "columns.csv"), DataFrame(name = [join(c, ", ") for (c, _, _) in COLUMNS],
            unit = getindex.(COLUMNS, 2), source = [src_of(c) for (c, _, _) in COLUMNS], description = last.(COLUMNS)))
        write_parquet(replace(out, r"\.csv$" => ".parquet"), df, COLUMNS, src_of)
        tp = tds_pulses(ep)
        if !isnothing(tp)
            CSV.write(joinpath(dirname(out), "pulse_columns.csv"), DataFrame(name = [join(c, ", ") for (c, _, _) in PULSE_COLUMNS],
                unit = getindex.(PULSE_COLUMNS, 2), description = last.(PULSE_COLUMNS)))
            write_parquet(PULSE_TABLE, select(tp, reduce(vcat, first.(PULSE_COLUMNS))), PULSE_COLUMNS, _ -> "ours")
        end
    end
    df
end

"Typed table for pandas/polars; the column dictionary rides along as file metadata \"columns\" (JSON)."
function write_parquet(path, df, columns, src_of)
    unit_of(c, u, k) = (us = split(u, ", "); length(us) == length(c) ? us[k] : u)   # one unit per column of the group, or one for all
    dict = Dict(string(n) => Dict("unit" => unit_of(c, u, k), "source" => src_of(c), "description" => d)
        for (c, u, d) in columns for (k, n) in enumerate(c))
    Parquet2.writefile(path, df; metadata = Dict("columns" => JSON.json(dict), "pipeline" => pipeline_version()))
end

"Rebuild `n` days and print how they differ from the full table."
function check_sample(n = 24; rtol = 1e-6, atol = 1e-6)
    out = joinpath(@__DIR__, "..", "data", "cache", "sample_tds_tha.csv")
    mkpath(dirname(out))
    @time build(; out, sample = n)
    compare_tables(CSV.read(FULL_TABLE, DataFrame), CSV.read(out, DataFrame); rtol, atol)
end

if abspath(PROGRAM_FILE) == @__FILE__()
    arg = findfirst(startswith("--sample"), ARGS)
    isnothing(arg) ? build() : show(check_sample(parse(Int, something(match(r"=(\d+)", ARGS[arg]), (; captures = ["24"])).captures[1])); allrows = true, allcols = true)
end
