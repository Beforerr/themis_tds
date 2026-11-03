using CSV, DataFrames, Dates

const EVENT_COLS = [:start_time, :end_time, :duration, :radius, :gsm_x, :gsm_y, :gsm_z, :mlt, :mlat,
    :lang_max_Ez, :lang_rms_Ew, :lang_Bx, :lang_By, :lang_Bz, :lang_Vx, :lang_Vy, :lang_Vz]

const TFMT = dateformat"yyyy-mm-dd/HH:MM:SS"

"""Read Lang Qin's per-day event files. `*-QLP.txt` re-processes the same day (with full V vectors) and replaces the plain file."""
function read_events(root = joinpath(@__DIR__, "..", "data", "raw", "AllData"))
    files = String[]
    for (dir, _, fs) in walkdir(root), f in fs
        endswith(f, ".txt") && push!(files, joinpath(dir, f))
    end
    qlp = Set(replace(f, "-QLP" => "") for f in files if occursin("-QLP", f))
    filter!(f -> !(f in qlp), files)
    dfs = map(sort(files)) do f
        df = CSV.read(f, DataFrame; header = EVENT_COLS, delim = ' ', ignorerepeated = true,
            types = Dict(:start_time => String, :end_time => String), on_error = :collect)  # blank trailer rows; dropped below
        df.source .= basename(f)
        df
    end
    ev = dropmissing(reduce(vcat, dfs), :start_time)
    ev.start_time = DateTime.(ev.start_time, TFMT)
    ev.end_time = DateTime.(ev.end_time, TFMT)
    sort!(ev, :start_time)
    unique!(ev, [:start_time, :end_time])
    ev.id = 1:nrow(ev)
    ev
end
