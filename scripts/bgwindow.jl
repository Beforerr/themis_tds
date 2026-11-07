# Background-field window test (after build.jl) → data/bgwindow.csv
# For each event and half-window W, how far b̂_W is from the event-window field, from the adopted ±5-min b̂,
# and between the two halves of the window (a rotation there means the window straddles a boundary or flapping).
using CSV, DataFrames, Dates, LinearAlgebra
include("themis.jl")

const HALF_WINDOWS = [15.0, 30.0, 60.0, 120.0, 300.0, 600.0, 1200.0]

angle(a, b) = acosd(clamp(dot(a, b) / (norm(a) * norm(b)), -1, 1))

function bgwindow(df = CSV.read(joinpath(@__DIR__, "..", "data", "tds_tha.csv"), DataFrame))
    days = collect(groupby(transform(df, :start_time => ByRow(Date) => :date), :date))
    out = Vector{Vector{NamedTuple}}(undef, length(days))
    Threads.@threads for i in eachindex(days)
        g = days[i]; d = ProbeDay("tha", g.date[1], dayfile)
        out[i] = NamedTuple[]
        isnothing(input(d, :fgm)) && continue
        tcs = [(datetime2unix(DateTime(r.start_time)) + datetime2unix(DateTime(r.end_time))) / 2 for r in eachrow(g)]
        Wmax = maximum(HALF_WINDOWS)
        fgs = readseries(d, :fgm, "tha_fgs_time", "tha_fgs_gsm", (minimum(tcs) - Wmax, maximum(tcs) + Wmax))
        for (r, tc) in zip(eachrow(g), tcs)
            B = [r.B_x, r.B_y, r.B_z]; V = [r.Vi_x, r.Vi_y, r.Vi_z]
            ref, = wmed(fgs, tc - BG_HALF, tc + BG_HALF)
            for W in HALF_WINDOWS
                bw, n = wmed(fgs, tc - W, tc + W)
                n == 0 && continue
                pre, = wmed(fgs, tc - W, tc); post, = wmed(fgs, tc, tc + W)
                push!(out[i], (id = r.id, W = W, ang_event = angle(bw, B), ang_ref = angle(bw, ref),
                    ang_halves = angle(pre, post), V_par = dot(V, bw) / norm(bw), region = r.region))
            end
        end
    end
    DataFrame(reduce(vcat, out))
end

abspath(PROGRAM_FILE) == @__FILE__() && CSV.write(joinpath(@__DIR__, "..", "data", "bgwindow.csv"), bgwindow())
