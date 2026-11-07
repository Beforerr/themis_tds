# data/dates.txt for scripts/download.sh: event days plus the neighbouring days reached by windows around events
# (read through `spanvars`; the widest is bgwindow.jl's ±20 min).
using Dates
include("events.jl")

const CONTEXT = Minute(20)

ev = read_events()
days = sort!(unique([Date.(ev.start_time .- CONTEXT); Date.(ev.start_time); Date.(ev.end_time .+ CONTEXT)]))
write(joinpath(@__DIR__, "..", "data", "dates.txt"), join(Dates.format.(days, "yyyymmdd"), "\n") * "\n")
println(length(days), " days")
