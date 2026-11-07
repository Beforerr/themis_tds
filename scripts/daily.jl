# Per-(probe, day) products: a product is a function of one probe-day's files and an argument (e.g. that day's events).
# `daily` runs every stale product of a day in one task over that day's files and caches each result under
# data/cache/<name>/<probe>/<yyyy>/, keyed by this file, the product's source files, the day's input file names and the argument.
# Input files are versioned and never rewritten, so their names stand for their content. Package versions are not
# part of the key: clear data/cache/<name> after an upgrade that changes results.
using CDFDatasets, DataFrames, Dates, Serialization

const DAILY_CACHE = joinpath(@__DIR__, "..", "data", "cache")
const REPO = joinpath(@__DIR__, "..")

"""
One probe-day. `files(kind, date, probe)` resolves a file kind (e.g. `:fgm`) to a path or `nothing`. Whole variables
read with `dayvar` are decoded once and shared by all products of the day, so they must not be mutated.
"""
struct ProbeDay
    probe::String
    date::Date
    files::Function
    paths::Dict{Symbol, Union{Nothing, String}}
    sets::Dict{Symbol, Any}
    vars::Dict{Tuple{Symbol, String}, Any}
    nb::Dict{Int, ProbeDay}
end
ProbeDay(probe, date, files) = ProbeDay(probe, date, files, Dict(), Dict(), Dict(), Dict())
neighbour(d::ProbeDay, k) = get!(() -> ProbeDay(d.probe, d.date + Day(k), d.files), d.nb, k)

input(d::ProbeDay, kind) = get!(() -> d.files(kind, d.date, d.probe), d.paths, kind)
"Lazy dataset of a file kind (for partial reads), or `nothing`."
dayset(d::ProbeDay, kind) = get!(d.sets, kind) do
    p = input(d, kind)
    isnothing(p) ? nothing : CDFDataset(p)
end
"Whole decoded variable, or `nothing` if the file or variable is absent or empty."
dayvar(d::ProbeDay, kind, name) = get!(d.vars, (kind, name)) do
    ds = dayset(d, kind)
    isnothing(ds) || !haskey(ds, name) || length(ds[name]) == 0 ? nothing : Array(ds[name])
end

"""
    spanvars(d, kind, tvar, vars, lo, hi) -> (t, x1, x2, ...) or `nothing`

`dayvar`s joined with the neighbouring days' records in [lo, hi] (units of `tvar`; records along the last dimension)
that lie outside the day file's own time span, for windows crossing a file boundary (THEMIS day files run from minutes
after midnight to minutes after the next). The product must set `neighbours = true` so those days' files key its cache.
"""
function spanvars(d::ProbeDay, kind, tvar, vars, lo, hi)
    t = dayvar(d, kind, tvar); xs = map(v -> dayvar(d, kind, v), vars)
    (isnothing(t) || any(isnothing, xs)) && return nothing
    a, b = extrema(t)
    parts = [(t, xs)]
    for (k, keep) in ((-1, x -> lo <= x < a), (1, x -> b < x <= hi))
        (k < 0 ? lo < a : hi > b) || continue
        n = neighbour(d, k); tn = dayvar(n, kind, tvar); xn = map(v -> dayvar(n, kind, v), vars)
        (isnothing(tn) || any(isnothing, xn)) && continue
        r = findall(keep, tn)
        isempty(r) && continue
        part = (tn[r], map(x -> x[ntuple(_ -> Colon(), ndims(x) - 1)..., r], xn))
        k < 0 ? pushfirst!(parts, part) : push!(parts, part)
    end
    length(parts) == 1 && return (t, xs...)
    (reduce(vcat, first.(parts)), (reduce((a, b) -> cat(a, b; dims = ndims(a)), getindex.(last.(parts), i)) for i in eachindex(vars))...)
end

"""`f(day, arg)` → result (`nothing` is a valid result); `code`: repo files the result depends on; `inputs`: file kinds read;
`neighbours`: whether `f` reads the adjacent days through `spanvars`."""
struct DailyProduct
    name::String
    code::Vector{String}
    inputs::Vector{Symbol}
    f::Function
    neighbours::Bool
end
DailyProduct(name, code, inputs, f; neighbours = false) = DailyProduct(name, code, inputs, f, neighbours)

"`product => items` with one item `(probe, date, arg)` per day of `dates`; `arg`: that day's rows of `df` as a NamedTuple of vectors."
function perday(p::DailyProduct, df, dates; probe = "tha")
    g = groupby(insertcols(df, :_date => dates), :_date; sort = true)
    p => [(probe, k._date, (; (c => collect(s[!, c]) for c in propertynames(df))...)) for (k, s) in pairs(g)]
end

"""
    daily(jobs; files) -> one result vector per job, aligned with its items

`jobs`: `product => items` (see `perday`). Threaded over probe-days; a product that throws is logged, gives `nothing`
and is not cached.
"""
function daily(jobs; files)
    codekeys = [hash((read(@__FILE__), read.(joinpath.(REPO, p.code)))) for (p, _) in jobs]   # this file: read semantics
    units = Dict{Tuple{String, Date}, Vector{Tuple{Int, Int}}}()
    for (j, (_, items)) in enumerate(jobs), (i, (probe, date, _)) in enumerate(items)
        push!(get!(units, (probe, date), Tuple{Int, Int}[]), (j, i))
    end
    res = [Vector{Any}(nothing, length(items)) for (_, items) in jobs]
    ran = [Threads.Atomic{Int}(0) for _ in jobs]
    Threads.@threads :greedy for (u, refs) in collect(units)
        d = ProbeDay(u..., files)
        for (j, i) in refs
            p, items = jobs[j]
            res[j][i] = cached_day(p, d, items[i][3], codekeys[j], ran[j])
        end
    end
    for (j, (p, items)) in enumerate(jobs)
        @info "daily $(p.name)" days = length(items) computed = ran[j][]
    end
    res
end

function cached_day(p, d, arg, codekey, ran)
    io = IOBuffer(); serialize(io, arg)
    days = p.neighbours ? (neighbour(d, -1), d, neighbour(d, 1)) : (d,)
    names = [(f = input(x, k); isnothing(f) ? "" : basename(f)) for k in p.inputs for x in days]
    dir = joinpath(DAILY_CACHE, p.name, d.probe, string(year(d.date)))
    stem = Dates.format(d.date, "yyyymmdd")
    path = joinpath(dir, "$stem-$(string(hash((codekey, names, take!(io))); base = 16)).jls")
    isfile(path) && return deserialize(path)
    r = try
        p.f(d, arg)
    catch e
        @warn "$(p.name) failed" d.probe d.date exception = e
        return nothing
    end
    Threads.atomic_add!(ran, 1)
    mkpath(dir)
    foreach(f -> startswith(f, stem) && rm(joinpath(dir, f)), readdir(dir))
    serialize(path * ".tmp", r); mv(path * ".tmp", path)
    r
end
