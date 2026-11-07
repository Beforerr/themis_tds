using Dates

const CDFDIR = joinpath(@__DIR__, "..", "data", "cdf")

"Latest version of a THEMIS day file of product `prod` (e.g. \"fgm\" or :fgm), or `nothing`; also the `ProbeDay` resolver."
function dayfile(prod, date, probe = "tha")
    prod = string(prod)
    dir = joinpath(CDFDIR, prod)
    pre = "$(probe)_l$(prod == "efw" ? 1 : 2)_$(prod)_$(Dates.format(date, "yyyymmdd"))_v"
    fs = readdir(dir)
    i = findlast(f -> startswith(f, pre) && endswith(f, ".cdf"), fs)
    isnothing(i) ? nothing : joinpath(dir, fs[i])
end
