using DataFrames

"""
    compare_tables(a, b; key = :id, rtol = 1e-9, atol = 1e-9, skip = Symbol[]) -> DataFrame

Differences of table `b` from `a` over the rows they share by `key`: one row per differing column, with counts,
largest absolute and relative differences and the first differing key. NaN equals NaN and missing equals missing.
"""
function compare_tables(a, b; key = :id, rtol = 1e-9, atol = 1e-9, skip = Symbol[])
    only_a, only_b = setdiff(propertynames(a), propertynames(b)), setdiff(propertynames(b), propertynames(a))
    ia = Dict(k => i for (i, k) in enumerate(a[!, key]))
    shared = [(ia[k], j) for (j, k) in enumerate(b[!, key]) if haskey(ia, k)]
    rows = NamedTuple[]
    for c in setdiff(intersect(propertynames(a), propertynames(b)), [key], skip)
        n, dabs, drel, first_key = 0, 0.0, 0.0, nothing
        for (i, j) in shared
            x, y = a[i, c], b[j, c]
            if x isa Real && y isa Real && !(x isa Bool || y isa Bool)
                (isnan(x) && isnan(y)) && continue
                d = abs(float(x) - float(y))
                !isnan(d) && d <= max(atol, rtol * max(abs(float(x)), abs(float(y)))) && continue
                isnan(d) || (dabs = max(dabs, d); drel = max(drel, d / max(abs(float(x)), abs(float(y)))))
            elseif isequal(x, y)
                continue
            end
            n += 1
            isnothing(first_key) && (first_key = b[j, key])
        end
        n > 0 && push!(rows, (; column = c, ndiff = n, max_abs = dabs, max_rel = drel, first_key))
    end
    println("compared ", length(shared), " rows and ", ncol(a) - length(only_a), " columns (", nrow(b) - length(shared), " rows only in b): ", length(rows), " columns differ")
    isempty(only_a) || println("only in a: ", join(only_a, ", "))
    isempty(only_b) || println("only in b: ", join(only_b, ", "))
    isempty(rows) ? DataFrame(column = Symbol[], ndiff = Int[], max_abs = Float64[], max_rel = Float64[], first_key = Any[]) : DataFrame(rows)
end
