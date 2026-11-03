using CairoMakie

"""
Waveform panel for burst `b` into grid position `gl`: high-passed E∥ and Ez over the burst (Lang windows
shaded, accepted spikes red ticks, rejected candidates grey) and a ±`zw` s zoom (Ex, Ey, Ez, E∥) at `tz`.
"""
function burst_panel!(gl, b::Burst, sp, lw; tz = nothing, zw = 0.02, title = "", p = PARAMS)
    E = copy(b.E); for c in 1:3; fillnan!(view(E, c, :)); end
    hp = [first(highpass_chunks(E[c, :], b.t, b.fs, p.τ)) for c in 1:3]
    par = any(isnan, b.bh) ? fill(NaN, length(b.t)) : vec(sum(reduce(hcat, hp)' .* b.bh, dims = 1))
    tt = b.t .- b.t[1]
    axs = [Axis(gl[i, 1], ylabel = l, xticklabelsvisible = i == 2, titlesize = 11, ylabelsize = 10,
               title = i == 1 ? title : "") for (i, l) in enumerate(["E∥ hp", "Ez hp"])]
    for (ax, y) in zip(axs, (par, hp[3]))
        for (a, c) in lw; vspan!(ax, a - b.t[1], c - b.t[1], color = (:orange, 0.2)); end
        lines!(ax, tt, y, color = :black, linewidth = 0.3)
        r = sp[.!sp.tdslike, :]; isempty(r) || scatter!(ax, r.t .- b.t[1], fill(maximum(abs, y), nrow(r)), color = (:gray, 0.5), marker = :vline, markersize = 6)
        r = sp[sp.tdslike, :]; isempty(r) || scatter!(ax, r.t .- b.t[1], fill(maximum(abs, y), nrow(r)), color = :red, marker = :vline, markersize = 8)
    end
    linkxaxes!(axs...)
    k = tz === nothing ? argmax(abs.(isnan(par[1]) ? hp[3] : par) .* .!masked.(b.φ)) : searchsortedfirst(b.t, tz)
    r = max(1, k - round(Int, zw * b.fs)):min(length(tt), k + round(Int, zw * b.fs))
    az = Axis(gl[1:2, 2], xlabel = "ms", titlesize = 10, title = isnan(b.φ[k]) ? "φ unknown" : "φ = $(round(Int, b.φ[k]))°")
    for (y, c, l) in zip((hp[1], hp[2], hp[3], par), (:gray60, :green, :blue, :red), ("Ex", "Ey", "Ez", "E∥"))
        lines!(az, 1000 .* (tt[r] .- tt[k]), y[r], color = c, linewidth = 0.8, label = l)
    end
    colsize!(gl, 2, Relative(0.3))
    az
end

lang_windows(b, ev) = [(datetime2unix(r.start_time), datetime2unix(r.end_time) + 1) for r in eachrow(ev)
                       if datetime2unix(r.start_time) <= b.t[end] && datetime2unix(r.end_time) + 1 >= b.t[1]]

"Stack of burst panels for (date, burst, title) rows; spikes from `S`."
function gallery(rows, S, ev, file; getday = load_day)
    fig = Figure(size = (1200, 280 * length(rows)), fontsize = 11)
    for (i, (d, k, title)) in enumerate(rows)
        b = getday(d)[k]
        sp = S[(S.date .== d) .& (S.burst .== k), :]
        tz = any(sp.tdslike) ? sp.t[sp.tdslike][argmax(abs.(sp.a[sp.tdslike]))] : nothing
        az = burst_panel!(GridLayout(fig[i, 1]), b, sp, lang_windows(b, ev); tz, title = "$d #$k  " * title)
        i == 1 && axislegend(az, position = :lt, labelsize = 8, patchsize = (8, 4), padding = 2, rowgap = 0)
    end
    save(file, fig)
end

"Grid of ±zw s zooms (Ex, Ey, Ez, E∥ high-passed) around pulses `rows` (date, burst, t, title)."
function pulse_grid(rows, file; zw = 0.012, ncol = 4, p = PARAMS)
    fig = Figure(size = (300ncol, 220cld(length(rows), ncol)), fontsize = 9)
    cache = Dict{Date, Vector{Burst}}()
    for (n, (d, k, t, title)) in enumerate(rows)
        b = get!(() -> load_day(d), cache, d)[k]
        E = copy(b.E); for c in 1:3; fillnan!(view(E, c, :)); end
        hp = [first(highpass_chunks(E[c, :], b.t, b.fs, p.τ)) for c in 1:3]
        i = searchsortedfirst(b.t, t); r = max(1, i - round(Int, zw * b.fs)):min(length(b.t), i + round(Int, zw * b.fs))
        ax = Axis(fig[(n - 1) ÷ ncol + 1, (n - 1) % ncol + 1], title = title, titlesize = 8)
        par = any(isnan, b.bh) ? nothing : sum(hp[c][r] .* b.bh[c, r] for c in 1:3)
        for (y, c) in zip((hp[1][r], hp[2][r], hp[3][r]), (:gray60, :green, :blue)); lines!(ax, 1000 .* (b.t[r] .- t), y, color = c, linewidth = 0.7); end
        par === nothing || lines!(ax, 1000 .* (b.t[r] .- t), par, color = :red, linewidth = 0.9)
    end
    save(file, fig)
end
