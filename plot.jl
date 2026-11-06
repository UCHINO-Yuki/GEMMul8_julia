import Pkg

let required_packages = ["CSV", "DataFrames", "Plots", "LaTeXStrings"]
    missing_packages = filter(pkg -> isnothing(Base.find_package(pkg)), required_packages)
    isempty(missing_packages) || Pkg.add(missing_packages)
end

using CSV
using DataFrames
using Plots
using LaTeXStrings

error_label(stat, method) =
    string(latexstring("\\mathrm{err}_{\\mathrm{", stat, "}}"), " of ", method)

function series_specs(suffix, prec)
    return [
        (func="$(suffix)GEMM-err", label_max=error_label("max", "native FP$prec $(suffix)GEMM"),
            label_min=nothing, color="#333333", marker=:none, linestyle=:dot),
        (func="$(suffix)GEMM-est", label_max=latexstring("ku_{", prec, "}|A||B|"),
            label_min=nothing, color="#777777", marker=:none, linestyle=:dashdot),
        (func="OS2-$suffix-Prev-err", label_max=error_label("max", "Ozaki-II-prev"),
            label_min=error_label("min", "Ozaki-II-prev"), color="#009E73", marker=:square, linestyle=:solid),
        (func="OS2-$suffix-err", label_max=error_label("max", "Ozaki-II"),
            label_min=error_label("min", "Ozaki-II"), color="#0072B2", marker=:circle, linestyle=:solid),
        (func="OS2-$suffix-est", label_max="max. of (4.34)",
            label_min="min. of (4.34)", color="#7B3294", marker=:utriangle, linestyle=:solid),
        (func="OS2-$suffix-est-rough", label_max="max. of (4.35)",
            label_min="min. of (4.35)", color="#B8860B", marker=:dtriangle, linestyle=:solid),
        (func="OS2-$suffix-ef-err", label_max=error_label("max", "Ozaki-II (exact CRT)"),
            label_min=error_label("min", "Ozaki-II (exact CRT)"), color="#D55E00", marker=:diamond, linestyle=:solid),
        (func="OS2-$suffix-ef-est", label_max="max. of (4.36)",
            label_min="min. of (4.36)", color="#0099B4", marker=:hexagon, linestyle=:solid),
    ]
end

function validate_results(df)
    required = [:func, :m, :n, :k, :phi, :num_moduli, :err_max, :err_min]
    all(in(propertynames(df)), required) ||
        error("results.csv must contain $(join(required, ", "))")
    isempty(df) && error("results.csv is empty")
    nrow(unique(df[:, [:m, :n, :k]])) == 1 ||
        error("Plot one matrix shape (m, n, k) at a time")
    for g in groupby(df, [:func, :phi, :num_moduli])
        if nrow(g) != 1
            error("Duplicate series: func=$(g.func[1]), phi=$(g.phi[1]), " *
                "N=$(g.num_moduli[1]). Store exact-CRT estimates as " *
                "OS2-D-ef-est or OS2-S-ef-est, separately from OS2-*-est.")
        end
    end
    return nothing
end

function positive_errors(values)
    errors = Float64.(values)
    errors[.!isfinite.(errors) .| (errors .<= 0)] .= NaN
    return errors
end

function error_limits(g, specs)
    errors = Float64[]
    for spec in specs
        rows = g[g.func .== spec.func, :]
        append!(errors, positive_errors(rows.err_max))
        spec.label_min === nothing || append!(errors, positive_errors(rows.err_min))
    end
    filter!(isfinite, errors)
    isempty(errors) && error("No positive finite values for phi=$(g.phi[1])")
    lo = log10(minimum(errors)) - 0.4
    hi = log10(maximum(errors)) + 0.4
    step = max(1, ceil(Int, (hi - lo) / 6))
    exponents = (ceil(Int, lo/step)*step):step:floor(Int, hi)
    return (10.0^lo, 10.0^hi), 10.0 .^ exponents
end

function validate_series(g, specs)
    moduli = sort(unique(g.num_moduli[g.num_moduli .> 0]))
    isempty(moduli) && error("No Ozaki-II data for phi=$(g.phi[1])")
    for spec in specs
        rows = g[g.func .== spec.func, :]
        isempty(rows) && error("Missing $(spec.func) data for phi=$(g.phi[1])")
        if spec.label_min === nothing
            nrow(rows) == 1 || error("Expected one $(spec.func) row for phi=$(g.phi[1])")
        else
            sort(rows.num_moduli) == moduli ||
                error("Missing or invalid moduli for $(spec.func), phi=$(g.phi[1])")
        end
    end
    return moduli
end

function series_style(spec, column)
    spec.label_min === nothing && return spec.linestyle
    return column === :err_max ? :solid : :dash
end

function plot_error(g, specs; ylim, yticks, moduli, ylabel)
    phi = g.phi[1]
    tick_positions = unique(vcat(moduli[2:4:end], moduli[end]))
    plt = plot(
        yscale=:log10,
        title=latexstring("\\phi = ", phi),
        xlabel="Number of moduli",
        ylabel=ylabel,
        ylim=ylim, yticks=yticks,
        xlim=(first(moduli) - 0.5, last(moduli) + 0.5), xticks=tick_positions,
        legend=false, tickfontsize=13, guidefontsize=15, titlefontsize=16,
        gridalpha=0.15, framestyle=:box,
        left_margin=8Plots.mm, right_margin=4Plots.mm,
        top_margin=2Plots.mm, bottom_margin=6Plots.mm,
    )
    omitted_nonfinite = false
    for spec in specs
        rows = sort(g[g.func .== spec.func, :], :num_moduli)
        for (column, label) in ((:err_max, spec.label_max), (:err_min, spec.label_min))
            label === nothing && continue
            raw = Float64.(rows[!, column])
            omitted_nonfinite |= any(x -> !isfinite(x), raw)
            values = positive_errors(raw)
            linestyle = series_style(spec, column)
            if spec.label_min === nothing
                hline!(plt, values; label=label, color=spec.color,
                    linestyle=linestyle, linewidth=2)
            else
                plot!(plt, rows.num_moduli, values; label=label, color=spec.color,
                    linewidth=2, linestyle=linestyle)
                if spec.marker !== :none
                    scatter!(plt, rows.num_moduli, values;
                        label="", marker=spec.marker, markersize=5,
                        markercolor=column === :err_max ? spec.color : :white,
                        markerstrokecolor=spec.color, markerstrokewidth=1)
                end
            end
        end
    end
    if omitted_nonfinite
        y = 10.0^(log10(ylim[2]) - 0.06 * log10(ylim[2] / ylim[1]))
        annotate!(plt, last(moduli) - 0.5, y,
            text("Nonfinite values omitted", 9, :right, :gray))
    end
    return plt
end

function legend_panel(specs)
    nlabels = sum(s -> s.label_min === nothing ? 1 : 2, specs)
    plt = plot(; axis=false, grid=false, framestyle=:none,
        xlim=(0, 1), ylim=(0, 1), legend=:top, legendfontsize=14,
        legend_column=min(4, cld(nlabels, 2)),
        foreground_color_legend=:transparent, background_color_legend=:transparent,
        left_margin=0Plots.mm, right_margin=0Plots.mm,
        top_margin=0Plots.mm, bottom_margin=0Plots.mm)
    ordered_specs = vcat(filter(s -> s.label_min !== nothing, specs),
        filter(s -> s.label_min === nothing, specs))
    for spec in ordered_specs
        for (column, label) in ((:err_max, spec.label_max), (:err_min, spec.label_min))
            label === nothing && continue
            plot!(plt, [NaN], [NaN]; label=label, color=spec.color,
                linestyle=series_style(spec, column), linewidth=2, marker=spec.marker,
                markersize=5, markercolor=column === :err_max ? spec.color : :white,
                markerstrokecolor=spec.color, markerstrokewidth=1)
        end
    end
    return plt
end

function make_figure(df, suffix, prec, kind::Symbol)
    kind in (:error, :bound) || throw(ArgumentError("kind must be :error or :bound"))
    specs = series_specs(suffix, prec)
    filter!(s -> endswith(s.func, "-err") == (kind === :error), specs)
    if suffix == "S" && !any(df.func .== "OS2-S-Prev-err")
        filter!(s -> s.func != "OS2-S-Prev-err", specs)
    end
    data = df[in.(df.func, Ref([s.func for s in specs])), :]
    isempty(data) && error("No $(suffix)GEMM data")
    groups = collect(groupby(data, :phi; sort=true))
    panels = Plots.Plot[]
    for g in groups
        moduli = validate_series(g, specs)
        limits, ticks = error_limits(g, specs)
        ylabel = kind === :error ? "Absolute error" : "Absolute error bound"
        push!(panels, plot_error(g, specs;
            ylim=limits, yticks=ticks, moduli=moduli, ylabel=ylabel))
    end
    ncols = min(3, length(groups))
    nrows = cld(length(groups), ncols)
    while length(panels) < nrows * ncols
        push!(panels, plot(; axis=false, grid=false, framestyle=:none))
    end
    push!(panels, legend_panel(specs))
    figure_layout = @layout [grid(nrows, ncols); legend{0.05h}]
    return plot(panels...;
        layout=figure_layout,
        size=(560 * ncols, 440 * nrows + 50), dpi=300, overwrite_figure=false)
end

df = CSV.read(joinpath(@__DIR__, "results.csv"), DataFrame)
validate_results(df)
for (suffix, prec) in (("S", "32"), ("D", "64"))
    for kind in (:error, :bound)
        plt = make_figure(df, suffix, prec, kind)
        output = joinpath(@__DIR__, "plot_$(suffix)_$(kind).png")
        savefig(plt, output)
        println("Saved: ", output)
    end
end
