import Pkg

let required_packages = ["CSV", "DataFrames", "Plots", "LaTeXStrings"]
    missing_packages = filter(pkg -> isnothing(Base.find_package(pkg)), required_packages)
    if !isempty(missing_packages)
        Pkg.add(missing_packages)
    end
end

using CSV
using DataFrames
using Plots
using LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "results.csv"), DataFrame)
isempty(df) && error("results.csv is empty")

df.max_rel_err = Float64.(df.max_rel_err)
invalid = .!isfinite.(df.max_rel_err) .| (df.max_rel_err .<= 0)
if any(invalid)
    @warn "Omitting nonpositive or nonfinite errors from the log plots" count=count(invalid)
    df.max_rel_err[invalid] .= NaN
end
valid_errors = filter(isfinite, df.max_rel_err)
isempty(valid_errors) && error("No positive finite errors to plot")

ymin = minimum(valid_errors) * 0.5
ymax = maximum(valid_errors) * 2
yexponents = ceil(Int, log10(ymin)):2:floor(Int, log10(ymax))
ytick_values = 10.0 .^ yexponents

function plot_precision(g, suffix, crt_label; ylim=(ymin, ymax))
    phi = g.phi[1]
    oz2 = sort(g[g.func .== "OS2-$suffix", :], :num_moduli)
    oz2_CRT_exact = sort(g[g.func .== "OS2-$suffix-CRT-exact", :], :num_moduli)
    GEMM = g[g.func .== "$(suffix)GEMM", :]

    isempty(oz2) && error("Missing OS2-$suffix data for phi = $phi")
    isempty(oz2_CRT_exact) && error("Missing OS2-$suffix-CRT-exact data for phi = $phi")
    nrow(GEMM) == 1 || error("Expected one $(suffix)GEMM row for phi = $phi")

    moduli = sort(unique(vcat(oz2.num_moduli, oz2_CRT_exact.num_moduli)))
    tick_positions = unique(vcat(moduli[2:4:end], moduli[end]))
    plt = plot(
        yscale=:log10,
        title=latexstring("\\phi = ", phi),
        xlabel="Number of moduli",
        ylabel="Maximum relative error",
        ylim=ylim,
        yticks=ytick_values,
        xlim=(first(moduli) - 0.5, last(moduli) + 0.5),
        xticks=tick_positions,
        legend=:bottomleft,
        legendfontsize=9,
        tickfontsize=9,
        guidefontsize=11,
        titlefontsize=12,
        gridalpha=0.2,
        left_margin=10Plots.mm,
        right_margin=4Plots.mm,
        top_margin=0Plots.mm,
        bottom_margin=6Plots.mm,
    )
    plot!(plt, oz2.num_moduli, oz2.max_rel_err;
        label="Ozaki-II ($crt_label)", color=:royalblue,
        linewidth=2, marker=:circle, markersize=3)
    plot!(plt, oz2_CRT_exact.num_moduli, oz2_CRT_exact.max_rel_err;
        label="Ozaki-II (Exact CRT)", color=:darkorange,
        linewidth=2, linestyle=:dash, marker=:diamond, markersize=3)
    hline!(plt, [GEMM.max_rel_err[1]];
        label="native $(suffix)GEMM", color=:black, linestyle=:dot, linewidth=2)
    return plt
end

plot_S(g; kwargs...) = plot_precision(g, "S", "FP64 CRT"; kwargs...)
plot_D(g; kwargs...) = plot_precision(g, "D", "Pseudo-DD CRT"; kwargs...)

groups = collect(groupby(df, :phi; sort=true))
ncols = min(3, length(groups))
nrows = cld(length(groups), ncols)
for (plotter, filename) in ((plot_S, "plot_S.png"), (plot_D, "plot_D.png"))
    plt = plot([plotter(g) for g in groups]...;
        layout=(nrows, ncols), size=(500 * ncols, 400 * nrows))
    for i in eachindex(groups)
        plot!(plt; subplot=i, legend=(i == length(groups) ? :bottomleft : false))
    end
    output = joinpath(@__DIR__, filename)
    savefig(plt, output)
    println("Saved: ", output)
end
