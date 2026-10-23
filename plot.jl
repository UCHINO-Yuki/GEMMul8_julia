using CSV
using DataFrames
using Plots

df = CSV.read(joinpath(@__DIR__, "results.csv"), DataFrame)

ymin = minimum(df.max_rel_err) * 0.5
ymax = maximum(df.max_rel_err) * 2

function plot_S(g; ylim=(ymin, ymax))
    phi = g.phi[1]
    plt = plot(yscale=:log, title="phi = $phi", ylim=ylim)
    oz2 = g[g.func .== "OS2-S", :]
    oz2_CRT_exact = g[g.func .== "OS2-S-CRT-exact", :]
    GEMM = g[g.func .== "SGEMM", :]
    xticks = sort(unique(oz2.num_moduli))
    plot!(plt, oz2.num_moduli, oz2.max_rel_err, label="Ozaki-II (FP64 CRT)")
    plot!(plt, oz2_CRT_exact.num_moduli, oz2_CRT_exact.max_rel_err, label="Ozaki-II (Exact CRT)")
    plot!(plt, xticks, [GEMM.max_rel_err[1] for _ in xticks], label="native SGEMM", xticks=xticks)
    plt
end

for g in groupby(df, :phi)
    plot_S(g)
end

plt = plot([plot_S(g, ylim=(ymin, ymax)) for g in groupby(df, :phi)]..., layout = (1, length(unique(df.phi))), size = (1400, 400), yticks = (10.0 .^ (0:-2:-16), ["10^{$e}" for e in 0:-2:-16]))
savefig(plt, joinpath(@__DIR__, "plot_S.png"))


function plot_D(g; ylim=(ymin, ymax))
    phi = g.phi[1]
    plt = plot(yscale=:log, title="phi = $phi", ylim=ylim)
    oz2 = g[g.func .== "OS2-D", :]
    oz2_CRT_exact = g[g.func .== "OS2-D-CRT-exact", :]
    GEMM = g[g.func .== "DGEMM", :]
    xticks = sort(unique(oz2.num_moduli))
    plot!(plt, oz2.num_moduli, oz2.max_rel_err, label="Ozaki-II (Pseudo-DD CRT)")
    plot!(plt, oz2_CRT_exact.num_moduli, oz2_CRT_exact.max_rel_err, label="Ozaki-II (Exact CRT)")
    plot!(plt, xticks, [GEMM.max_rel_err[1] for _ in xticks], label="native DGEMM", xticks=xticks)
    plt
end

for g in groupby(df, :phi)
    plot_D(g)
end

plt = plot([plot_D(g, ylim=(ymin, ymax)) for g in groupby(df, :phi)]..., layout = (1, length(unique(df.phi))), size = (1400, 400), yticks = (10.0 .^ (0:-2:-16), ["10^{$e}" for e in 0:-2:-16]))
savefig(plt, joinpath(@__DIR__, "plot_D.png"))
