import Pkg

let required_packages = ["MultiFloats", "CSV", "DataFrames"]
    missing_packages = filter(pkg -> isnothing(Base.find_package(pkg)), required_packages)
    if !isempty(missing_packages)
        Pkg.add(missing_packages)
    end
end

using LinearAlgebra
using MultiFloats
using Random
using CSV, DataFrames, Dates
using InteractiveUtils

include(joinpath(@__DIR__, "gemmul8.jl"))

gen_rand(m::Int64, n::Int64) = rand(Float64, m, n)
gen_randn(m::Int64, n::Int64) = randn(Float64, m, n)
gen_mat(m::Int64, n::Int64, phi::AbstractFloat) = (gen_rand(m, n) .- 0.5) .* exp.(gen_randn(m, n) .* phi)

function calc_err(C::AbstractMatrix{<:AbstractFloat}, C_ref::AbstractMatrix{Float64x2})
    abserr = abs.(Float64x2.(C) .- Float64x2.(C_ref))
    relerr = abserr ./ Float64x2.(abs.(C_ref))
    Float64(maximum(relerr))
end

println("-----")
versioninfo()
println("-----")

Random.seed!(1)
m = 8
n = 8
k = 8192

rows = NamedTuple[]
for phi_i in [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]
    phi = Float64(phi_i)

    A = gen_mat(m, k, phi)
    B = gen_mat(k, n, phi)
    C_64 = A * B
    C_dd = Float64x2.(A) * Float64x2.(B)
    err_64 = calc_err(C_64, C_dd)

    push!(rows,
        (func="DGEMM", m=m, n=n, k=k, phi=phi, num_moduli=-1, max_rel_err=err_64),
    )

    for num_moduli in 3:40
        C_emu = os2(A, B, num_moduli)
        err_64_emu = calc_err(C_emu, C_dd)

        push!(rows,
            (func="OS2-D", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, max_rel_err=err_64_emu),
        )
        
        C_emu = oz2_err_free_crt(A, B, num_moduli)
        err_64_err_free = calc_err(C_emu, C_dd)
        push!(rows,
            (func="OS2-D-CRT-exact", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, max_rel_err=err_64_err_free),
        )
    end

    A32 = Float32.(A)
    B32 = Float32.(B)
    C_32 = A32 * B32
    C_dd = Float64x2.(A32) * Float64x2.(B32)
    err_32 = calc_err(C_32, C_dd)

    push!(rows,
        (func="SGEMM", m=m, n=n, k=k, phi=phi, num_moduli=-1, max_rel_err=err_32),
    )

    for num_moduli in 3:40
        C_emu = os2(A32, B32, num_moduli)
        err_32_emu = calc_err(C_emu, C_dd)

        push!(rows,
            (func="OS2-S", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, max_rel_err=err_32_emu),
        )
        
        C_emu = oz2_err_free_crt(A32, B32, num_moduli)
        err_32_err_free = calc_err(C_emu, C_dd)
        push!(rows,
            (func="OS2-S-CRT-exact", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, max_rel_err=err_32_err_free),
        )
    end
end

df = DataFrame(rows)
CSV.write(joinpath(@__DIR__, "results.csv"), df)
