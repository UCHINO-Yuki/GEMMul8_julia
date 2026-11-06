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
    err = Float64.(abs.(Float64x2.(C) .- Float64x2.(C_ref)))
    (maximum(err), minimum(err))
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
    err_64_max, err_64_min = calc_err(C_64, C_dd)

    push!(rows,
        (func="DGEMM-err", m=m, n=n, k=k, phi=phi, num_moduli=-1, err_max=err_64_max, err_min=err_64_min),
    )

    C_64_abs = Float64.(Float64x2.(abs.(A)) * Float64x2.(abs.(B)))
    est_64_max = maximum(k .* 2.0^-53 .* C_64_abs)
    push!(rows,
        (func="DGEMM-est", m=m, n=n, k=k, phi=phi, num_moduli=-1, err_max=est_64_max, err_min=est_64_max),
    )

    for num_moduli in 3:40

        bounds_64 = error_bound(A, B, num_moduli)

        C_emu_prev = os2_prev(A, B, num_moduli)
        err_64_emu_prev_max, err_64_emu_prev_min = calc_err(C_emu_prev, C_dd)
        push!(rows,
            (func="OS2-D-Prev-err", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=err_64_emu_prev_max, err_min=err_64_emu_prev_min),
        )

        C_emu = os2(A, B, num_moduli)
        err_64_emu_max, err_64_emu_min = calc_err(C_emu, C_dd)
        push!(rows,
            (func="OS2-D-err", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=err_64_emu_max, err_min=err_64_emu_min),
        )

        est_64_emu_max = maximum(bounds_64.err1)
        est_64_emu_min = minimum(bounds_64.err1)
        push!(rows,
            (func="OS2-D-est", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_64_emu_max, err_min=est_64_emu_min),
        )

        est_64_emu_max = maximum(bounds_64.err2)
        est_64_emu_min = minimum(bounds_64.err2)
        push!(rows,
            (func="OS2-D-est-rough", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_64_emu_max, err_min=est_64_emu_min),
        )

        C_emu_ef = oz2_err_free_crt(A, B, num_moduli)
        err_64_ef_max, err_64_ef_min = calc_err(C_emu_ef, C_dd)
        push!(rows,
            (func="OS2-D-ef-err", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=err_64_ef_max, err_min=err_64_ef_min),
        )

        est_64_ef_max = maximum(bounds_64.err_exact_crt)
        est_64_ef_min = minimum(bounds_64.err_exact_crt)
        push!(rows,
            (func="OS2-D-ef-est", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_64_ef_max, err_min=est_64_ef_min),
        )
    end

    A32 = Float32.(A)
    B32 = Float32.(B)
    C_32 = A32 * B32
    C_dd = Float64x2.(A32) * Float64x2.(B32)
    err_32_max, err_32_min = calc_err(C_32, C_dd)

    push!(rows,
        (func="SGEMM-err", m=m, n=n, k=k, phi=phi, num_moduli=-1, err_max=err_32_max, err_min=err_32_min),
    )

    C_32_abs = Float64.(Float64x2.(abs.(A32)) * Float64x2.(abs.(B32)))
    est_32_max = maximum(k .* 2.0^-24 .* C_32_abs)
    push!(rows,
        (func="SGEMM-est", m=m, n=n, k=k, phi=phi, num_moduli=-1, err_max=est_32_max, err_min=est_32_max),
    )

    for num_moduli in 3:40

        bounds_32 = error_bound(A32, B32, num_moduli)

        C_emu = os2(A32, B32, num_moduli)
        err_32_emu_max, err_32_emu_min = calc_err(C_emu, C_dd)
        push!(rows,
            (func="OS2-S-err", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=err_32_emu_max, err_min=err_32_emu_min),
        )

        est_32_emu_max = maximum(bounds_32.err1)
        est_32_emu_min = minimum(bounds_32.err1)
        push!(rows,
            (func="OS2-S-est", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_32_emu_max, err_min=est_32_emu_min),
        )

        est_32_emu_max = maximum(bounds_32.err2)
        est_32_emu_min = minimum(bounds_32.err2)
        push!(rows,
            (func="OS2-S-est-rough", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_32_emu_max, err_min=est_32_emu_min),
        )

        C_emu_ef = oz2_err_free_crt(A32, B32, num_moduli)
        err_32_ef_max, err_32_ef_min = calc_err(C_emu_ef, C_dd)
        push!(rows,
            (func="OS2-S-ef-err", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=err_32_ef_max, err_min=err_32_ef_min),
        )

        est_32_ef_max = maximum(bounds_32.err_exact_crt)
        est_32_ef_min = minimum(bounds_32.err_exact_crt)
        push!(rows,
            (func="OS2-S-ef-est", m=m, n=n, k=k, phi=phi, num_moduli=num_moduli, err_max=est_32_ef_max, err_min=est_32_ef_min),
        )
    end
end

df = DataFrame(rows)
CSV.write(joinpath(@__DIR__, "results.csv"), df)

