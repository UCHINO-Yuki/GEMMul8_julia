import Pkg

let required_packages = ["CSV", "DataFrames"]
    missing_packages = filter(pkg -> isnothing(Base.find_package(pkg)), required_packages)
    if !isempty(missing_packages)
        Pkg.add(missing_packages)
    end
end

using LinearAlgebra
using Random
using CSV, DataFrames

const BITS = 512
const MODULI = Int[
    256, 255, 253, 251, 247, 241, 239, 233, 229, 227,
    223, 217, 211, 199, 197, 193, 191, 181, 179, 173,
    167, 163, 157, 151, 149, 139, 137, 131, 127, 113,
    109, 107, 103, 101, 97, 89, 83, 79, 73, 71,
    67, 61, 59, 53, 47, 43, 41, 37, 29]

down(f) = setrounding(f, BigFloat, RoundDown)
up(f) = setrounding(f, BigFloat, RoundUp)
work(f) = setprecision(() -> setrounding(f, BigFloat, RoundNearest), BigFloat, BITS)

function largest_indices(v, width::Int)
    width in 1:4 || throw(ArgumentError("star_width must be in 1:4"))
    isempty(v) && throw(ArgumentError("cannot select from an empty vector"))
    q = min(width, length(v))
    idx = zeros(Int, q)
    used = 0
    for i in eachindex(v)
        pos = used + 1
        while pos > 1 && v[i] > v[idx[pos-1]]
            pos <= q && (idx[pos] = idx[pos-1])
            pos -= 1
        end
        pos <= q && (idx[pos] = i)
        used = min(q, used + 1)
    end
    Tuple(idx)
end

function partial_product_lower(X, Y; width::Int=3)
    m, k = size(X)
    size(Y, 1) == k || throw(DimensionMismatch("incompatible matrices"))
    n = size(Y, 2)
    width in 1:4 || throw(ArgumentError("star_width must be in 1:4"))
    Q = zeros(BigFloat, m, n)
    down() do
        for h in 1:k
            ir = largest_indices(view(X, :, h), width)
            jc = largest_indices(view(Y, h, :), width)
            for i in ir, j in 1:n
                Q[i, j] += X[i, h] * Y[h, j]
            end
            for j in jc, i in 1:m
                i in ir && continue
                Q[i, j] += X[i, h] * Y[h, j]
            end
        end
    end
    Q
end

function minvector_lower(X, Y; star_width::Int=3)
    a = vec(maximum(X; dims=2))
    b = vec(maximum(Y; dims=1))
    xi = down(() -> vec(minimum(X ./ a; dims=1)))
    eta = down(() -> vec(minimum(Y ./ transpose(b); dims=2)))
    rankone_upper = up(() -> a .* transpose(xi))
    residual_lower = down(() -> max.(0, X .- rankone_upper))
    r = down(() -> residual_lower * eta)
    c = down(() -> transpose(Y) * xi)
    base = down(() -> r .* transpose(b) .+ a .* transpose(c))
    residualB_lower = down(() -> max.(0, Y .- up(() -> eta .* transpose(b))))
    correction = partial_product_lower(residual_lower, residualB_lower; width=star_width)
    down(() -> base .+ correction)
end

function prepare_selection(
    A::AbstractMatrix{T},
    B::AbstractMatrix{T};
    Cbar=nothing,
    star_width::Int=3
) where {T<:Union{Float32,Float64}}
    star_width in 1:4 || throw(ArgumentError("star_width must be in 1:4"))
    m, k = size(A)
    size(B, 1) == k || throw(DimensionMismatch("A and B are incompatible"))
    n = size(B, 2)
    m > 0 && n > 0 && 1 <= k <= 2^17 || throw(ArgumentError("require m,n > 0 and 1 <= k <= 2^17"))
    all(isfinite, A) && all(isfinite, B) || throw(ArgumentError("inputs must be finite"))

    amax = vec(maximum(abs, A; dims=2))
    bmax = vec(maximum(abs, B; dims=1))
    all(>(0), amax) && all(>(0), bmax) || throw(ArgumentError("zero rows of A or zero columns of B are excluded"))

    return work() do
        X, Y = abs.(BigFloat.(A)), abs.(BigFloat.(B))
        alpha, beta = exponent.(amax), exponent.(bmax)
        mu0, nu0 = 5 .- alpha, 5 .- beta
        Ab = ceil.(Int64, ldexp.(X, mu0))
        Bb = ceil.(Int64, ldexp.(Y, transpose(nu0)))
        if Cbar !== nothing
            size(Cbar) == (m, n) || throw(DimensionMismatch("wrong Cbar size"))
            eltype(Cbar) <: Integer || throw(ArgumentError("Cbar must be integer"))
        end
        C = Cbar === nothing ? Ab * Bb : Int64.(Cbar)
        all(z -> 0 <= z <= 2^29, C) || throw(ArgumentError("Cbar out of range"))
        Ah = ldexp.(BigFloat.(Ab), -mu0)
        Bh = ldexp.(BigFloat.(Bb), transpose(-nu0))
        Gh = ldexp.(BigFloat.(C), -mu0 .- transpose(nu0))

        sigmaA, sigmaB, Delta = up() do
            sa = vec(sum(X; dims=2))
            sb = vec(sum(Y; dims=1))
            sha = vec(sum(Ah; dims=2))
            shb = vec(sum(Bh; dims=1))
            da = vec(maximum(Ah .- X; dims=2))
            db = vec(maximum(Bh .- Y; dims=1))
            dg_old = min.(da .* transpose(sb) .+ sha .* transpose(db), da .* transpose(shb) .+ sa .* transpose(db))

            countA = vec(sum(Ab .> 0; dims=2))
            countB = vec(sum(Bb .> 0; dims=1))
            vartheta = max.(0, countA .+ transpose(countB) .- k)

            correction = down(() -> vartheta .* da .* transpose(db))
            dg_shift = da .* transpose(shb) .+ sha .* transpose(db) .- correction
            sa, sb, min.(dg_old, dg_shift)
        end
        G_lower = down(() -> max.(0, Gh .- Delta))
        G_lower = max.(G_lower, minvector_lower(X, Y; star_width))
        rowmax = max.(1, vec(maximum(C; dims=2)))
        colmax = max.(1, vec(maximum(C; dims=1)))
        x = up(() -> ldexp.(sqrt.(BigFloat.(rowmax)), alpha))
        y = up(() -> ldexp.(sqrt.(BigFloat.(colmax)), beta))
        (; k, b=(T === Float32 ? 32 : 64), mu0, nu0, rowmax, colmax,
            sigmaA, sigmaB, x, y, G_lower, Cbar=C)
    end
end

function scaling_shifts(data, P)
    plo = down(() -> (log2(BigFloat(P - 1)) - 1) / 2)
    phi = up(() -> (log2(BigFloat(P - 1)) - 1) / 2)
    p32 = Float32(plo, RoundDown)
    p32 == Float32(phi, RoundDown) || error("P' rounding is unresolved; increase BITS")
    coef = Float32(down(() -> -BigFloat(2)^23 / (BigFloat(2)^24 - 4)), RoundDown)
    function shift(z)
        e = log2(Float32(z, RoundUp))
        s = BigFloat(coef) * BigFloat(e) + BigFloat(p32)
        floor(Int, Float32(s, RoundDown))
    end
    data.mu0 .+ shift.(data.rowmax), data.nu0 .+ shift.(data.colmax)
end

function candidate_bounds(data, N)
    P = prod(BigInt.(MODULI[1:N]))
    rho = sum(fld(p, 2) for p in MODULI[1:N])
    mu, nu = scaling_shifts(data, P)
    u = ldexp(BigFloat(1), data.b == 32 ? -24 : -53)
    u64 = ldexp(BigFloat(1), -53)
    chi = data.b == 32 ? 1 : 3
    t2 = up(() -> inv(BigFloat(32 * (P - 1))))
    t = up(() -> sqrt(t2))
    S, Tn, Gamma, zeta, uprime = up() do
        S = t2 .* data.x .* transpose(data.y)
        Tn = t .* data.sigmaA .* transpose(data.y) .+ t .* data.x .* transpose(data.sigmaB) .+ data.k .* S
        Gamma = ldexp.(S, mu .+ transpose(nu))
        zeta = data.b == 32 ?
               (1 + u) * (N + 2) * u64 * rho * BigFloat(P) :
               (1 + 3u) * BigFloat(2)^(1 + ndigits(rho - 1; base=2)) * (N + 2) * u64^2 * rho * BigFloat(P)
        S, Tn, Gamma, zeta, u / (1 + u)
    end
    (; S, Tn, Gamma, zeta, uprime, u, chi, mu, nu)
end

function check_selection(data, N::Integer; kprime=data.k)
    2 <= N <= length(MODULI) || throw(ArgumentError("N must be in 2:49"))
    return work() do
        kp = down() do
            kprime === :sqrt ? sqrt(BigFloat(data.k)) :
            kprime === :log2 ? log2(BigFloat(data.k)) : BigFloat(kprime)
        end
        isfinite(kp) && 0 < kp <= data.k ||
            throw(ArgumentError("require 0 < kprime <= k"))
        c = candidate_bounds(data, N)

        crt_factor = up(() -> c.chi .* c.Gamma)
        budget = down(() -> c.u .* (kp .- crt_factor))
        lhs = up(() -> c.Tn .+ c.zeta .* c.S)
        rhs = down(() -> budget .* data.G_lower)
        finite = all(>(0), budget) && all(lhs .<= rhs)

        budget_exact = down(() -> kp * c.u - c.uprime)
        lhs_exact = up(() -> (1 + c.uprime) .* c.Tn)
        rhs_exact = down(() -> budget_exact .* data.G_lower)
        exact = budget_exact > 0 && all(lhs_exact .<= rhs_exact)
        (; finite, exact)
    end
end

function select_moduli(data; kprime=data.k, Ns=2:49)
    finite_N = exact_N = nothing
    for N in sort(unique(collect(Ns)))
        ok = check_selection(data, N; kprime)
        isnothing(finite_N) && ok.finite && (finite_N = N)
        isnothing(exact_N) && ok.exact && (exact_N = N)
        !isnothing(finite_N) && !isnothing(exact_N) && break
    end
    (finite=finite_N, exact=exact_N)
end

gen_rand(m::Int64, n::Int64) = rand(Float64, m, n)
gen_randn(m::Int64, n::Int64) = randn(Float64, m, n)
gen_mat(m::Int64, n::Int64, phi::AbstractFloat) = (gen_rand(m, n) .- 0.5) .* exp.(gen_randn(m, n) .* phi)

function run_experiments(; output=joinpath(@__DIR__, "select_moduli.csv"), star_width::Int=3)
    Random.seed!(1)
    m, n, k = 8, 8, 8192
    rows = NamedTuple[]
    for phi in (0.0, 1.0, 2.0, 4.0, 8.0, 16.0)
        A = gen_mat(m, k, phi)
        B = gen_mat(k, n, phi)
        for FT in (Float64, Float32)
            prefix = FT === Float64 ? "DGEMM" : "SGEMM"
            excluded = FT === Float32 && phi == 16.0
            data = excluded ? nothing : prepare_selection(FT.(A), FT.(B); star_width)
            for (tag, kp) in (("k", k), ("sqrtk", :sqrt))
                result = excluded ? (finite="N/A", exact="N/A") :
                         select_moduli(data; kprime=kp, Ns=2:49)
                push!(rows, (func=prefix * "-" * tag,
                    m=m, n=n, k=k, phi=phi,
                    finite=result.finite, exact=result.exact))
            end
        end
    end
    df = DataFrame(rows)
    CSV.write(output, df;
        transform=(_, value) -> something(value, missing), missingstring="-")
    return df
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiments()
end
