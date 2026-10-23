using LinearAlgebra

const moduli_list = Int[
    256, 255, 253, 251, 247, 241, 239, 233, 229, 227,
    223, 217, 211, 199, 197, 193, 191, 181, 179, 173,
    167, 163, 157, 151, 149, 139, 137, 131, 127, 113,
    109, 107, 103, 101, 97, 89, 83, 79, 73, 71
]

const CRT_ERR_FREE_PRECISION = 512

struct Table
    moduli::Vector{Int}
    P_hi::Float64
    P_lo::Float64
    P_inv::Float64
    log2P::Float32
    s_hi::Vector{Float64}
    s_lo::Vector{Float64}
    s::Vector{Float64}
    s_err_free::Vector{BigFloat}
    P_err_free::BigFloat
end

@inline function ceil_log2(x::Integer)
    x <= 1 && return 0
    return ndigits(x - 1, base=2)
end

function gen_const(NUM_MODULI::Int)
    1 <= NUM_MODULI <= length(moduli_list) ||
        throw(ArgumentError("NUM_MODULI must be in 1:$(length(moduli_list))"))

    moduli = moduli_list[1:NUM_MODULI]
    P_int = prod(BigInt.(moduli))
    P_hi = Float64(P_int)
    P_lo = Float64(P_int - BigInt(P_hi))

    P_inv, log2P = setprecision(BigFloat, 512) do
        setrounding(BigFloat, RoundNearest) do
            pinv = Float64(inv(BigFloat(P_int)), RoundNearest)
            lbound = setrounding(BigFloat, RoundDown) do
                Float32((log2(BigFloat(P_int - 1)) - 1) / 2, RoundDown)
            end
            (pinv, lbound)
        end
    end

    Pi = P_int .÷ moduli
    q = invmod.(Pi, moduli)
    s_mp = q .* Pi
    sum_half_moduli = sum(fld.(moduli, 2))
    bits = 53 - ceil_log2(sum_half_moduli)
    ufp = maximum(ndigits(abs(si), base=2) - 1 for si in s_mp)
    sft = bits - ufp - 1
    if sft < 0
        shift = -sft
        s_hi_int = (s_mp .>> shift) .<< shift
    else
        s_hi_int = copy(s_mp)
    end
    s_lo_int = s_mp .- s_hi_int
    s_hi = Float64.(s_hi_int)
    s_lo = Float64.(s_lo_int)
    s = Float64.(s_mp)

    s_err_free = [BigFloat(si, RoundNearest; precision=CRT_ERR_FREE_PRECISION)
                  for si in s_mp]
    P_err_free = BigFloat(P_int, RoundNearest; precision=CRT_ERR_FREE_PRECISION)

    return Table(moduli, P_hi, P_lo, P_inv, log2P, s_hi, s_lo, s,
        s_err_free, P_err_free)
end

function calc_sft(amax::AbstractVector{<:Integer}, log2P::Float32)
    all(x -> x >= 0, amax) || throw(ArgumentError("amax must be nonnegative"))
    amax_f = Float32.(amax, RoundUp)
    log2amax = log2.(amax_f)
    mh4u_ru = Float32(0.5 / (1.0 - 2.0^-22), RoundUp)
    sft_d = Float64(log2P) .- Float64(mh4u_ru) .* Float64.(log2amax)
    return [iszero(amax[i]) ? Int32(0) :
            floor(Int32, Float32(sft_d[i], RoundDown)) for i in eachindex(amax)]
end

ufp_exponent(x::AbstractFloat) = iszero(x) ? 0 : exponent(x)

function symmod(A::AbstractMatrix{<:Real}, p::Integer)
    1 <= p <= 256 ||
        throw(ArgumentError("p must be in 1:256"))

    lo, hi = -fld(p, 2), cld(p, 2)

    return map(A) do x
        r = rem(x, p)
        if r < lo
            r += p
        elseif r >= hi
            r -= p
        end
        Int8(r)
    end
end

function _int8_matmul(A::AbstractMatrix{Int8}, B::AbstractMatrix{Int8})
    k = size(A, 2)
    k == size(B, 1) || throw(DimensionMismatch("incompatible matrix dimensions"))

    C = Float64.(A) * Float64.(B)
    return Int32.(C)
end

function scaling(A::AbstractMatrix{TA}, B::AbstractMatrix{TB}, c::Table) where
{TA<:Union{Float32,Float64},TB<:Union{Float32,Float64}}
    m, k = size(A)
    k == size(B, 1) || throw(DimensionMismatch("size(A, 2) must equal size(B, 1)"))
    n = size(B, 2)
    all(isfinite, A) && all(isfinite, B) ||
        throw(ArgumentError("A and B must contain only finite values"))
    if m == 0 || n == 0 || k == 0
        return zeros(Float64, m, k), zeros(Float64, k, n),
        zeros(Int16, m), zeros(Int16, n)
    end
    log2P = c.log2P

    ufpA = vec(ufp_exponent.(maximum(abs, A; dims=2)))
    mu1 = Int16.(5 .- ufpA)
    Au = Int8.(ceil.(ldexp.(abs.(Float64.(A)), mu1)))

    ufpB = vec(ufp_exponent.(maximum(abs, B; dims=1)))
    nu1 = Int16.(5 .- ufpB)
    Bu = Int8.(ceil.(ldexp.(abs.(Float64.(B)), transpose(nu1))))

    Cu = _int8_matmul(Au, Bu)
    mu = Int16.(mu1 .+ calc_sft(vec(maximum(Cu; dims=2)), log2P))
    nu = Int16.(nu1 .+ calc_sft(vec(maximum(Cu; dims=1)), log2P))

    Aint = trunc.(ldexp.(Float64.(A), mu))
    Bint = trunc.(ldexp.(Float64.(B), transpose(nu)))
    return Aint, Bint, mu, nu
end

function crt(Aint::AbstractMatrix{<:AbstractFloat},
    Bint::AbstractMatrix{<:AbstractFloat}, c::Table, ::Type{Float64})
    size(Aint, 2) == size(Bint, 1) ||
        throw(DimensionMismatch("incompatible matrix dimensions"))
    m, n = size(Aint, 1), size(Bint, 2)
    C1 = zeros(Float64, m, n)
    C2 = zeros(Float64, m, n)
    for i in eachindex(c.moduli)
        p = c.moduli[i]
        Ai = symmod(Aint, p)
        Bi = symmod(Bint, p)
        Ci = _int8_matmul(Ai, Bi)
        Wi = symmod(Ci, p)
        C1 = fma.(c.s_hi[i], Float64.(Wi), C1)
        C2 = fma.(c.s_lo[i], Float64.(Wi), C2)
    end
    Q = round.(c.P_inv .* C1)
    return fma.(-Q, c.P_lo, fma.(-Q, c.P_hi, C1) .+ C2)
end

function crt(
    Aint::AbstractMatrix{<:AbstractFloat},
    Bint::AbstractMatrix{<:AbstractFloat},
    c::Table, ::Type{Float32}
)
    size(Aint, 2) == size(Bint, 1) ||
        throw(DimensionMismatch("incompatible matrix dimensions"))
    m, n = size(Aint, 1), size(Bint, 2)
    C1 = zeros(Float64, m, n)
    for i in eachindex(c.moduli)
        p = c.moduli[i]
        Ai = symmod(Aint, p)
        Bi = symmod(Bint, p)
        Ci = _int8_matmul(Ai, Bi)
        Wi = symmod(Ci, p)
        C1 = fma.(c.s[i], Float64.(Wi), C1)
    end
    Q = round.(c.P_inv .* C1)
    return fma.(-Q, c.P_hi, C1)
end

function crt(Aint::AbstractMatrix{TA}, Bint::AbstractMatrix{TB}, c::Table) where
{TA<:AbstractFloat,TB<:AbstractFloat}
    T = TA === Float32 && TB === Float32 ? Float32 : Float64
    return crt(Aint, Bint, c, T)
end

function crt_err_free(Aint::AbstractMatrix{TA}, Bint::AbstractMatrix{TB}, c::Table) where
{TA<:Union{Float32,Float64},TB<:Union{Float32,Float64}}
    k = size(Aint, 2)
    k == size(Bint, 1) ||
        throw(DimensionMismatch("incompatible matrix dimensions"))
    k <= (Int64(1) << 39) ||
        throw(OverflowError("Float64 exactness bound exceeded"))

    return setprecision(BigFloat, CRT_ERR_FREE_PRECISION) do
        setrounding(BigFloat, RoundNearest) do
            m, n = size(Aint, 1), size(Bint, 2)
            C1 = zeros(BigFloat, m, n)
            for i in eachindex(c.moduli)
                p = c.moduli[i]
                Ai = symmod(Aint, p)
                Bi = symmod(Bint, p)
                Ci = _int8_matmul(Ai, Bi)
                Wi = symmod(Ci, p)
                C1 .= fma.(c.s_err_free[i], BigFloat.(Wi), C1)
            end

            P = c.P_err_free
            halfP = P / 2
            return map(C1) do x
                r = mod(x, P)
                if r >= halfP
                    r -= P
                end
                Float64(r, RoundNearest)
            end
        end
    end
end

function undo_scaling(Cint::AbstractMatrix{Float64},
    mu::AbstractVector{<:Integer}, nu::AbstractVector{<:Integer})
    size(Cint) == (length(mu), length(nu)) ||
        throw(DimensionMismatch("scale lengths must match the result dimensions"))
    return ldexp.(Cint, (-Int.(mu)) .- transpose(Int.(nu)))
end

function os2(A::AbstractMatrix{TA}, B::AbstractMatrix{TB}, NUM_MODULI::Int) where
{TA<:Union{Float32,Float64},TB<:Union{Float32,Float64}}
    c = gen_const(NUM_MODULI)
    Aint, Bint, mu, nu = scaling(A, B, c)
    T = TA === Float32 && TB === Float32 ? Float32 : Float64
    Cint = crt(Aint, Bint, c, T)
    C = undo_scaling(Cint, mu, nu)
    if T === Float32
        return Float32.(C)
    end
    return C
end

function oz2_err_free_crt(A::AbstractMatrix{TA}, B::AbstractMatrix{TB},
    NUM_MODULI::Int) where
{TA<:Union{Float32,Float64},TB<:Union{Float32,Float64}}
    c = gen_const(NUM_MODULI)
    Aint, Bint, mu, nu = scaling(A, B, c)
    Cint = crt_err_free(Aint, Bint, c)
    C = undo_scaling(Cint, mu, nu)
    if TA === Float32 && TB === Float32
        return Float32.(C)
    end
    return C
end
