@inline function anm(::FullyNormalized, n::Int, m::Int)
        return sqrt(((2*n - 1)*(2*n + 1))/((n-m)*(n+m)))
end 

@inline function bnm(::FullyNormalized, n::Int, m::Int)
    num = (2*n + 1)*(n+m-1)*(n-m-1)
    denom = (n-m)*(n+m)*(2*n-3)
    return sqrt(num/denom)
end

@inline function coeffm(::FullyNormalized, m::Int)
    return sqrt((2*m + 1)/(2*m))
end

@inline function Fnm(::FullyNormalized, n::Int, m::Int)
    n == m && return 0
    return sqrt(((n-m) * (n+m) * (2*n + 1)) / (2*n - 1))
end

struct AssociatedLegendrePlan{N, A, B, V}
    lmax::Int
    normalization::N
    anm::A
    bnm::B
    coeffm::V
end

function AssociatedLegendrePlan(::Type{T}, lmax::Int, 
    norm::N=FullyNormalized()) where {
        T<:AbstractFloat,
        N<:AbstractSHNormalization,
    }
    lmax ≥ 0 || throw(ArgumentError("lmax must be > 0"))
    coeff_anm = zeros(T, lmax+1, lmax+1)
    coeff_bnm = zeros(T, lmax+1, lmax+1)
    sectorial_coeffm = zeros(T, lmax+1)

    @inbounds for order in 1:lmax
        sectorial_coeffm[order+1] = coeffm(norm, order)
    end

    @inbounds for order in 0:lmax
        for degree in order+1:lmax
            coeff_anm[degree+1, order+1] = anm(norm, degree, order)
        end

        for degree in order+2:lmax
            coeff_bnm[degree+1, order+1] = bnm(norm, degree, order)
        end
    end
    return AssociatedLegendrePlan(lmax, norm, coeff_anm, coeff_bnm, sectorial_coeffm)
end

function allocate_plm(plan::AssociatedLegendrePlan{T}) where {T}
    return Matrix{T}(undef, plan.lmax+1, plan.lmax+1)
end

function plm!(plm_mat::AbstractMatrix{T}, plan::AssociatedLegendrePlan{T}, θ::Real) where {T}
    error("plm! not implemented yet...")
end
