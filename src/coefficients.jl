abstract type AbstractSHCoeffRepresentation end

abstract type AbstractSHNormaliziation end

struct FullyNormalized <: AbstractSHNormaliziation end
struct UnNormalized <: AbstractSHNormaliziation end

struct SHCoefficients{
        T,
        A<:AbstractArray{T}, 
        L<:AbstractSHCoeffRepresentation, 
        N<:AbstractSHNormaliziation
    }
    coeffs::A,
    lmax::Int,
    repr:L,
    normaliziation::N
end

Base.eltype(::SHCoefficients{T}) where {T} = T
maxdegree(field::SHCoefficients) = field.lmax
repr(field::SHCoefficients) = field.repr
normaliziation(field::SHCoefficients) = field.normaliziation

@inline function check_lm(field::SHCoefficients, degree::Int, order::Int)
    lmax = maxdegree(field)
    0 ≤ order ≤ degree ≤ lmax || throw(BoundsError(field, (degree, order)))
    return nothing
end


struct CS <: AbstractSHCoeffRepresentation end
struct SC <: AbstractSHCoeffRepresentation end
struct Orderwise <: AbstractSHCoeffRepresentation end

@inline function C(field::SHCoefficients{T, A, CS}, degree::Int, order::Int) where {T, A}
    

end