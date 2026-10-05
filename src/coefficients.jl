abstract type AbstractSHCoeffRepresentation end

abstract type AbstractSHNormalization end

struct FullyNormalized <: AbstractSHNormalization end
struct UnNormalized <: AbstractSHNormalization end

struct SHCoefficients{
        T,
        A<:AbstractArray{T}, 
        L<:AbstractSHCoeffRepresentation, 
        N<:AbstractSHNormalization
    }
    data::A
    lmax::Int
    repr::L
    normaliziation::N
end

Base.eltype(::SHCoefficients{T}) where {T} = T
maxdegree(coeffs::SHCoefficients) = coeffs.lmax
repr(coeffs::SHCoefficients) = coeffs.repr
normaliziation(coeffs::SHCoefficients) = coeffs.normaliziation

@inline function check_lm(coeffs::SHCoefficients, degree::Int, order::Int)
    lmax = maxdegree(coeffs)
    0 ≤ order ≤ degree ≤ lmax || throw(BoundsError(coeffs, (degree, order)))
    return nothing
end


struct CS <: AbstractSHCoeffRepresentation end

# allocation
function allocate_coefficients(::Type{T}, lmax::Int, repr::CS, 
            norm::AbstractSHNormalization=FullyNormalized()) where {T}
    lmax >= 0 || throw(ArgumentError("max. degree (lmax) must be > 0"))
    allocated_coeff = Matrix{T}(undef, lmax+1, lmax+1)
    return SHCoefficients(allocated_coeff, lmax, repr, norm)
end

#getters
@inline function C(coeffs::SHCoefficients{T, A, CS}, degree::Int, order::Int) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)

    @inbounds return coeffs.coeffs[degree+1, order+1]
end

@inline function S(coeffs::SHCoefficients{T, A, CS}, degree::Int, order::Int) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)
    order == 0 && return throw(ArgumentError("order == 0 not defined"))

    @inbounds return coeffs.coeffs[order, degree+1]
end

# setters
@inline function setC!(coeffs::SHCoefficients{T, A, CS}, degree::Int, order::Int, val) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)

    @inbounds coeffs.coeffs[degree+1, order+1] = val
    return coeffs
end

@inline function setS!(coeffs::SHCoefficients{T, A, CS}, degree::Int, order::Int, val) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)
    if order == 0
        throw(ArgumentError("order == 0 not defined"))
        return coeffs
    end
    @inbounds coeffs.coeffs[order, degree+1] = val
    return coeffs
end



struct SC <: AbstractSHCoeffRepresentation end

# allocation
function allocate_coefficients(::Type{T}, lmax::Int, repr::SC,
                norm::AbstractSHNormalization=FullyNormalized()) where {T}
    lmax >= 0 || throw(ArgumentError("max. degree (lmax)  must be > 0"))
    return SHCoefficients(fill(T(NaN), lmax+1, 2*lmax + 1), lmax, repr, norm)
end 

# getters
@inline function C(
    coeffs::SHCoefficients{T, A, SC},
    degree::Int,
    order::Int,
) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)
    @inbounds return coeffs.coeffs[degree+1, coeffs.lmax+1+order]

end

@inline function S(
    coeffs::SHCoefficients{T, A, SC},
    degree::Int,
    order::Int,
) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)
    order == 0 && throw(ArgumentError("order == 0 is invalid"))
    @inbounds return coeffs.coeffs[degree+1, coeffs.lmax+1-order]
end

# setters
@inline function setC!(
    coeffs::SHCoefficients{T, A, SC}, 
    degree::Int, 
    order::Int, 
    val) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)

    @inbounds coeffs.coeffs[degree+1, coeffs.lmax+1+order] = val
    return coeffs
end


@inline function setS!(
    coeffs::SHCoefficients{T, A, SC}, 
    degree::Int, 
    order::Int, 
    val) where {T, A}
    @boundscheck check_lm(coeffs, degree, order)
    order == 0 && throw(ArgumentError("order == 0 is invalid"))
    @inbounds coeffs.coeffs[degree+1, coeffs.lmax+1-order] = val
    return coeffs
end

struct Orderwise <: AbstractSHCoeffRepresentation end

coeff = allocate_coefficients(Float64, 5, SC(), FullyNormalized())

coeff.data
