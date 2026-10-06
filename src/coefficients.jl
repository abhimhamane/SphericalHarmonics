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
    normalization::N
end

Base.eltype(::SHCoefficients{T}) where {T} = T
maxdegree(data::SHCoefficients) = data.lmax
repr(data::SHCoefficients) = data.repr
normalization(data::SHCoefficients) = data.normalization

@inline function check_lm(data::SHCoefficients, degree::Int, order::Int)
    lmax = maxdegree(data)
    0 ≤ order ≤ degree ≤ lmax || throw(BoundsError(data, (degree, order)))
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
@inline function C(data::SHCoefficients{T, A, CS}, degree::Int, order::Int) where {T, A}
    @boundscheck check_lm(data, degree, order)

    @inbounds return data.data[degree+1, order+1]
end

@inline function S(data::SHCoefficients{T, A, CS}, degree::Int, order::Int) where {T, A}
    @boundscheck check_lm(data, degree, order)
    order == 0 && return throw(ArgumentError("order == 0 not defined"))

    @inbounds return data.data[order, degree+1]
end

# setters
@inline function setC!(data::SHCoefficients{T, A, CS}, degree::Int, order::Int, val) where {T, A}
    @boundscheck check_lm(data, degree, order)

    @inbounds data.data[degree+1, order+1] = val
    return data
end

@inline function setS!(data::SHCoefficients{T, A, CS}, degree::Int, order::Int, val) where {T, A}
    @boundscheck check_lm(data, degree, order)
    if order == 0
        throw(ArgumentError("order == 0 not defined"))
        return data
    end
    @inbounds data.data[order, degree+1] = val
    return data
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
    data::SHCoefficients{T, A, SC},
    degree::Int,
    order::Int,
) where {T, A}
    @boundscheck check_lm(data, degree, order)
    @inbounds return data.data[degree+1, data.lmax+1+order]

end

@inline function S(
    data::SHCoefficients{T, A, SC},
    degree::Int,
    order::Int,
) where {T, A}
    @boundscheck check_lm(data, degree, order)
    order == 0 && throw(ArgumentError("order == 0 is invalid"))
    @inbounds return data.data[degree+1, data.lmax+1-order]
end

# setters
@inline function setC!(
    data::SHCoefficients{T, A, SC}, 
    degree::Int, 
    order::Int, 
    val) where {T, A}
    @boundscheck check_lm(data, degree, order)

    @inbounds data.data[degree+1, data.lmax+1+order] = val
    return data
end


@inline function setS!(
    data::SHCoefficients{T, A, SC}, 
    degree::Int, 
    order::Int, 
    val) where {T, A}
    @boundscheck check_lm(data, degree, order)
    order == 0 && throw(ArgumentError("order == 0 is invalid"))
    @inbounds data.data[degree+1, data.lmax+1-order] = val
    return data
end