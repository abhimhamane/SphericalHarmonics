abstract type AbstractSHCoeffOrdering end

struct Orderwise <: AbstractSHCoeffOrdering end
struct OrderwiseParityaware <: AbstractSHCoeffOrdering end
struct Degreewise <: AbstractSHCoeffOrdering end

struct SHCoeffParameterVector{
    T,
    V<:AbstractVector{T},
    O<:AbstractSHCoeffOrdering,
    N<:AbstractSHNormalization,
}
    data::V
    lmin::Int
    lmax::Int
    ordering::O
    normalization::N
end

@inline num_coeffs(lmin::Int, lmax::Int) = (lmax+1)^2 - lmin^2

@inline function check_degree_range(coeffs::SHCoefficients, lmin::Int, lmax::Int)
    0<=lmin<=lmax<=coeffs.lmax || throw(ArgumentError("must satisfy 0 ≤ lmin ≤ lmax ≤ $(coeffs.lmax)"))
    return nothing
end

function pack(coeffs::SHCoefficients, ordering::AbstractSHCoeffOrdering; lmin::Int=0, lmax::Int=coeffs.lmax)
    check_degree_range(coeffs, lmin, lmax)
    data = Vector{eltype(coeffs)}(undef, num_coeffs(lmin, lmax))
    pack!(data, coeffs, ordering; lmin, lmax)
    return SHCoeffParameterVector(data, lmin, lmax, ordering, coeffs.normalization)
end

function pack!(data::AbstractVector, coeffs::SHCoefficients, ::Orderwise; lmin::Int=0, lmax::Int=coeffs.lmax)

    idx = firstindex(data)
    @inbounds for order in 0:lmax
        for degree in max(lmin, order):lmax
            data[idx] = C(coeffs, degree, order)
            idx += 1
        end
    end
   
    @inbounds for order in 1:lmax
        for degree in max(lmin, order):lmax
            data[idx] = S(coeffs, degree, order)
            idx += 1
        end
    end
    return data
end

function pack!(data::AbstractVector, coeffs::SHCoefficients, ::OrderwiseParityaware; lmin::Int=0, lmax::Int=coeffs.lmax)
    idx = firstindex(data)
    @inbounds for order in 0:lmax
        for degree in max(lmin, order):2:lmax
            data[idx] = C(coeffs, degree, order)
            idx += 1
        end
    end
    @inbounds for order in 0:lmax
        for degree in max(lmin, order)+1:2:lmax
            data[idx] = C(coeffs, degree, order)
            idx += 1
        end
    end
   
    @inbounds for order in 1:lmax
        for degree in max(lmin, order):2:lmax
            data[idx] = S(coeffs, degree, order)
            idx += 1
        end
    end

    @inbounds for order in 1:lmax
        for degree in max(lmin, order)+1:2:lmax
            data[idx] = S(coeffs, degree, order)
            idx += 1
        end
    end
    return data
end

function unpack(coeffs::SHCoefficients, params::SHCoeffParameterVector)
    params.normalization == coeffs.normalization || throw(ArgumentError("normaliziation mismatch b/w SHCoefficients and parameter vector"))
    params.lmax ≤ coeffs.lmax || throw(DimensionMismatch("coefficients maximum degree too small"))
    unpack!(coeffs, params)
end

function unpack!(coeffs::SHCoefficients, params::SHCoeffParameterVector{T, V, Orderwise},) where {T, V}
    data = params.data
    lmin = params.lmin
    lmax = params.lmax

    idx = firstindex(data)
    @inbounds for order in 0:lmax
        for degree in max(lmin, order):lmax
            setC!(coeffs, degree, order, data[idx])
            idx += 1
        end
    end
    @inbounds for order in 1:lmax
        for degree in max(lmin, order):lmax
            setS!(coeffs, degree, order, data[idx])
            idx += 1
        end
    end
end

function unpack!(coeffs::SHCoefficients, params::SHCoeffParameterVector{T, V, OrderwiseParityaware},) where {T, V}
    data = params.data
    lmin = params.lmin
    lmax = params.lmax

    throw(ArgumentError("Not Implemented"))
end


