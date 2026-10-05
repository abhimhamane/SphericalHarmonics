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

function pack(coeffs::SHCoefficients, ordering::AbstractSHCoeffOrdering, lmin::Int=0, lmax::Int=coeffs.lmax)
    data = Vector{eltype(coeffs)}(undef, num_coeffs(lmin, lmax))
    pack!(data, coeffs, ordering; lmin, lmax)
    return SHCoeffParameterVector(data, lmin, lmax, ordering, coeffs.normalization)
end

function pack!(data::AbstractVector, coeffs::SHCoefficients, ::Orderwise, lmin::Int=0, lmax::Int=coeffs.lmax)
    idx = firstindex(data)
    @inbounds for order in 0:lmax
        for degree in max(lmin, order):lmax
            data[idx] = C(coeffs, degree, order)
            idx += 1
        end
    end
   
    @inbounds for order in 1:lmax
        for degree in max(lmin, order):lmax
            data[idx] = S[coeffs, degree, order]
            idx += 1
        end
    end
    return data
end

function pack!(data::AbstractVector, coeffs::SHCoefficients, ::OrderwiseParityaware, lmin::Int=0, lmax::Int=coeffs.lmax)
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
            data[idx] = S[coeffs, degree, order]
            idx += 1
        end
    end

    @inbounds for order in 1:lmax
        for degree in max(lmin, order)+1:2:lmax
            data[idx] = S[coeffs, degree, order]
            idx += 1
        end
    end
    return data
end
