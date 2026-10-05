abstract type AbstractSHCoeffOrdering end

struct Orderwise <: AbstractSHCoeffOrdering end
struct OrderwiseParityaware <: AbstractSHCoeffOrdering end
struct Degreewise <: AbstractSHCoeffOrdering end

struct SHCoeffParamVector{
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

function pack(coeffs::SHCoefficients, ::Orderwise, lmin::Int=0, lmax::Int=coeffs.lmax)
    data = Vector{eltype(coeffs)}(undef, num_coeffs(lmin, lmax))

    
end