"""    anm(n::Int, m::Int)
computes the coefficient for normalizing associated legndre polynomials.

# Reference:

    1 - ``W_{lm} = \\sqrt{\\frac{(2l-1)(2l+1)}{(l+m)(l-m)}}`` from Prof. Nico Sneeuw, Lecture Notes on Physical Geodesy, 
        Chapter 6 - Solving Laplace's equation, Pg. 88, from - 15th June 2006, 
        URL (as of 16 Nov. 2025) - https://www.gis.uni-stuttgart.de/lehre/campus-docs/LNErdm.pdf
    
    2 - Holmes and Featherstone, “A Unified Approach to the Clenshaw Summation and the Recursive Computation 
        of Very High Degree and Order Normalised Associated Legendre Functions.” Journal of Geodesy. 
        DOI-10.1007/s00190-002-0216-2

"""
@inline function anm(::FullyNormalized, n::Int, m::Int)
        return sqrt(((2*n - 1)*(2*n + 1))/((n-m)*(n+m)))
end 

"""
    bnm(n::Int, m::Int)
computes the coefficient for normalizing associated legndre polynomials.

# Reference: 
    
    1 - simplified ``W_{lm}*W_{l-1,m}^{-1}`` from Prof. Nico Sneeuw, Lecture Notes on Physical Geodesy, 
        Chapter 6 - Solving Laplace's equation, Pg. 88, from - 15th June 2006, 
        URL (as of 16 Nov. 2025) https://www.gis.uni-stuttgart.de/lehre/campus-docs/LNErdm.pdf
    
    2 - Holmes and Featherstone, “A Unified Approach to the Clenshaw Summation and the Recursive Computation 
        of Very High Degree and Order Normalised Associated Legendre Functions.” Journal of Geodesy. 
        DOI-10.1007/s00190-002-0216-2
"""
@inline function bnm(::FullyNormalized, n::Int, m::Int)
    num = (2*n + 1)*(n+m-1)*(n-m-1)
    denom = (n-m)*(n+m)*(2*n-3)
    return sqrt(num/denom)
end

"""
    coeffm(m::Int)
computes the coefficient for normalizing sectorial recursion part of PLM computation.

# Reference: 

    1 - ``W_{mm} = \\sqrt{\\frac{2m+1}{2m}}`` from Prof. Nico Sneeuw, Lecture Notes on Physical Geodesy,
        Chapter 6 - Solving Laplace's equation, Pg. 88, from - 15th June 2006, 
        URL (as of 16 Nov. 2025) https://www.gis.uni-stuttgart.de/lehre/campus-docs/LNErdm.pdf

    2 - Holmes and Featherstone, “A Unified Approach to the Clenshaw Summation and the Recursive Computation 
        of Very High Degree and Order Normalised Associated Legendre Functions.” Journal of Geodesy. 
        DOI-10.1007/s00190-002-0216-2
"""
@inline function coeffm(::FullyNormalized, m::Int)
    return sqrt((2*m + 1)/(2*m))
end

"""
    Plm(l::Int, m::Int, theta::Float64)
Computes associated legendre polynomials.

Approach follows column recursion formulas in a 3-step process (more details see reference).
    1. for a given degree (`m`), compute `P_{mm}` using sectorial recursion
    2. compute the next off-diagonal element `P_{l1}^m` such that `l1 = m + 1` using special case of generic recursion
    3. use the 2 point recursive formula `P_{l-1}^m, P_{l-2}^m` along the column.
    4. also handles special cases and invlid inputs

    # Inputs
    - ::AbstractSHNormalization (currently only FullyNormalized() implemented)
    - `l`: degree 
    - `m`: order
    - `theta`: co-latitude [radians]

# Notes
    - no input check performed on theta
    - not defined for `l < m`

# Reference
    1 - Prof. Nico Sneeuw, Lecture Notes on Physical Geodesy, 
        Chapter 6 - Solving Laplace's equation, Pg. 88, from - 15th June 2006, 
        URL (as of 16 Nov. 2025) https://www.gis.uni-stuttgart.de/lehre/campus-docs/LNErdm.pdf

    2 - Holmes and Featherstone, “A Unified Approach to the Clenshaw Summation and the Recursive Computation 
        of Very High Degree and Order Normalised Associated Legendre Functions.” Journal of Geodesy. 
        DOI-10.1007/s00190-002-0216-2
"""
function Plm(::FullyNormalized, l::Int, m::Int, theta::Float64)
    if l < abs(m)
        error("Invalid Input: order cannot be greater than degree.")
    end
    
    
    #TODO - dealing with negative order

    P = 1.0
    Pmm = 1.0
    l == 0 && m == 0 && return P

    u = sin(theta)
    t = cos(theta)
    Pmm = sqrt(3)*u
    l == 1 && m == 1 && return Pmm

    # sectoral Plm - diagonal
    @inbounds for mm = 2:m
        Pmm = u * coeffm(FullyNormalized(), mm) * Pmm
    end
    l == m && return Pmm

    Pmin2 = Pmm
    if m == 0
        Pmin2 = 1
    end
    Pmin1 =  anm(FullyNormalized(), m+1, m)* t * Pmin2
    l == m+1 && return Pmin1

    @inbounds for l1 in m+2:l
        P = anm(FullyNormalized(), l1, m)*t*Pmin1 - bnm(FullyNormalized(), l1, m)*Pmin2
        Pmin2, Pmin1 = Pmin1, P
    end
    return P
    
end

#Plm(FullyNormalized(), 2,0, π/3)

struct AssociatedLegendrePlan{T, N}
    lmax::Int
    normalization::N
    anm::Matrix{T}
    bnm::Matrix{T}
    coeffm::Vector{T}
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

function allocate_plm(plan::AssociatedLegendrePlan{T, N}) where {T<:AbstractFloat, N<:AbstractSHNormalization}
    return zeros(T, plan.lmax+1, plan.lmax+1)
end

function plm!(plm_mat::AbstractMatrix{T}, plan::AssociatedLegendrePlan{T}, θ::Real) where {T}
    u, t = sincos(θ)
    plm_mat[1, 1] = 1.0
    plm_mat[2, 1] = sqrt(3)*t
    plm_mat[2, 2] = sqrt(3)*u

    # sectoral Plm - diagonal
    @inbounds for m = 2:plan.lmax
        plm_mat[m+1, m+1] = u*plan.coeffm[m+1]*plm_mat[m,m]
        plm_mat[m+1, m] = plan.anm[m+1,m]*t*plm_mat[m,m]
    end

    @inbounds for m = 0:plan.lmax-2
        @inbounds for n = m+2:plan.lmax
            plm_mat[n+1, m+1] = plan.anm[n+1,m+1]*t*plm_mat[n, m+1] - plan.bnm[n+1,m+1]*plm_mat[n-1, m+1] 
        end
    end
end


@inline function Fnm(::FullyNormalized, n::Int, m::Int)
    n == m && return 0
    return sqrt(((n-m) * (n+m) * (2*n + 1)) / (2*n - 1))
end

