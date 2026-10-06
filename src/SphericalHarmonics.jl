module SphericalHarmonics

include("coefficients.jl")
include("orderings.jl")
include("legendre.jl")

export
    # normalization
    AbstractSHNormalization,
    FullyNormalized,
    UnNormalized,

    # coefficient storage
    AbstractSHCoeffRepresentation,
    SHCoefficients,
    SC,
    CS,
    allocate_coefficients,
    C,
    S,
    setC!,
    setS!,

    # parameter ordering
    AbstractSHCoeffOrdering,
    Orderwise,
    OrderwiseParityAware,
    Degreewise,
    SHCoeffParameterVector,
    pack,
    pack!,
    unpack,
    unpack!,

    # Legendre
    AssociatedLegendrePlan,
    anm, bnm, coeffm,
    allocate_plm,
    plm!
end