using Documenter
using SphericalHarmonics

makedocs(
    modules = [SphericalHarmonics],
    sitename = "SphericalHarmonics.jl",
    pages = [
        "Home" => "index.md",
        #"Coefficients" => "coefficients.md",
        #"Parameter orderings" => "orderings.md",
        #"Legendre functions" => "legendre.md",
    ],
)

deploydocs(
    repo = "github.com/abhimhamane/SphericalHarmonics.git",
    devbranch = "main",
)