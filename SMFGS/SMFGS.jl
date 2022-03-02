module SMFGS

using LinearAlgebra
using Kronecker
using Statistics
using StaticArrays
using Quadmath, DoubleFloats
using QuadGK
using ForwardDiff
using FiniteDifferences
using Optim, Roots

export Sz_CG,
       Sz_QG,
       Sz_CWK, Sx_CWK,
       Sz_QWK, Sx_QWK,
       Sz_CUS, Sx_CUS,
       Sz_QUS, Sx_QUS,
       S_CMF,  S_CMF_T0, S_CMF_ground, Pdens_CMF,
       S_QRC, precompute_operators_QRC,
       Sz_QGV, Sx_QGV, QGV_bound,
       GeomRange,
       LorentzianSD, OhmicSD, DrudeLorentzSD

include("auxmath.jl")

include("adim.jl")

include("SpectralDensity.jl")
include("SDIntegrals.jl")

include("CG.jl")
include("QG.jl")

include("CWK.jl")
include("QWK.jl")

include("CUS.jl")
include("QUS.jl")

include("CMF.jl")

include("QRC.jl")
include("QGV.jl")

end
