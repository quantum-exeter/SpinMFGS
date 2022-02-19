using NPZ
using ProgressBars

include("../SMFGS/SMFGS.jl")
using .SMFGS

filename = "../data/CMF_Pdens_sweep.npz"
T = [1.0]
ζ = [0.0 ; 1.0]

θ0 = π/4

θ = LinRange(0.0,  π, 200)
φ = LinRange(0.0, 2π, 200)

Pdens = zeros((length(ζ), length(T), length(θ), length(φ)))

idx = [id for id in Iterators.product(1:length(θ), 1:length(φ), 1:length(ζ), 1:length(T))]
Threads.@threads for id in ProgressBar(idx)
  n,m,j,k = id
  Pdens[j,k,n,m] = Pdens_CMF(θ[n],φ[m],1/T[k],ζ[j],θ0)
end

npzwrite(filename,
  Dict("T" => T, "zeta" => ζ, "theta" => θ0, "sphtheta" => θ, "sphphi" => φ, "Pdens" => Pdens));

