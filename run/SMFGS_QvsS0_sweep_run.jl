#using Plots
using NPZ
using ProgressBars

include("../SMFGS/SMFGS.jl")
using .SMFGS

# figQtoC
filename="../data/CMF_QvsS0_sweep_figQtoC.npz"
T = GeomRange(1e-3, 1e+4, 7*25)
ζ = [ 0.004 ; 0.006 ; 0.008 ; 0.01 ; 0.02 ; 0.04 ; 0.06 ; 0.08 ; 0.1 ]
S0 = [ 1//2 ; 1 ; 5//2 ; 100 ; 125 ; 250 ; 500 ; 1000 ]

ζ0 = 1
ω0, Γ = 7.0, 5.0
α = (2*ω0^2)*(ζ./ζ0)
θ0 = π/4

println("Running for")
println(min(T...), " ", max(T...), " ", length(T))
println(S0)

s_cmf = zeros((length(ζ), length(T), 3))
s_cwk = zeros((length(ζ), length(T), 3))
s_cgs = zeros((length(T), 3))
s_cus = zeros((length(T), 3))
s_css = zeros((length(ζ), 3))
s_gnd = zeros((length(ζ), 3))
s_qwk = zeros((length(S0), length(ζ), length(T), 3))
s_qgs = zeros((length(S0), length(T), 3))
s_qgv = zeros((length(ζ), length(T), 3))

idx = [id for id in Iterators.product(1:length(ζ), 1:length(T))]
Threads.@threads for d in ProgressBar(idx)
  n,m = d
  s_cmf[n,m,:] = S_CMF(1/T[m], ζ[n], θ0)
  s_cwk[n,m,1] = Sx_CWK(1/T[m], ζ[n], θ0)
  s_cwk[n,m,3] = Sz_CWK(1/T[m], ζ[n], θ0)
  J = LorentzianSD(α[n], ω0, Γ)
  for j in 1:length(S0)
    s_qwk[j,n,m,1] = Sx_QWK(1/T[m], S0[j], J, ζ0, θ0)
    s_qwk[j,n,m,3] = Sz_QWK(1/T[m], S0[j], J, ζ0, θ0)
  end
  s_qgv[n,m,1] = Sx_QGV(1/T[m], ζ[n], θ0)
  s_qgv[n,m,3] = Sz_QGV(1/T[m], ζ[n], θ0)
end

for m in 1:length(T)
  s_cgs[m,3] = Sz_CG(1/T[m])
  s_cus[m,1] = Sx_CUS(1/T[m], θ0)
  s_cus[m,3] = Sz_CUS(1/T[m], θ0)
  for j in 1:length(S0)
    s_qgs[j,m,3] = Sz_QG(1/T[m], S0[j])
  end
end

for n in 1:length(ζ)
  s_css[n,:] = S_CMF_T0(ζ[n], θ0)
  s_gnd[n,:] = S_CMF_ground(ζ[n], θ0)
end
 
npzwrite(filename,
  Dict("S0" => Float64.(S0), "T" => T, "zeta" => ζ, "zeta0" => ζ0,
       "omega0" => ω0, "Gamma" => Γ, "alpha" => α, "theta" => θ0,
       "sqwk" => s_qwk, "sqgs" => s_qgs,
       "scmf" => s_cmf, "scwk" => s_cwk, "scgs" => s_cgs, "scus" => s_cus,
       "scss" => s_css, "sgnd" => s_gnd));

