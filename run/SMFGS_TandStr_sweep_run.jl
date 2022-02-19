using NPZ
using ProgressBars

include("../SMFGS/SMFGS.jl")
using .SMFGS

# figconv(lin)
filename = "../data/CMF_TandStr_sweep_figconvlin.npz"
T = GeomRange(1e-3, 1e+3, 6*25)
T = LinRange(1e-4, 20, 500)
ζ = [ 0.3 ; 0.6 ; 1.2 ; 1.0 ; 2.0 ; 4.0 ; 2.5 ; 5.0 ; 10.0 ; 20.0 ; 40.0 ]
S0 = 1//2

# figcss
#filename = "../data/CMF_TandStr_sweep_figcss.npz"
#T = GeomRange(1e-3, 1e+4, 7*25)
#ζ = [ 0.01 ; 0.02 ; 0.04 ; 0.08 ;
#      0.015 ; 0.03 ; 0.06 ;
#      0.1 ; 0.25 ; 0.5 ;
#      1.0 ; 2.0 ; 4.0 ;
#      0.75 ; 1.5 ; 3.0 ;
#      7.0 ; 14.0 ; 28.0 ;
#      500.0 ; 1000.0 ; 2000.0 ]
#S0 = 1//2

# hdapt2
#filename = "../data/CMF_TandStr_sweep_hdapt2.npz"
#T = GeomRange(1e-3, 1e+4, 7*15)
#ζ = GeomRange(1e-6, 1e+5, 11*15)
#S0 = 1000

ζ0 = 1
ω0, Γ = 7.0, 5.0
α = (2*ω0^2)*(ζ./ζ0)
θ0 = π/4

println("Running for")
println(min(T...), " ", max(T...), " ", length(T))
println(min(ζ...), " ", max(ζ...), " ", length(ζ))

s_cmf = zeros((length(ζ), length(T), 3))
s_cwk = zeros((length(ζ), length(T), 3))
s_cgs = zeros((length(T), 3))
s_cus = zeros((length(T), 3))
s_css = zeros((length(ζ), 3))
s_gnd = zeros((length(ζ), 3))
s_qgs = zeros((length(T), 3))
s_qwk = zeros((length(ζ), length(T), 3))
s_qgv = zeros((length(ζ), length(T), 3))

idx = [id for id in Iterators.product(1:length(ζ), 1:length(T))]
Threads.@threads for d in ProgressBar(idx)
  n,m = d
  s_cmf[n,m,:] = S_CMF(1/T[m], ζ[n], θ0)
  s_cwk[n,m,1] = Sx_CWK(1/T[m], ζ[n], θ0)
  s_cwk[n,m,3] = Sz_CWK(1/T[m], ζ[n], θ0)
  s_qgv[n,m,1] = Sx_QGV(1/T[m], ζ[n])
  s_qgv[n,m,3] = Sz_QGV(1/T[m], ζ[n])
  J = LorentzianSD(α[n], ω0, Γ)
  s_qwk[n,m,1] = Sx_QWK(1/T[m], S0, J, ζ0, θ0)
  s_qwk[n,m,3] = Sz_QWK(1/T[m], S0, J, ζ0, θ0)
end

Threads.@threads for m in 1:length(T)
  s_cgs[m,3] = Sz_CG(1/T[m])
  s_cus[m,1] = Sx_CUS(1/T[m], θ0)
  s_cus[m,3] = Sz_CUS(1/T[m], θ0)
  s_qgs[m,3] = Sz_QG(1/T[m], S0)
end

Threads.@threads for n in 1:length(ζ)
  #s_css[n,:] = S_CMF_T0(ζ[n], θ0)
  s_gnd[n,:] = S_CMF_ground(ζ[n], θ0)
end
 
npzwrite(filename,
  Dict("T" => T, "zeta" => ζ, "theta" => θ0,
       "scmf" => s_cmf, "scwk" => s_cwk, "scgs" => s_cgs, "scus" => s_cus,
       "scss" => s_css, "sgnd" => s_gnd,
       "sqgs" => s_qgs, "sqwk" => s_qwk, "sqgv" => s_qgv));

