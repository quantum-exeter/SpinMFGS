#using Plots
using NPZ
using ProgressBars

include("../SMFGS/SMFGS.jl")
using .SMFGS

#filename="../data/QRCTest.npz"
#T = GeomRange(1e-3, 1e+3, 6*3)
#ζ = GeomRange(1e-3, 1e+4, 7*3)
#S0 = [ 1//2 ; 2//2 ; 3//2 ; 4//2 ; 5//2 ; 8//2 ; 12//2 ; 20//2 ; 30//2 ; 40//2 ; 80//2 ]
#Ncutoff = [ 2 ; 3 ; 5 ; 10 ; 20 ; 30 ; 40 ; 50 ; 70 ; 80 ]

filename="../data/QRCTest_vShalf.npz"
T = GeomRange(1e+0, 1e+3, 3*3)
ζ = GeomRange(1e+0, 1e+4, 4*4)
S0 = [ 1//2 ]
Ncutoff = 2:30:1230

ζ0 = 1
ω0, Γ = 7.0, 0.01
α = (2*ω0^2)*(ζ./ζ0)
θ0 = π/4

println("Running for")
println(min(T...), " ", max(T...), " ", length(T))
println(S0)
println(Ncutoff)

s_cmf = zeros((length(ζ), length(T), 3))
s_cwk = zeros((length(ζ), length(T), 3))
s_cgs = zeros((length(T), 3))
s_cus = zeros((length(T), 3))
s_css = zeros((length(ζ), 3))
s_gnd = zeros((length(ζ), 3))
s_qwk = zeros((length(S0), length(ζ), length(T), 3))
s_qgs = zeros((length(S0), length(T), 3))
s_qrc = zeros((length(Ncutoff), length(S0), length(ζ), length(T), 3))
n_qrc = zeros((length(Ncutoff), length(S0), length(ζ), length(T)))
s_qgv = zeros((length(ζ), length(T), 3))

idx = [id for id in Iterators.product(1:length(ζ), 1:length(T))]
Threads.@threads for d in ProgressBar(idx)
  n,m = d
  #s_cmf[n,m,:] = S_CMF(1/T[m], ζ[n], θ0)
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

for j in 1:length(S0)
  for k in ProgressBar(1:length(Ncutoff))
    operators = precompute_operators_QRC(S0[j], Ncutoff[k])
    Threads.@threads for d in idx
      n,m = d
      v = S_QRC(1/T[m], S0[j], α[n], ω0, θ0; operators=operators, computen=true)
      s_qrc[k,j,n,m,:] = v[1:3]
      n_qrc[k,j,n,m] = v[4]
    end
  end
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
  #s_css[n,:] = S_CMF_T0(ζ[n], θ0)
  s_gnd[n,:] = S_CMF_ground(ζ[n], θ0)
end
 
npzwrite(filename,
  Dict("S0" => Float64.(S0), "T" => T, "zeta" => ζ, "zeta0" => ζ0,
       "omega0" => ω0, "Gamma" => Γ, "alpha" => α, "theta" => θ0,
       "sqwk" => s_qwk, "sqgs" => s_qgs, "sqrc" => s_qrc, "nqrc" => n_qrc, "Ncutoff" => Ncutoff,
       "scmf" => s_cmf, "scwk" => s_cwk, "scgs" => s_cgs, "scus" => s_cus,
       "scss" => s_css, "sgnd" => s_gnd));

