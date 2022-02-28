#### Classical MFGS ####

## Aux: type selection ##
function seltype(β, ζ)
  mexponent = β*(3 + ζ)
  if mexponent < 1e2
    return Float64
  elseif mexponent < 4e3
    return Float128
  else
    return BigFloat
  end
end

## Classical HMF ##
function H_CMF(θ, φ, ζ, θ0, B::AbstractArray=[0;0;1])
  return B[1]*xsph(θ,φ) + B[2]*ysph(θ,φ) + B[3]*zsph(θ,φ) - ζ*(zsph(θ,φ)*cos(θ0) + xsph(θ,φ)*sin(θ0))^2
end

H_CMF(S, ζ, θ0, B::AbstractArray=[0;0;1]) = H_CMF(θφsph(S)..., ζ, θ0, B)

## Classical HMF parititon function ##
function Z_CMF(β, ζ, θ0)
  T = seltype(β, ζ)
  I(θ_,φ_) = exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)
  return dblquadgk(I, [0.0, 0.0], [π, 2π])
end

## Classical HMF spin expectation values ##
## integrating ##
function S_CMF_pfint(β, ζ, θ0; Z=nothing)
  if isnothing(Z)
    Z = Z_CMF(β, ζ, θ0)
  end
  T = seltype(β, ζ)
  Ix(θ_,φ_) = xsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Iy(θ_,φ_) = ysph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Iz(θ_,φ_) = zsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Sx = Float64(dblquadgk(Ix, [0.0, 0.0], [π, 2π]))
  Sy = Float64(dblquadgk(Iy, [0.0, 0.0], [π, 2π]))
  Sz = Float64(dblquadgk(Iz, [0.0, 0.0], [π, 2π]))
  return SA[Sx ; Sy ; Sz]
end

function S_CMF_pfint_fast(β, ζ, θ0; Z=nothing)
  if isnothing(Z)
    Z = Z_CMF(β, ζ, θ0)
  end
  T = seltype(β, ζ)
  Ix(θ_,φ_) = xsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Iz(θ_,φ_) = zsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Sx = Float64(dblquadgk(Ix, [0.0, 0.0], [π, 2π]))
  Sz = Float64(dblquadgk(Iz, [0.0, 0.0], [π, 2π]))
  return [Sx ; 0.0 ; Sz]
end

function Sz_CMF_pfint(β, ζ, θ0; Z=nothing)
  if isnothing(Z)
    Z = Z_CMF(β, ζ, θ0)
  end
  T = seltype(β, ζ)
  Iz(θ_,φ_) = zsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Sz = Float64(dblquadgk(Iz, [0.0, 0.0], [π, 2π]))
  return Sz
end

function Sx_CMF_pfint(β, ζ, θ0; Z=nothing)
  if isnothing(Z)
    Z = Z_CMF(β, ζ, θ0)
  end
  T = seltype(β, ζ)
  Ix(θ_,φ_) = xsph(θ_,φ_)*exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0))))*dΩ(θ_,φ_)/Z
  Sx = Float64(dblquadgk(Ix, [0.0, 0.0], [π, 2π]))
  return Sx
end

## Classical HMF spin expectation values ##
## derivating the partition function ##
function S_CMF_pfdrv(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0,B))))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = grad(central_fdm(12, 1), logZ, [0.0 ; 0.0 ; 1.0])[1]
  s = -(1/β)*dlogZ
  return Float64.(s)
end

function S_CMF_pfdrv_fast(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0,[B[1];0.0;B[2]]))))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = grad(central_fdm(12, 1), logZ, [0.0, 1.0])[1]
  s = -(1/β)*dlogZ
  return Float64[s[1] ; 0.0 ; s[2]]
end

function Sz_CMF_pfdrv(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0,[0.0;0.0;B]))))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = central_fdm(12, 1)(logZ, 1.0)[1]
  s = -(1/β)*dlogZ
  return Float64(s)
end

function Sx_CMF_pfdrv(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-T(β*(H_CMF(θ_,φ_,ζ,θ0,[B;0.0;1.0]))))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = central_fdm(12, 1)(logZ, 0.0)[1]
  s = -(1/β)*dlogZ
  return Float64(s)
end

## Classical HMF spin expectation values ##
## derivating the partition function with automatic differentiation ##
function S_CMF_pfad(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-β*(H_CMF(θ_,φ_,ζ,θ0,B)))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = ForwardDiff.gradient(logZ, T[0.0; 0.0; 1.0])
  s = -(1/β)*dlogZ
  return Float64.(s)
end

function S_CMF_pfad_fast(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-β*(H_CMF(θ_,φ_,ζ,θ0,[B[1];0.0;B[2]])))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = ForwardDiff.gradient(logZ, T[0.0; 1.0])
  s = -(1/β)*dlogZ
  return Float64[s[1] ; 0.0 ; s[2]]
end

function Sz_CMF_pfad(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-β*(H_CMF(θ_,φ_,ζ,θ0,[0.0;0.0;B])))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = ForwardDiff.derivative(logZ, T(1.0))
  s = -(1/β)*dlogZ
  return Float64(s)
end

function Sx_CMF_pfad(β, ζ, θ0)
  T = seltype(β, ζ)
  logZ(B) = log(dblquadgk((θ_,φ_) -> exp(-β*(H_CMF(θ_,φ_,ζ,θ0,[B;0.0;1.0])))*dΩ(θ_,φ_), [0.0, 0.0], [π, 2π]))
  dlogZ = ForwardDiff.derivative(logZ, T(0.0))
  s = -(1/β)*dlogZ
  return Float64(s)
end

S_CMF = S_CMF_pfad_fast

## Classical HMF T=0 expectation values ##
function S_CMF_T0(ζ, θ0)
  if sin(θ0) ≈ 0.0
    E1 = H_CMF([0.0,0.0,1.0],ζ,θ0)
    E2 = H_CMF([0.0,0.0,-1.0],ζ,θ0)
    if E1 < E2
      return [0.0, 0.0, 1.0]
    else
      return [0.0, 0.0, -1.0]
    end
  end

  # v1
  # sz(N) = -(2*N^2*ζ*sin(θ0)^2 + N)/(2*N*ζ + 1)
  # sx(N) = N^2*ζ*sin(2*θ0)/(2*N*ζ + 1)
  # v2
  # sz(N) = N*(1 - 2*N*ζ*sin(θ0)^2)/(1 - 2*N*ζ)
  # sx(N) = 2*N^2*ζ*sin(θ0)*cos(θ0)/(1 - 2*N*ζ)
  # v3
  sx(N) = 2*N^2*ζ*sin(θ0)*cos(θ0)/(1 - 2*N*ζ)
  sz(N) = N + cot(θ0)*sx(N)

  function norm_cond(x)
    N  = x
    s  = [sx(N) ; 0.0 ; sz(N)]
    s2 = norm(s)^2
    Δ  = s2 - 1.0
    return Δ
  end

  guess = sign(ζ)
  sol_N = find_zero(norm_cond, guess)
  sol_S = [sx(sol_N) ; 0.0 ; sz(sol_N)]

  S_T0 = sol_S
  E_T0 = H_CMF(S_T0,ζ,θ0)

  tst_S = [-sx(sol_N) ; 0.0 ; sz(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [sx(sol_N) ; 0.0 ; -sz(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [-sx(sol_N) ; 0.0 ; -sz(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [sz(sol_N) ; 0.0 ; sx(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [sz(sol_N) ; 0.0 ; -sx(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [-sz(sol_N) ; 0.0 ; sx(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  tst_S = [-sz(sol_N) ; 0.0 ; -sx(sol_N)]
  tst_E = H_CMF(tst_S,ζ,θ0)
  if tst_E < E_T0
    E_T0 = tst_E
    S_T0 = tst_S
  end

  return S_T0
end

## Classical HMF gorund state ##
function S_CMF_ground(ζ, θ0; debug=false)
  # v1
  # H(θ,φ) = H_CMF(θ,φ,ζ,θ0)
  # res = optimize(x -> H(x[1], x[2]), [π, 0], Nelder-Mead())
  # res = optimize(x -> H(x[1], x[2]), [π, 0], BFGS(); autodiff=:forward)
  # res = optimize(x -> H(x[1], x[2]), [π, 0], SimulatedAnnealing())
  # res = optimize(x -> H(x[1], x[2]), [π, 0], ParticleSwarm())
  # θc, φc = Optim.minimizer(res)
  # v2
  H1(θ) = H_CMF(θ,0.0,ζ,θ0)
  H2(θ) = H_CMF(θ,π,ζ,θ0)
  res1 = optimize(x -> H1(x), 0.0, π, Brent())# autodiff=:forward)
  res2 = optimize(x -> H2(x), 0.0, π, Brent())# autodiff=:forward)
  E1 = Optim.minimum(res1)
  E2 = Optim.minimum(res1)
  if E1 < E2
    θc, φc = Optim.minimizer(res1), 0.0
  else
    θc, φc = Optim.minimizer(res2), π
  end

  return Xsph([θc ; φc])

  θg, φg = θc, φc
  Eg = H(θc, φc)

  θt, φt = π - θc, π + φc
  Et = H(θt, φt)
  if Et < Eg
    Eg = Et
    θg, φg = θt, φt
  end

  θt, φt = π - θc, φc
  Et = H(θt, φt)
  if Et < Eg
    Eg = Et
    θg, φg = θt, φt
  end

  θt, φt = θc, π + φc
  Et = H(θt, φt)
  if Et < Eg
    Eg = Et
    θg, φg = θt, φt
  end

  θt, φt = θc, -φc
  Et = H(θt, φt)
  if Et < Eg
    Eg = Et
    θg, φg = θt, φt
  end

  S_ground = Xsph([θg ; φg])
  if debug
    return S_ground, (θg, φg), res
  else
    return S_ground
  end
end

## Classical HMF probability density ##
function Pdens_CMF(θ, φ, β, ζ, θ0; Z=nothing)
  if isnothing(Z)
     Z = Z_CMF(β, ζ, θ0)
  end
  T = seltype(β, ζ)
  #return exp(-T(β*(H_CMF(θ,φ,ζ,θ0))))*dΩ(θ,φ)/Z
  return exp(-T(β*(H_CMF(θ,φ,ζ,θ0))))/Z
end
