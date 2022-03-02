#### Quantum Reaction-Coordinate Mapping for LorentzianSDs ####

## Quantum spin operators ##
function sz_Op(S0::Union{Int,Rational})
  Diagonal([m/S0 for m in -S0:1:S0])
end

function sx_Op(S0::Union{Int,Rational})
  SymTridiagonal(zeros(Int(2*S0+1)),
                 [sqrt(S0*(S0+1)-m*(m+1))/2/S0 for m in -S0:1:S0-1])
end

function sy_Op(S0::Union{Int,Rational})
  Tridiagonal([sqrt(S0*(S0+1)-m*(m-1))/2im/S0 for m in -S0+1:1:S0],
               zeros(Int(2*S0+1)),
              [-sqrt(S0*(S0+1)-m*(m+1))/2im/S0 for m in -S0:1:S0-1])
end

## Quantum harmonic oscillator operators ##
function N_Op(Ncutoff::Int)
  Diagonal(0:Ncutoff)
end

function X_Op(Ncutoff::Int)
  SymTridiagonal(zeros(Ncutoff+1),
                [sqrt(n+1) for n in 0:Ncutoff-1])
end

## Pre-compute all operators ##
function precompute_operators_QRC(S0::Union{Int,Rational}, Ncutoff::Int; computey=false)
  sz, sx, n, x = sz_Op(S0), sx_Op(S0), N_Op(Ncutoff), X_Op(Ncutoff)
  szx = kronecker(sz,x)
  sxx = kronecker(sx,x)
  szI = kronecker(sz,I(Ncutoff+1))
  sxI = kronecker(sx,I(Ncutoff+1))
  In = kronecker(I(Int(2*S0+1)),n)
  if computey
    syI = kronecker(sy_Op(S0),I(Ncutoff+1))
    return [szI, sxI, szx, sxx, In, syI]
  else
    return [szI, sxI, szx, sxx, In]
  end
end

## Type selection ##
function seltype_QRC(β, E)
  mexponent = β*E
  if mexponent < 1e2
    return Float64
  elseif mexponent < 4e3
    return Float128
  else
    return BigFloat
  end
end

## Reaction Coordinate Hamiltonian ##
function H_QRC(S0::Union{Int,Rational}, α, ω0, θ0; Ncutoff::Int=20, operators=nothing)
  if isnothing(operators)
    operators = precompute_operators_QRC(S0, Ncutoff; computey=false)
  end
  λ = sqrt(α/ω0/S0/2)
  return Symmetric(operators[1] + λ*(cos(θ0)*operators[3] + sin(θ0)*operators[4]) + (ω0/S0)*operators[5])
end

H_QRC(S0::Union{Int,Rational}, J::LorentzianSD, θ0; Ncutoff::Int=20, operators=nothing) = H_QRC(S0, J.α, J.ω0, θ0; Ncutoff=Ncutoff, operators=operators)

## Spin expectation values from the Reaction Coordinate thermal state ##
function S_QRC(β, S0::Union{Int,Rational}, α, ω0, θ0; Ncutoff::Int=20, operators=nothing, computey=false, computen=false, method="eigen")
  if isnothing(operators)
    operators = precompute_operators_QRC(S0, Ncutoff; computey=computey)
  end
  Hrc = H_QRC(S0, α, ω0, θ0; Ncutoff=Ncutoff, operators=operators)
  if method == "schur"
    F = schur(Hrc)
  elseif method == "eigen"
    F = eigen(Hrc)
  else
    throw(DomainError("'method' argument of S_QRC must be one of: ['schur', 'eigen']"))
  end
  E, B_E = F.values, F.vectors
  T = seltype_QRC(β, abs(minimum(E)))
  weights = [exp(-T(β*E[k])) for k in 1:length(E)]
  Z = sum(weights)
  sz_vals = [real(adjoint(B_E[:,k])*operators[1]*B_E[:,k]) for k in 1:length(E)]
  sx_vals = [real(adjoint(B_E[:,k])*operators[2]*B_E[:,k]) for k in 1:length(E)]
  szavg = sum([sz_vals[k]*weights[k]/Z for k in 1:length(E)])
  sxavg = sum([sx_vals[k]*weights[k]/Z for k in 1:length(E)])
  if computey
    sy_vals = [real(adjoint(B_E[:,k])*operators[6]*B_E[:,k]) for k in 1:length(E)]
    syavg = sum([sy_vals[k]*weights[k]/Z for k in 1:length(E)])
  else
    syavg = 0.0
  end
  if computen
    n_vals = [real(adjoint(B_E[:,k])*operators[5]*B_E[:,k]) for k in 1:length(E)]
    navg = sum([n_vals[k]*weights[k]/Z for k in 1:length(E)])
    return Float64[sxavg ; syavg ; szavg ; navg]
  end
  return Float64[sxavg ; syavg ; szavg]
end

function S_QRC(β, S0::Union{Int,Rational}, J::LorentzianSD, θ0; Ncutoff::Int=20, operators=nothing, computey=false, computen=false, method="eigen")
  return S_QRC(β, S0, J.α, J.ω0, θ0; Ncutoff=Ncutoff, operators=operators, computey=computey, computen=computen, method=method)
end
