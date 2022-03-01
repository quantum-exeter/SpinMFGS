#### Quantum MFGS weak (2nd order) coupling ####

## Quantum Sz weak coupling ##
## using QG expectation values ##
function Sz_QWK(β, S0::Union{Int,Rational}, J::SpectralDensity, ζ0, θ0)
  Sz   = Sz_QG(β, S0)
  Szz  = Szz_QG(β, S0)
  Szzz = Szzz_QG(β, S0)

  Σ1 = ζ0*Σ(J)
  Σ2 = ζ0*Σ′(J)
  Δ1 = ζ0*Δ(J, β, S0)/S0
  Δ2 = ζ0*Δ′(J, β, S0)/S0
  Q  = ζ0*reorgenergy(J)

  invS0 = 1/S0

  prefact1 = -(sin(θ0))^2
  term11 = (1 + invS0 - Szz)*Σ2
  term12 = Sz*Δ2
  prefact2 = β*(sin(θ0))^2
  term21 = (Szz - Sz^2)*Δ1
  term22 = -(Szzz - Sz*Szz)*Σ1
  prefact3 = β*(cos(θ0))^2
  term31 = (Szzz - Sz*Szz)*Q
  Sz = Sz + prefact1*(term11 + term12) + prefact2*(term21 + term22) + prefact3*term31
  return Sz
end

## Quantum Sx weak coupling ##
## using QG expectation values ##
function Sx_QWK(β, S0::Union{Int,Rational}, J::SpectralDensity, ζ0, θ0)
  Sz   = Sz_QG(β, S0)
  Szz  = Szz_QG(β, S0)

  Σ1 = ζ0*Σ(J)
  Δ1 = ζ0*Δ(J, β, S0)/S0
  Q  = ζ0*reorgenergy(J)

  invS0 = 1/S0

  prefact = -sin(2*θ0)
  term1 = (1 + invS0 - Szz)*Σ1
  term2 = Sz*Δ1
  term3 = -Szz*Q
  Sx = prefact*(term1 + term2 + term3)
  return Sx
end
