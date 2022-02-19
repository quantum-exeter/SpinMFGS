#### Quantum MFGS weak (2nd order) coupling ####

## ... ##
function Σ(J::SpectralDensity; ε=1e-10)
  I = ω -> 0.5*sd(J,ω)*ω/(ω^2 - 1)
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

## ... ##
function Σ′(J::SpectralDensity; ε=1e-10)
  I = ω -> sd(J,ω)*ω/(ω^2 - 1)^2
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

## ... ##
function Δ(J::SpectralDensity, β, S0::Rational; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I = ω -> 0.5*sd(J,ω)*coth(β*ω*invtwoS0)/(ω^2 - 1)
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

Δ(J::SpectralDensity, β, S0::Int; ε=1e-10) = Δ(J, β, S0//1, ε=ε)

## ... ##
function Δ′(J::SpectralDensity, β, S0::Rational; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I = ω -> 0.5*sd(J,ω)*(ω^2 + 1)*coth(β*ω*invtwoS0)/(ω^2 - 1)^2
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

Δ′(J::SpectralDensity, β, S0::Int; ε=1e-10) = Δ'(J, β, S0//1, ε=ε)

## Quantum Sz weak coupling ##
## using QG expectation values ##
function Sz_QWK(β, S0::Rational, J::SpectralDensity, ζ0, θ0; ε=1e-10)
  Sz   = Sz_QG(β, S0)
  Szz  = Szz_QG(β, S0)
  Szzz = Szzz_QG(β, S0)

  Σ1 = ζ0*Σ(J, ε=ε)
  Σ2 = ζ0*Σ′(J, ε=ε)
  Δ1 = ζ0*Δ(J, β, S0, ε=ε)/S0
  Δ2 = ζ0*Δ′(J, β, S0, ε=ε)/S0
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

Sz_QWK(β, S0::Int, J::SpectralDensity, ζ0, θ0; ε=1e-10) = Sz_QWK(β, S0//1, J, ζ0, θ0, ε=ε)

## Quantum Sx weak coupling ##
## using QG expectation values ##
function Sx_QWK(β, S0::Rational, J::SpectralDensity, ζ0, θ0; ε=1e-10)
  Sz   = Sz_QG(β, S0)
  Szz  = Szz_QG(β, S0)

  Σ1 = ζ0*Σ(J, ε=ε)
  Δ1 = ζ0*Δ(J, β, S0, ε=ε)/S0
  Q  = ζ0*reorgenergy(J)

  invS0 = 1/S0

  prefact = -sin(2*θ0)
  term1 = (1 + invS0 - Szz)*Σ1
  term2 = Sz*Δ1
  term3 = -Szz*Q
  Sx = prefact*(term1 + term2 + term3)
  return Sx
end

Sx_QWK(β, S0::Int, J::SpectralDensity, ζ0, θ0; ε=1e-10) = Sx_QWK(β, S0//1, J, ζ0, θ0, ε=ε)
