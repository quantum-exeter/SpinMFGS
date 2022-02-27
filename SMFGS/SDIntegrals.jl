#### Spectral density integrals ####

## Σ SD integral used in the WK approximation ##
function Σ(J::SpectralDensity; ε=1e-10)
  I = ω -> 0.5*sd(J,ω)*ω/(ω^2 - 1)
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

## Σ′ SD integral used in the WK approximation ##
function Σ′(J::SpectralDensity; ε=1e-10)
  I = ω -> sd(J,ω)*ω/(ω^2 - 1)^2
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

## Δ SD integral used in the WK approximation ##
function Δ(J::SpectralDensity, β, S0::Rational; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I = ω -> 0.5*sd(J,ω)*coth(β*ω*invtwoS0)/(ω^2 - 1)
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

Δ(J::SpectralDensity, β, S0::Int; ε=1e-10) = Δ(J, β, S0//1, ε=ε)

## Δ′ SD integral used in the WK approximation ##
function Δ′(J::SpectralDensity, β, S0::Rational; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I = ω -> 0.5*sd(J,ω)*(ω^2 + 1)*coth(β*ω*invtwoS0)/(ω^2 - 1)^2
  Il = quadgk(I, 0.0, 1.0 - ε)[1]
  Ir = quadgk(I, 1.0 + ε, Inf)[1]
  return Il + Ir
end

Δ′(J::SpectralDensity, β, S0::Int; ε=1e-10) = Δ'(J, β, S0//1, ε=ε)
