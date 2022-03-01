#### Spectral density integrals ####

## Σ SD integral used in the WK approximation ##
function Σ(J::SpectralDensity)
  I(ω) = 0.5*sd(J,ω)*ω/(ω + 1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end

## Σ′ SD integral used in the WK approximation ##
function Σ′(J::SpectralDensity)
  I(ω) = sd(J,ω)*ω/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end

## Δ SD integral used in the WK approximation ##
function Δ(J::SpectralDensity, β, S0::Union{Int,Rational})
  invtwoS0 = 1/(2*S0)
  I(ω) = 0.5*sdinvω(J,ω)*xcoth(β*ω*invtwoS0)/(β*invtwoS0)/(ω + 1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end

## Δ′ SD integral used in the WK approximation ##
function Δ′(J::SpectralDensity, β, S0::Union{Int,Rational})
  invtwoS0 = 1/(2*S0)
  I(ω) = 0.5*sdinvω(J,ω)*(ω^2 + 1)*xcoth(β*ω*invtwoS0)/(β*invtwoS0)/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end
