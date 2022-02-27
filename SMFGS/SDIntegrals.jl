#### Spectral density integrals ####

## Σ SD integral used in the WK approximation ##
function Σ(J::SpectralDensity; ε=1e-10)
  I(ω) = 0.5*sd(J,ω)*ω/(ω + 1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end

## Σ′ SD integral used in the WK approximation ##
function Σ′(J::SpectralDensity; ε=1e-10)
  I(ω) = sd(J,ω)*ω/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end

## Δ SD integral used in the WK approximation ##
function Δ(J::SpectralDensity, β, S0::Union{Int,Rational}; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I(ω) = 0.5*sd(J,ω)*coth(β*ω*invtwoS0)/(ω + 1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end

## Δ′ SD integral used in the WK approximation ##
function Δ′(J::SpectralDensity, β, S0::Union{Int,Rational}; ε=1e-10)
  invtwoS0 = 1/(2*S0)
  I(ω) = 0.5*sd(J,ω)*(ω^2 + 1)*coth(β*ω*invtwoS0)/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end