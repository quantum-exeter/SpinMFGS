#### Spectral densities types ####

abstract type SpectralDensity end

function reorgenergy(J::SpectralDensity)
  return quadgk(ω -> sd(J,ω)/ω, 0.0, Inf)[1]
end

## Lorentzian Spectral Density ##

struct LorentzianSD{T<:Real} <: SpectralDensity
  α::T
  ω0::T
  Γ::T
end

function sd(J::LorentzianSD, ω)
  iszero(ω) ? zero(ω) : (J.α*J.Γ/π)*ω/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)
end

reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2

## Ohmic Spectral Density ##

struct OhmicSD{T<:Real} <: SpectralDensity
  α::T
  ωc::T
end

function sd(J::OhmicSD, ω)
  return J.α*ω*exp(-ω/J.ωc)
end

reorgenergy(J::OhmicSD) = J.α*J.ωc

## Ohmic Spectral Density ##

struct PolySD{T<:Real} <: SpectralDensity
  n::Int
  α::T
  ωc::T
end

function sd(J::PolySD, ω)
  return α*((ω^n)/(ωc^(n-1)))*exp(-ω/ωc)
end

reorgenergy(J::PolySD) = J.α*J.ωc*factorial(J.n-1)

## Drude-Lorentz Spectral Density ##

struct DrudeLorentzSD{T<:Real} <: SpectralDensity
    γ::T
    ωD::T
end

function sd(J::DrudeLorentzSD, ω)
  return (2/π)*J.γ*J.ωD*ω/(ω^2 + J.ωD^2)
end

reorgenergy(J::DrudeLorentzSD) = J.γ
