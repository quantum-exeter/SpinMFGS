#### Spectral densities types ####

abstract type SpectralDensity end

function sd(J::SpectralDensity, ω)
  return sdinvω(J,ω)*ω
end

function sdinvω(J::SpectralDensity, ω)
  return sd(J,ω)/ω
end

function reorgenergy(J::SpectralDensity)
  return quadgk(ω -> sdinvω(J,ω), 0.0, Inf)[1]
end

## Lorentzian Spectral Density ##

struct LorentzianSD{T<:Real} <: SpectralDensity
  α::T
  ω0::T
  Γ::T
end

function sdinvω(J::LorentzianSD, ω)
  return (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)
end

reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2

## Ohmic Spectral Density ##

struct OhmicSD{T<:Real} <: SpectralDensity
  α::T
  ωc::T
end

function sdinvω(J::OhmicSD, ω)
  return J.α*exp(-ω/J.ωc)
end

reorgenergy(J::OhmicSD) = J.α*J.ωc

## Ohmic Spectral Density ##

struct PolySD{T<:Real} <: SpectralDensity
  n::Int
  α::T
  ωc::T
end

function sdinvω(J::PolySD, ω)
  return J.α*((ω/J.ωc)^(J.n-1))*exp(-ω/J.ωc)
end

reorgenergy(J::PolySD) = J.α*J.ωc*factorial(J.n-1)

## Drude-Lorentz Spectral Density ##

struct DrudeLorentzSD{T<:Real} <: SpectralDensity
    γ::T
    ωD::T
end

function sdinvω(J::DrudeLorentzSD, ω)
  return (2/π)*J.γ*J.ωD/(ω^2 + J.ωD^2)
end

reorgenergy(J::DrudeLorentzSD) = J.γ
