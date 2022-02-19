#### Classical MFGS weak (2nd order) coupling ####

## Classical Sz weak coupling ##
## using CG expectation values ##
function Sz_CWK_cgdep(β, ζ, θ0)
  Sz   = Sz_CG(β)
  Szz  = Szz_CG(β)
  Szzz = Szzz_CG(β)
  return Sz + 0.5*β*ζ*(3*(cos(θ0))^2 - 1)*(Szzz -Sz*Szz)
end

## Classical Sx weak coupling ##
## using CG expectation values ##
function Sx_CWK_cgdep(β, ζ, θ0)
  Sz   = Sz_CG(β)
  Szz  = Szz_CG(β)
  Szzz = Szzz_CG(β)
  return -0.5*sin(2*θ0)*β*ζ*(Sz - Szzz)
end

## Classical Sz weak coupling ##
## using explicit expressions ##
function Sz_CWK_explicit(β, ζ, θ0)
  term1 = -coth(β) + 1/β
  prefact2 = -ζ*(3*(cos(θ0))^2 - 1)
  term21 = (csch(β))^2
  term22 = coth(β)/β
  term23 = -2/β^2
  Sz = term1 + prefact2*(term21 + term22 + term23)
  return Sz
end

## Classical Sx weak coupling ##
## using explicit expressions ##
function Sx_CWK_explicit(β, ζ, θ0)
  return ζ*sin(2*θ0)*(1 - 3*(coth(β)/β - 1/β^2))
end

Sz_CWK = Sz_CWK_cgdep
Sx_CWK = Sx_CWK_cgdep
