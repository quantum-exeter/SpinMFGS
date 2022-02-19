#### Classical Gibbs state ###

## Classical Gibbs partition function ##
## using numerical integration ##
function Z_CG_numeric(β)
  I(θ) = sin(θ)*exp(-β*cos(θ))
  Ieval = quadgk(I, 0.0, π)[1]
  return Ieval
end

## Classical Gibbs Sz^n expectation value ##
## using numerical integration ##
function Szn_CG_numeric(n::Int, β)
  Z = Z_CG_numeric(β)
  I(θ) = sin(θ)*(cos(θ)^n)*exp(-β*cos(θ))
  Ieval = quadgk(I, 0.0, π)[1]
  return Ieval/Z
end

## Classical Gibbs Sz expectation value ##
## using numerical integration ##
Sz_CG_numeric(β) = Szn_CG_numeric(1, β)

## Classical Gibbs Sz^2 expectation value ##
## using numerical integration ##
Szz_CG_numeric(β) = Szn_CG_numeric(2, β)

## Classical Gibbs Sz^3 expectation value ##
## using numerical integration ##
Szzz_CG_numeric(β) = Szn_CG_numeric(3, β)

## Classical Gibbs partition function ##
## using analytical expressions ##
function Z_CG_analytic(β)
  return sinhc(β)
end

## Classical Gibbs Sz expectation value ##
## using analytical expressions ##
function Sz_CG_analytic(β)
  return -cothminv(β)
end

## Classical Gibbs Sz^2 expectation value ##
## using analytical expressions ##
function Szz_CG_analytic(β)
  return 1 - 2*cothinvminv2(β)
end

## Classical Gibbs Sz^3 expectation value ##
## using analytical expressions ##
function Szzz_CG_analytic(β)
  return iszero(β) ? zero(β) : (2/β)*(one(β) - 3*cothinvminv2(β)) - cothminv(β)
end

Sz_CG   = Sz_CG_analytic
Szz_CG  = Szz_CG_analytic
Szzz_CG = Szzz_CG_analytic
