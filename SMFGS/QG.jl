#### Quantum Gibbs state ###

## Quantum Gibbs Sz^n expectation value ##
## using numerical integration ##
function Szn_QG_numeric(n::Int, β, S0::Rational)
  βweights = reverse([exp(-β*m/S0) for m = 0:Int(2*S0)])
  Z  = sum(βweights)
  Pz = βweights ./ Z
  Sz = sum(sort([(m/S0 - 1)^n*Pz[m+1] for m = 0:Int(2*S0)]))
  return Sz
end

Szn_QG_numeric(n::Int, β, S0::Int) = Szn_QG_numeric(n, β, S0//1)

## Quantum Gibbs Sz expectation value ##
## using numerical integration ##
Sz_QG_numeric(β, S0) = Szn_QG_numeric(1, β, S0)

## Quantum Gibbs Sz^2 expectation value ##
## using numerical integration ##
Szz_QG_numeric(β, S0) = Szn_QG_numeric(2, β, S0)
#
## Quantum Gibbs Sz^3 expectation value ##
## using numerical integration ##
Szzz_QG_numeric(β, S0) = Szn_QG_numeric(3, β, S0)

## Quantum Gibbs Sz expectation value ##
## using analytic expressions ##
function Sz_QG_analytic(β, S0::Rational)
  invtwoS0 = 1//(2*S0)
  return invtwoS0*coth(β*invtwoS0) - (1 + invtwoS0)*coth(β*(1 + invtwoS0))
end

Sz_QG_analytic(β, S0::Int) = Sz_QG_analytic(β, S0//1)

## Quantum Gibbs Sz^2 expectation value ##
## using analytic expressions ##
function Szz_QG_analytic(β, S0::Rational)
  invtwoS0 = 1//(2*S0)
  invS0 = 1//S0
  term1 = (1 + invtwoS0)^2
  term2 = -(1 + invtwoS0)*invS0*coth(β*invtwoS0)*coth(β*(1 + invtwoS0))
  term3 = (invtwoS0^2)*(2*(coth(β*invtwoS0))^2 - 1)
  Szz   = term1 + term2 + term3
  return Szz
end

Szz_QG_analytic(β, S0::Int) = Szz_QG_analytic(β, S0//1)

## Quantum Gibbs Sz^3 expectation value ##
## using analytic expressions ##
function Szzz_QG_analytic(β, S0::Rational)
  invtwoS0 = 1//(2*S0)
  invS0 = 1//S0
  term1 = 12*invS0*((1 + invtwoS0)^2)*coth(β*invtwoS0)
  term2 = (invS0*coth(β*invtwoS0))^3
  term3 = -6*(1 + invtwoS0)*((invS0*coth(β*invtwoS0))^2)*coth((1 + invtwoS0)*β)
  term4 = -2*(1 + invtwoS0)*coth((1 + invtwoS0)*β)*(4*(1 + invtwoS0)^2 + 3*(invS0*csch(β*invtwoS0))^2)
  term5 = 2.5*(invS0^3)*((csch(β*invtwoS0))^4)*sinh(β*invS0)
  Szzz  = 0.125*(term1 + term2 + term3 + term4 + term5)
  return Szzz
end

Szzz_QG_analytic(β, S0::Int) = Szzz_QG_analytic(β, S0//1)

Sz_QG   = Sz_QG_analytic
Szz_QG  = Szz_QG_analytic
Szzz_QG = Szzz_QG_analytic
