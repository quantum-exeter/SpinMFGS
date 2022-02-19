#### Gelzinis and Valkunas MFGS approximation for spin 1/2 ####

## Sz expectation value from the GV approximation ##
function Sz_QGV(β, ζ)
  z = 1 + exp(-β*ζ/12)
  d = 1 + exp(-β*ζ/6)
  return -z*tanh(β*sqrt(d/2))/sqrt(2*d)
end

## Sx expectation value from the GV approximation ##
function Sx_QGV(β, ζ)
  x = -expm1(-β*ζ/12)
  d = 1 + exp(-β*ζ/6)
  return x*tanh(β*sqrt(d/2))/sqrt(2*d)
end

## Estimated bound of validity of the GV approximation ##
function QGV_bound(ζ)
  return 2*ζ
end
