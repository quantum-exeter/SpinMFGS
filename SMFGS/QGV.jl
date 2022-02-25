#### Gelzinis and Valkunas MFGS approximation for spin 1/2 ####

## Sz expectation value from the GV approximation ##
function Sz_QGV(β, ζ, θ0)
  g = β*ζ/3
  η = sqrt(cos(θ0)^2 + exp(-2*g)*sin(θ0)^2)
  z = cos(θ0)^2 + exp(-g)*sin(θ0)^2
  return -(z/η)*tanh(β*η)
end

## Sx expectation value from the GV approximation ##
function Sx_QGV(β, ζ, θ0)
  g = β*ζ/3
  η = sqrt(cos(θ0)^2 + exp(-2*g)*sin(θ0)^2)
  x = cos(θ0)*sin(θ0)*(1 - exp(-g))
  return -(x/η)*tanh(β*η)
end

## Estimated bound of validity of the GV approximation ##
function QGV_bound(ζ)
  return 2*ζ
end
