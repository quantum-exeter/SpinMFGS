#### Adimentionalization of various parameters ####
#### Most functions in the module are desined ####
#### to be used only with these (non) units ####

## Adimentionalize the inverse temperature  ##
function β_admin(β, S0, ωL)
  return β*S0*ωL
end

## Adimentionalize temperature ##
function T_admin(T, S0, ωL, kB=1.0)
  return kB*T/(S0*ωL)
end

## Adimentionalize spin ##
function S_admin(S, ℏ=1.0)
  return S/ℏ
end

## Adimentionalize time ##
function t_admin(t, ωL)
  return t*ωL
end

## Adimentionalize frequency ##
function ω_admin(ω, ωL)
  return ω/ωL
end
