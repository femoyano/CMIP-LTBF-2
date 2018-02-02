# flux_functions.r

## Functions to calculate diffusion depending on options -----
get_g_sm <- function(moist, ps, Dth, n, m) {
  if (moist <= Dth) {g_sm <- 0} else
   g_sm <- (ps - Dth)^n * ((moist - Dth) / (ps - Dth))^m # n=1.5, m=2.5
}

## Reaction kinetics ---------
ReactionMM <- function(S, E, V, K, depth) {
  S <- S/(depth)
  E <- E/(depth)
  Flux <- (V * E * S)/(K + S) * depth
}

Reaction2nd <- function(S, E, V, depth) {
  S <- S/(depth)
  E <- E/(depth)
  Flux <- (V * S * E) * depth
}

Reaction1st <- function(S, V) {
  Flux <- V * S
}

F_adsorp <- function (C_D, C_A, Amax, k_ads) {
  A <- Amax - C_A
  return(C_D * A * k_ads)
}

F_desorp <- function (C_A, k_des) {
  return(C_A * k_des)
}

## Temperature response (Arrhenius)
TempRespEq <- function(k_ref, T, T_ref, E, R) {
  k_ref * exp(-E/R * (1/T - 1/T_ref))
}
