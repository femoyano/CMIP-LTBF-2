### Equilibrium solutions

CostEquil <- function(pars_new, pars = pars, C_obs = C_obs, site_data = site_data) {

  # Add or replace parameters from the list of optimized parameters
  # ----------------------
  pars <- c(pars, site_data)
  for(n in names(pars_new)) pars[[n]] <- pars_new[[n]]
  pars <- ParsCalc(pars)

  eq <- GetEquil(pars = pars)
  C_mod <- (eq[['C_P']] + eq[['C_A']]) * 10000 / 1000 # total converted to tons C per ha
  sqres <- (C_obs - C_mod)^2
  # print(sqres)
  return(sqres)
}
