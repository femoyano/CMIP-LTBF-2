### Equilibrium solutions

CostEquil <- function(par, pars_default = pars_default, C_obs = C_obs, site_data = site_data) {
# browser()
  # Add or replace parameters from the list of optimized parameters
  # ----------------------
  pars_default <- c(pars_default, site_data)
  for(n in names(par)) pars_default[[n]] <- par[[n]]
  pars <- ParsCalc(pars_default)

  eq <- GetEquil(pars = pars)
  C_mod <- (eq[['C_P']] + eq[['C_A']] + eq[['C_M']]) * 10000 / 1000 # total converted to tons C per ha
  sqres <- (C_obs - C_mod)^2
  if(sqres > 9999 | eq[['C_P']] < 0) sqres <- 9999
  print(sqres)
  return(sqres)
}
