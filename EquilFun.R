### Equilibrium solutions

EquilFun <- function(pars = pars, site_data = site_data) {
# browser()
  pars_new <- c('eq_litt' <- eq_litt, 'V_D_ref' <- V_D_ref, 'r_md_ref' <- r_md_ref)
  pars <- c(pars, site_data)
  for(n in names(pars_new)) pars[[n]] <- pars_new[[n]]
  pars <- ParsCalc(pars)

  eq <- GetEquil(pars = pars)
  C_mod <- (eq[['C_P']] + eq[['C_A']] + eq[['C_M']]) * 10000 / 1000 # total converted to tons C per ha
}
