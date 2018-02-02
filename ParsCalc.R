# ParsReplace.R
# Adds or replaces parameters from the list of optimized parameters ----------------------

ParsCalc <- function(pars) {
  new <- with(as.list(pars), {
    # psi_sat <- exp(6.5 - 1.3 * sand) / 1000 # [kPa] saturation water potential (Cosby et al. 1984 after converting their data from cm H2O to Pa)
    # b       <- 2.91 + 15.9 * clay  # [] b parameter (Campbell 1974) as in Cosby  et al. 1984 - Alternatively: obtain from land model.
    # # Rth     = ps * (psi_sat / psi_Rth)^(1 / b),  # [m3 m-3] Threshold relative water content for mic. respiration (water retention formula from Campbell 1984)
    # fc      <- ps * (psi_sat / psi_fc)^(1 / b)  # [m3 m-3] Field capacity relative water content (water retention formula from Campbell 1984) - Alternatively: obtain from land model.
    # Mad     <- ps * (psi_sat / psi_ad)^(1 / b)  # [m3 m-3] Water content at air dryness
    Amax    <- 200 * (100 * clay)^0.6 / 1000000 * bd * depth * Amax_mod  # [kg] max adsorption to minerals (Mayes et al. 2012)
    c(Amax = Amax)
  })
  return(c(pars, new))
}
