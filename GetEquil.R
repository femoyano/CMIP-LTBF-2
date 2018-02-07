### Equilibrium solutions
### For: dec_fun = 'MM', upt_fun = '1st', diff = 'hama'

GetEquil <- function(pars, site_data = NULL) {
  pars <- c(pars, site_data)
  initial_state <- with(as.list(pars), {

    # Temp function
    TempRespEq <- function(k_ref, T, T_ref, E, R) {
      k_ref * exp(-E/R * (1/T - 1/T_ref))
    }

    # Site data
    moist = eq_moist
    temp = eq_temp
    I_sl = eq_litt

    # Calculate end variables
    K_D = TempRespEq(K_D_ref, temp, T_ref, E_K, R)
    V_D = TempRespEq(V_D_ref, temp, T_ref, E_V, R)
    r_md = TempRespEq(r_md_ref, temp, T_ref, E_m, R)
    r_mr = TempRespEq(r_mr_ref, temp, T_ref, E_r, R)
    k_ads <- TempRespEq(k_ads_ref, temp, T_ref, E_V, R)
    k_des <- TempRespEq(k_des_ref, temp, T_ref, E_V, R)
    Amax = 200 * (100 * clay)^0.6 / 1000000 * bd * depth * Amax_mod # [kg] max adsorption to minerals (Mayes et al. 2012)

    C_P <-
      -K_D*depth*(2*I_sl*r_md - V_D*f_ug*min_md - 2*V_D*f_ug*r_mr + V_D*min_md + V_D*r_mr + V_D*sqrt(
        -4*I_sl*f_ug**2*r_md + 4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
          2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2) + 2*min_md*r_mr + 2*r_mr**2)/
      (2*I_sl*r_md + 2*V_D**2*f_ug**2 - 2*V_D**2*f_ug - 2*V_D*f_ug*min_md - 4*V_D*f_ug*r_mr +
         2*V_D*min_md + 2*V_D*r_mr + 2*min_md*r_mr + 2*r_mr**2)

    C_M <-
      (-f_ug*min_md + min_md + r_mr - sqrt(
        -4*I_sl*f_ug**2*r_md + 4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
          2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2))/(2*r_md*(f_ug - 1))

    C_A  <-
      Amax*C_P*k_ads/(C_P*k_ads + k_des)

    return(c(C_P = C_P, C_A = C_A, C_M = C_M))
  })
}
