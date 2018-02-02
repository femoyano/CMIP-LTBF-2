### Equilibrium solutions
### For: dec_fun = 'MM', upt_fun = '1st', diff = 'hama'

GetEquil <- function(pars) {
  initial_state <- with(as.list(pars), {

    # Temp function
    TempRespEq <- function(k_ref, T, T_ref, E, R) {
      k_ref * exp(-E/R * (1/T - 1/T_ref))
    }

    # Site data
    moist = eq_moist
    temp = eq_temp
    I_sl = eq_litt
    I_ml = eq_litt/10

    # Calculate intermediate variables
    b = 2.91 + 15.9 * clay
    # k_ads = TempRespEq(k_ads_ref, temp, T_ref, E_V, R) k_des
    # = TempRespEq(k_des_ref, temp, T_ref, E_V, R)
    psi_sat = exp(6.5 - 1.3 * sand)/1000
    fc = ps * (psi_sat/psi_fc)^(1/b)

    if (moist <= Dth) {g_sm <- 0} else { g_sm <- (ps - Dth)^n * ((moist - Dth) / (ps - Dth))^m }

    # Calculate end variables
    K_D = TempRespEq(K_D_ref, temp, T_ref, E_K, R)
    V_D = TempRespEq(V_D_ref, temp, T_ref, E_V, R)
    V_U = TempRespEq(V_U_ref, temp, T_ref, E_V, R)
    r_md = TempRespEq(r_md_ref, temp, T_ref, E_m, R)
    r_mr = TempRespEq(r_mr_ref, temp, T_ref, E_r, R)
    k_ads <- TempRespEq(k_ads_ref, temp, T_ref, E_V, R)
    k_des <- TempRespEq(k_des_ref, temp, T_ref, E_V, R)
    g = g_sm
    Amax = 200 * (100 * clay)^0.6 / 1000000 * bd * depth * Amax_mod # [kg] max adsorption to minerals (Mayes et al. 2012)

    C_P <-
      -K_D*depth*(2*I_ml**2*f_ug**2*r_md + 4*I_ml*I_sl*f_ug*r_md - I_ml*V_D*f_ug**2*min_md +
                    I_ml*V_D*f_ug*min_md - I_ml*V_D*f_ug*r_mr + I_ml*V_D*f_ug*sqrt(
                      -4*I_ml*f_ug**2*r_md + 4*I_ml*f_ug*r_md - 4*I_sl*f_ug**2*r_md +
                        4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
                        2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2) +
                    2*I_ml*f_ug*min_md*r_mr + 2*I_sl**2*r_md - I_sl*V_D*f_ug*min_md -
                    2*I_sl*V_D*f_ug*r_mr + I_sl*V_D*min_md + I_sl*V_D*r_mr + I_sl*V_D*sqrt(
                      -4*I_ml*f_ug**2*r_md + 4*I_ml*f_ug*r_md - 4*I_sl*f_ug**2*r_md +
                        4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
                        2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2) +
                    2*I_sl*min_md*r_mr + 2*I_sl*r_mr**2)/(
                      2*I_ml**2*f_ug**2*r_md + 4*I_ml*I_sl*f_ug*r_md + 2*I_ml*V_D**2*f_ug**2 -
                        2*I_ml*V_D**2*f_ug - 2*I_ml*V_D*f_ug**2*min_md + 2*I_ml*V_D*f_ug*min_md -
                        2*I_ml*V_D*f_ug*r_mr + 2*I_ml*f_ug*min_md*r_mr + 2*I_sl**2*r_md +
                        2*I_sl*V_D**2*f_ug**2 - 2*I_sl*V_D**2*f_ug - 2*I_sl*V_D*f_ug*min_md -
                        4*I_sl*V_D*f_ug*r_mr + 2*I_sl*V_D*min_md + 2*I_sl*V_D*r_mr +
                        2*I_sl*min_md*r_mr + 2*I_sl*r_mr**2)

    C_D <-
      (-I_ml*f_ug*r_md + I_ml*r_md - I_sl*f_ug*r_md + I_sl*r_md - f_ug*min_md*r_mr/2 +
         min_md*r_mr/2 + r_mr**2/2 - r_mr*sqrt(
           -4*I_ml*f_ug**2*r_md + 4*I_ml*f_ug*r_md - 4*I_sl*f_ug**2*r_md +
             4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
             2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2)/2)/
      (V_U*g*r_md*(f_ug**2 - 2*f_ug + 1))

    C_A <-
      Amax*C_D*k_ads/(C_D*k_ads + k_des)

    C_M <-
      (-f_ug*min_md + min_md + r_mr - sqrt(
        -4*I_ml*f_ug**2*r_md + 4*I_ml*f_ug*r_md - 4*I_sl*f_ug**2*r_md +
          4*I_sl*f_ug*r_md + f_ug**2*min_md**2 - 2*f_ug*min_md**2 -
          2*f_ug*min_md*r_mr + min_md**2 + 2*min_md*r_mr + r_mr**2))/(2*r_md*(f_ug - 1))

    return(c(C_P = C_P, C_D = C_D, C_A = C_A, C_M = C_M))
  })
}
