#### Model.R ================================================================

# must be defined as: func <- function(t, y, parms,...) for use with ode
Model_desolve <- function(t, initial_state, pars,
                          Approx_I_sl = Approx_I_sl,
                          Approx_I_ml = Approx_I_ml,
                          Approx_temp = Approx_temp,
                          Approx_moist = Approx_moist,
                          spinup, end) {

  with(as.list(c(initial_state, pars)), {

    # set time used for interpolating input data.
    t_i <- t
    if (spinup) t_i <- t%%end  # this causes spinups to repeat the input data

    # Calculate the input and forcing at time t
    I_sl <- Approx_I_sl(t_i)
    I_ml <- Approx_I_ml(t_i)
    temp <- Approx_temp(t_i)
    moist <- Approx_moist(t_i)

    # Calculate temporally changing variables
    K_D   <- TempRespEq(K_D_ref, temp, T_ref, E_K, R)
    if(upt_fun == "MM") K_U   <- TempRespEq(K_U_ref, temp, T_ref, E_K, R)
    V_D   <- TempRespEq(V_D_ref, temp, T_ref, E_V, R)
    V_U <- TempRespEq(V_U_ref, temp, T_ref, E_V, R)
    r_md  <- TempRespEq(r_md_ref, temp, T_ref, E_m , R)
    r_mr  <- TempRespEq(r_mr_ref, temp, T_ref, E_r , R)
    k_ads <- TempRespEq(k_ads_ref, temp, T_ref, E_ad, R)
    k_des <- TempRespEq(k_des_ref, temp, T_ref, E_ad, R)

    ## Diffusion calculations  --------------------------------------
    if (moist <= Dth) {g_sm <- 0} else
      g_sm <- (ps - Dth)^n * ((moist - Dth) / (ps - Dth))^m
    g <- g_0 * g_sm

    C_D_dif <- g * C_D

    ### Calculate all fluxes ------
    # Input rate
    F_slcp <- I_sl
    F_mlcd <- I_ml

    # Decomposition rate
    if(dec_fun == "MM") {
      F_cpcd <- ReactionMM(C_P, C_M, V_D, K_D, depth)
    }
    if(dec_fun == "2nd") {
      F_cpcd <- Reaction2nd(C_P, C_M, V_D, depth)
    }
    if(dec_fun == "1st") {
      F_cpcd <- Reaction1st(C_P, V_D)
    }

    # Calculate the uptake flux
    if(upt_fun == "MM") {
      Ucd <- ReactionMM(C_D_dif, C_M, V_U, K_U, depth)
    }
    if(upt_fun == "2nd") {
      Ucd <- Reaction2nd(C_D_dif, C_M, V_U, depth)
    }
    if(upt_fun == "1st") {
      Ucd <- Reaction1st(C_D_dif, V_U)
    }

	# Adsorption/desorption
    if(flag_ads) {
      F_cdca  <- F_adsorp(C_D, C_A, Amax, k_ads)
      F_cacd  <- F_desorp(C_A, k_des)
    } else {
      F_cdca <- 0
      F_cacd <- 0
    }

    # Microbial growth, mortality, respiration and enzyme production
    F_cdcm <- Ucd * f_ug
    F_cdcr <- Ucd * (1 - f_ug)

    ## linear rmd adjustment
    F_cmcp <- C_M * (C_M * r_md + min_md)
    F_cmcd <- 0
    F_cmcr <- C_M * r_mr

    ## Rate of change calculation for state variables ---------------
    dC_P  <- F_slcp + F_cmcp - F_cpcd
    dC_D  <- F_mlcd + F_cpcd  + F_cmcd - F_cdcm - F_cdcr + F_cacd - F_cdca
    dC_A  <- F_cdca - F_cacd
    dC_M  <- F_cdcm - F_cmcp - F_cmcr - F_cmcd
    dC_Rg <- F_cdcr
    dC_Rm <- F_cmcr

    return(list(c(dC_P, dC_D, dC_A, dC_M)))

  })  # end of with(...

}  # end of Model
