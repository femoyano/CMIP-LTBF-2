# Cost All function

CostAll <- function(pars_optim_all, pars_default = pars_default, site_data = site_data, C_obs = C_obs, input_spin, input_trans) {

  res_all <- NA


  for (i in 1:7) {

  site_data_ind <- site_data[i,]
  input_spin_ind <- input_spin[input_spin$climate == site_data_ind$climate,]
  input_trans_ind <- input_trans[input_trans$climate == site_data_ind$climate,]

  input_spin_ind$litt <- site_data_ind[['eq_litt']] # in kgC m-2 h-1
  input_spin_ind <- input_spin_ind[input_spin_ind$year < 11, ]
  input_spin_ind <- MonthlyInput(input_spin_ind)
  input_trans_ind <- MonthlyInput(input_trans_ind)
  pars_default["eq_litt"] <- site_data_ind[['eq_litt']] # in kgC m-2 h-1

  observ <- obs[obs$site == site_data_ind$site,]


  # run the spinup and get the last values for starting transient simulation
  out.sp <- as.data.frame(
    StartRun(input = input_spin_ind, site_data = site_data_ind, pars_default = pars_default, pars_new = pars_optim_all, initial_state = 'equil',
             spin.years = spin.years, spinup = 1)
  )

  n <- nrow(out.sp)
  initial <- c(C_P = out.sp$C_P[n], C_A = out.sp$C_A[n], C_M = out.sp$C_M[n])
  # browser()
  # run the transient simulation
  out <- StartRun(input = input_trans_ind, site_data = site_data_ind, pars_default = pars_default, pars_new = pars_optim_all, initial_state = initial, tsave = month)
  out <- out[1:nrow(out)%%12==0, ] # get values every 12 months, i.e. one per year
  years <- c(0:(nrow(out)-1))
  out <- out[years %in% observ$time, ]
  out <- out * 10000 / 1000  # conversion to tons C per ha
  soc <- apply(out[,-1], 1, sum) # out[['C_P']] + out[['C_A']] + out[['C_D']] + out[['C_M']]
  res <- observ$soc.t.ha - soc
  sqres <- mean(res^2)
  if(is.na(sqres)) sqres <- 9999
  res_all <- c(res_all, res)
  print(c(sqres, pars_optim_all))

  }
  return(res_all[-1])


  }
