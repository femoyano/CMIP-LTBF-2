# Cost function

CostTrans <- function(pars_optim, pars_default = pars_default, site_data = site_data, C_obs = C_obs, input_spin, input_trans) {

  input_spin$litt <- c(pars_optim, pars_default, site_data)[['eq_litt']] # in kgC m-2 h-1
  # run the spinup and get the last values for starting transient simulation
  out.sp <- as.data.frame(
    StartRun(input = input_spin, site_data = site_data, pars_default = pars_default, pars_new = pars_optim, initial_state = 'equil',
               spin.years = spin.years, spinup = 1)
    )

  n <- nrow(out.sp)
  initial <- c(C_P = out.sp$C_P[n], C_M = out.sp$C_M[n], C_A = out.sp$C_A[n])

  # run the transient simulation
  out <- StartRun(input = input_trans, pars_default = pars_default, pars_new = pars_optim, initial_state = initial, tsave = month)

  out <- out[1:nrow(out)%%12==0, ] # get values every 12 months, i.e. one per year
  year <- c(0:(nrow(out)-1))
  out <- out[year %in% obs$time, ]
  out <- out * 10000 / 1000  # conversion to tons C per ha
  soc <- apply(out[,-1], 1, sum) # out[['C_P']] + out[['C_A']] + out[['C_D']] + out[['C_M']]
  res <- obs$soc.t.ha - soc
  sqres <- mean(res^2)
  if(sqres > 9999 | out[['C_P']][1] < 0) sqres <- 9999
  print(c(sqres, pars_optim))
  return(res)
}
