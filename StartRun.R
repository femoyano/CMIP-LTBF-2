# Run_spinup.R

StartRun <- function(input, site_data, pars, pars_new = NULL, initial_state = "equil", spin.years = 0, spinup = 0, tsave = month) {

  pars <- c(pars, site_data)
  for(n in names(pars_new)) pars[[n]] <- pars_new[[n]]
  pars <- ParsCalc(pars)

  if (initial_state == "equil") {
    initial_state <- GetEquil(eq_temp = pars[['eq_temp']], eq_moist = pars[['eq_moist']],
                              eq_litt = pars[['eq_litt']], pars = pars)
  }

  end   <- tail(input$hour, 1)
  if(spinup) {finish <- spin.years * year / tstep} else {finish <- end}
  times <- seq(0, finish, by = tsave / tstep) # vector of output times

  # Define input interpolation functions
  Approx_I_sl  <- approxfun(input$hour, input$litt, method = "linear", rule = 2)
  Approx_I_ml  <- approxfun(input$hour, input$litt/10, method = "linear", rule = 2)
  Approx_temp  <- approxfun(input$hour, input$temp, method = "linear", rule = 2)
  Approx_moist <- approxfun(input$hour, input$moist, method = "linear", rule = 2)

  t0 <- Sys.time()
  out <- ode(initial_state, times, Model_desolve, pars,
             Approx_I_sl=Approx_I_sl, Approx_I_ml=Approx_I_ml,
             Approx_temp=Approx_temp, Approx_moist=Approx_moist,
             spinup = spinup, end = end)
  print(Sys.time() - t0)
  return(out)
}
