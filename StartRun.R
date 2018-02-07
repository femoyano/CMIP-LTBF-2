# Run_spinup.R

StartRun <- function(input, site_data,  pars_default = pars_default, pars_new = NULL, initial_state = "equil", spin.years = 10, spinup = 0, tsave = month) {

  pars_default <- c( pars_default, site_data)
  for(n in names(pars_new)) pars_default[[n]] <- pars_new[[n]]
  pars <- ParsCalc(pars_default)

  if (initial_state == "equil") {
    initial_state <- GetEquil(pars = pars)
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
