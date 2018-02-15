### ----------------------------------- ###
###       Run Script (All)              ###
### ----------------------------------- ###

rm(list=ls())

################################################################################
### Settings
################################################################################


opt_all <- 1
flag_ads  <- 1  # simulate adsorption to minerals
flag_lea  <- 0  # simulate leakage
dec_fun   <- "MM" # One of: 'MM', '2nd', '1st'
pars.default.file <-  "parsets/pars_M_eq.csv"

pars.optim.all.file <- "parsets/pars_optim_all.csv"

spin.years   <- 2    # years for spinup run
flag.cmi     <- 0     # use a constant mean input for spinup

t0 <- Sys.time()
starttime  <- format(t0, "%m%d-%H%M")

### Libraries =================================================================
library(deSolve)
library(reshape2)
library(FME)
library(readr)
library(dplyr)

### Sourced required files ====================================================
source("Model.R")
source("flux_functions.R")
source("GetEquil.R")
source("cost_all.R")
source("ParsCalc.R")
source("StartRun.R")
source("MonthlyInput.R")
source('EquilFun.R')
source('util.R')
source('Plots.R')

### === Parameter Options =====================
pars_default <- read.csv(pars.default.file, row.names = 1)
pars_default <- setNames(pars_default[[1]], row.names(pars_default))

pars_optim_all           <- read.csv(pars.optim.all.file, row.names = 1)
pars_optim_all_init      <- setNames(pars_optim_all[[1]],row.names(pars_optim_all))
pars_optim_all_lower     <- setNames(pars_optim_all[[2]],row.names(pars_optim_all))
pars_optim_all_upper     <- setNames(pars_optim_all[[3]],row.names(pars_optim_all))



### Variables =================================================================
year  <- 31104000      # seconds in a year
month <- 2592000
week  <- 604800
day   <- 86400
hour  <- 3600          # seconds in an hour
sec   <- 1             # seconds in a second!

### === Other Options =====================
tstep  <- hour

################################################################################
### Prepare site data
################################################################################

site_data <- suppressMessages(read_csv("./input/site_data.csv"))
site_data <- transform(site_data, eq_litt = c(0.0000085, 0.0000089, 0.000007, 0.0000061, 0.000016, 0.0000057, 0.000012))


obs <- suppressMessages(read_csv("./input/obs.csv", skip = 1))

trans.input.file  <- "./input/input_all_weather.csv"
input_trans <- suppressMessages(read_csv(trans.input.file, skip = 0))
input_spin <- input_trans
#input_spin <- MonthlyInput(input_spin)
#input_trans <- MonthlyInput(input_trans)

C_obs <-obs$soc.t.ha[1]  # observationsin tons per hectare


if (opt_all) {
  # fit_tr <- optim(fn =  CostTrans, par = pars_optim_init, pars_default = pars_default,
  #                 input_spin = input_spin, input_trans = input_trans,
  #                 method = "L-BFGS-B", upper = pars_optim_upper, lower = pars_optim_lower)
  fit_tr <- modFit(f = CostAll, p = pars_optim_all_init, pars_default = pars_default, C_obs = C_obs, site_data = site_data,
                   input_spin = input_spin, input_trans = input_trans, method = "Nelder-Mead",
                   upper = pars_optim_all_upper, lower = pars_optim_all_lower)
#  run_fittr <- as.data.frame(StartRun(pars = pars_default, pars_new = fit_tr$par,
#                                      site_data = site_data[5,], input = input_trans[input_trans$climate == site_data[5,1],]))
#  Plot1(run_fittr, obs[obs$site == site_data$site[5],])
#  out <- GetYearly(data = run_fittr, obs = obs, steps = year/month)
  #write_csv(out, path = paste0('S', sitenum, '_', site, '_trrun.csv'))
}

for (i in 1:7){
  site <- site_data[i,1]
  pars_default["eq_litt"] <- site_data[i,13]
input_trans_ind <- MonthlyInput(input_trans[input_trans$climate == site_data[i,2],])
  run_fittr <- as.data.frame(StartRun(pars = pars_default, pars_new = fit_tr$par,
                                    site_data = site_data[i,], input = input_trans_ind))
Plot1(run_fittr, obs[obs$site == site_data$site[i],])

}

#out <- GetYearly(data = run_fittr, obs = obs, steps = year/month)



# save.image(file = paste0("Optim_", site, "_", starttime, ".RData"))
