### ----------------------------------- ###
###       Run Script                    ###
### ----------------------------------- ###

rm(list=ls())

################################################################################
### Settings
################################################################################

sitenum <- 1  ##  1=askov_a, 2=askov_b, 3=grignon, 4=kursk, 5=rothamsted, 6=ultuna, 7=versailles
# sitenum <- commandArgs(trailingOnly = TRUE)
opt_eq <- 0
opt_tr <- 1
flag_ads  <- 1  # simulate adsorption to minerals
flag_lea  <- 0  # simulate leakage
dec_fun   <- "MM" # One of: 'MM', '2nd', '1st'
pars.default.file <-  "parsets/pars_M_eq.csv"
pars.optim.file <- "parsets/pars_optim_2.csv"
pars.optimeq.file <- "parsets/pars_optimeq_2.csv"

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
source("CostTrans.R")
source("CostEquil.R")
source("ParsCalc.R")
source("StartRun.R")
source("MonthlyInput.R")
source('EquilFun.R')
source('util.R')

### === Parameter Options =====================
pars_default <- read.csv(pars.default.file, row.names = 1)
pars_default <- setNames(pars_default[[1]], row.names(pars_default))

pars_optimeq       <- read.csv(pars.optimeq.file, row.names = 1)
pars_optimeq_init  <- setNames(pars_optimeq[[1]], row.names(pars_optimeq))
pars_optimeq_lower <- setNames(pars_optimeq[[2]], row.names(pars_optimeq))
pars_optimeq_upper <- setNames(pars_optimeq[[3]], row.names(pars_optimeq))

pars_optim       <- read.csv(pars.optim.file, row.names = 1)
pars_optim_init  <- setNames(pars_optim[[1]], row.names(pars_optim))
pars_optim_lower <- setNames(pars_optim[[2]], row.names(pars_optim))
pars_optim_upper <- setNames(pars_optim[[3]], row.names(pars_optim))

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

site_data <- suppressMessages(read_csv("../input/site_data.csv"))
site_data <- site_data[site_data$sitenum == sitenum, ]
site <- site_data$site
climate <- site_data$climate
site_data <- setNames(suppressMessages(as.numeric(site_data)), colnames(site_data))
obs <- suppressMessages(read_csv("../input/obs.csv", skip = 1))
obs <- obs[obs$site == site,]
trans.input.file  <- paste("../input/input_trans_" , climate, ".csv", sep="")
input_trans <- suppressMessages(read_csv(trans.input.file, skip = 2))
input_spin <- input_trans[input_trans$year < 11, ]
input_spin <- MonthlyInput(input_spin)
input_trans <- MonthlyInput(input_trans)

C_obs <-obs$soc.t.ha[1]  # observationsin tons per hectare

if (opt_eq) {
  fit_eq <- optim(par = pars_optimeq_init, fn = CostEquil, pars_default = pars_default, C_obs = C_obs, site_data = site_data,
                  method = "L-BFGS-B", lower = pars_optimeq_lower, upper = pars_optimeq_upper)
  print(paste0('value = ',fit_eq$value))
  print(fit_eq$par)

  run_fiteq <- as.data.frame(StartRun(pars = pars_default, pars_new = fit_eq$par,
                                      site_data = site_data, input = input_trans, tsave = day))
  source('plots.R')
  out <- GetYearly(data = run_fiteq, obs = obs, steps = year/day)
  write_csv(out, path = paste0('S', sitenum, '_', site, '_eqrun.csv'))
}

if (opt_tr) {
  # fit_tr <- optim(fn =  CostTrans, par = pars_optim_init, pars_default = pars_default,
  #                 input_spin = input_spin, input_trans = input_trans,
  #                 method = "L-BFGS-B", upper = pars_optim_upper, lower = pars_optim_lower)
  fit_tr <- modFit(f = CostTrans, p = pars_optim_init, pars_default = pars_default, C_obs = C_obs, site_data = site_data,
                   input_spin = input_spin, input_trans = input_trans, method = "SANN",
                   upper = pars_optim_upper, lower = pars_optim_lower)
  run_fittr <- as.data.frame(StartRun(pars = pars_default, pars_new = fit_tr$par,
                                      site_data = site_data, input = input_trans))
}

save.image(file = paste0("Optim_", site, "_", site, "_", starttime, ".RData"))
