### ----------------------------------- ###
###       Run Script                    ###
### ----------------------------------- ###

rm(list=ls())

t0 <- Sys.time()
starttime  <- format(t0, "%m%d-%H%M")
cmdarg <- commandArgs(trailingOnly = TRUE)
source(paste0(cmdarg, "settings.R"))

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

### === Parameter Options =====================
pars <- read.csv(pars.default.file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))

pars_optim       <- read.csv(pars.optim.file, row.names = 1)
pars_optim_init  <- setNames(pars_optim[[1]], row.names(pars_optim))
pars_optim_lower <- setNames(pars_optim[[2]], row.names(pars_optim))
pars_optim_upper <- setNames(pars_optim[[3]], row.names(pars_optim))

### Variables =================================================================
year  <- 31104000      # seconds in a year
month <- 2592000
week  <- 604800
hour  <- 3600          # seconds in an hour
sec   <- 1             # seconds in a second!

### === Other Options =====================
tstep  <- hour

################################################################################
### Prepare site data
################################################################################

obs <- read_csv("../input/obs.csv", skip = 1)
obs <- obs[obs$sitenum == sitenum,]
site <- obs$site[1]
site_data <- read_csv("../input/site_data.csv")
site_data <- setNames(as.numeric(site_data[site_data$sitenum == sitenum,-1]), colnames(site_data[,-1]))
trans.input.file  <- paste("../input/input_trans_" , site, ".csv", sep="")
input_trans <- read_csv(trans.input.file, skip = 2)
input_spin <- input_trans[input_trans$year < 11, ]
input_spin$litt <- 0.00005 # in kgC m-2 h-1
input_spin <- MonthlyInput(input_spin)
input_trans <- MonthlyInput(input_trans)

if (opt_eq) {
  C_obs <-obs$soc.t.ha[1]  # observationsin tons per hectare
  fit_eq <- optim(par = pars_optim_init, fn = CostEquil, pars = pars, C_obs = C_obs,
                  method = "L-BFGS-B", lower = pars_optim_lower, upper = pars_optim_upper)
  pars_optim_init <-fit_eq$par
}

run_fiteq <- StartRun(pars = pars, pars_new = fit_eq$par, site_data = site_data,
                        input = input_trans)

if (opt_tr) {
  # fit_tr <- optim(fn =  CostTrans, par = pars_optim_init, pars = pars,
  #                 input_spin = input_spin, input_trans = input_trans,
  #                 method = "L-BFGS-B", upper = pars_optim_upper, lower = pars_optim_lower)
  fit_tr <- modFit(f = CostTrans, p = pars_optim_init, pars = pars, C_obs = C_obs,
                   input_spin = input_spin, input_trans = input_trans, method = "SANN",
                   upper = pars_optim_upper, lower = pars_optim_lower)
}

run_fittr <- StartRun(pars = pars, pars_new = fit_tr$par, site_data = site_data,
                        input = input_trans)


save.image(file = paste0("Optim_", cmdarg, "_", site, "_", starttime, ".RData"))
