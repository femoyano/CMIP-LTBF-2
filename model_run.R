# run script

rm(list=ls())

sitenum <- 1 ##  1=askov_a, 2=askov_b, 3=grignon, 4=kursk, 5=rothamsted, 6=ultuna, 7=versailles
pars.default.file <- "parsets/pars_M2H_test.csv"
pars.new.file <- "parsets/pars_new_0.csv"
flag_ads  <- 1  # simulate adsorption to minerals
flag_lea  <- 0  # simulate leakage
diff_fun  <- "hama"  # Options: 'hama', 'cubic'
dec_fun   <- "MM" # One of: 'MM', '2nd', '1st'
upt_fun   <- "1st" # One of: 'MM', '2nd', '1st'

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
source("ParsCalc.R")
source("StartRun.R")
source("MonthlyInput.R")

### === Parameter Options =====================
pars <- read.csv(pars.default.file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))
pars_new  <- read.csv(pars.new.file, row.names = 1)
pars_new  <- setNames(pars_new[[1]], row.names(pars_new))

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

site_data <- suppressMessages(read_csv("../input/site_data.csv"))
site_data <- site_data[site_data$sitenum == sitenum, ]
site <- site_data$site
climate <- site_data$climate
site_data <- suppressMessages(setNames(as.numeric(site_data), colnames(site_data)))
obs <- suppressMessages(read_csv("../input/obs.csv", skip = 1))
obs <- obs[obs$site == site,]
trans.input.file  <- paste("../input/input_trans_" , climate, ".csv", sep="")
input_trans <- suppressMessages(read_csv(trans.input.file, skip = 2))
input_spin <- input_trans[input_trans$year < 11, ]
input_spin$litt <- 0.00005 # in kgC m-2 h-1
input_spin <- MonthlyInput(input_spin)
input_trans <- MonthlyInput(input_trans)

initial_state <- if (site == "askovb3")     c(C_P = 4.541, C_D = 0.040032, C_A = 0.580102, C_M = 0.095935)   else
  if (site == "askovb4")     c(C_P = 4.159, C_D = 0.043,    C_A = 0.58102, C_M = 0.095935) else
    if (site == "grignon")     c(C_P = 2.669, C_D = 0.048032, C_A = 1.400102, C_M = 0.045935) else
      if (site == "kursk")       c(C_P = 8.419, C_D = 0.098032, C_A = 1.380102, C_M = 0.105935) else
        if (site == "rothamsted")  c(C_P = 5.166, C_D = 0.048032, C_A = 0.900102, C_M = 0.095935) else
          if (site == "ultuna")      c(C_P = 2.606, C_D = 0.048032, C_A = 1.500102, C_M = 0.095935) else
            if (site == "versailles")  c(C_P = 5.539, C_D = 0.035958, C_A = 0.980102, C_M = 0.114297)

mod <- as.data.frame(StartRun(pars = pars, pars_new = pars_new, site_data = site_data,
         input = input_trans, initial_state = initial_state, tsave = year))

SOC <- apply(mod[,-1], MARGIN = 1, sum)
SOC <- SOC * 10000 / 1000  # conversion to tons C per ha
year <- mod$time/360/24
plot(SOC~year, ylim=c(0,max(SOC)), main = site, type = 'l')
points(obs$soc.t.ha~obs$time, col = 2)
