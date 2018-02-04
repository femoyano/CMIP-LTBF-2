### Calculate equilibrium state from a given set of parameters and site data.

sitenum <- 7 ##  1=askov_a, 2=askov_b, 3=grignon, 4=kursk, 5=rothamsted, 6=ultuna, 7=versailles
pars.default.file <- 'parsets/pars_M2H_test.csv' #"parsets/pars_M2H_test.csv"
pars.new.file <- "parsets/pars_new_0.csv"
flag_ads  <- 1  # simulate adsorption to minerals
flag_lea  <- 0  # simulate leakage
dec_fun   <- "MM" # One of: 'MM', '2nd', '1st'

library(readr)
source("GetEquil.R")

### === Parameter Options =====================
pars <- read.csv(pars.default.file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))
pars_new  <- read.csv(pars.new.file, row.names = 1)
pars_new  <- setNames(pars_new[[1]], row.names(pars_new))
site_data <- suppressMessages(read_csv("../input/site_data.csv"))
site_data <- site_data[site_data$sitenum == sitenum, ]
site <- site_data$site
climate <- site_data$climate
site_data <- suppressMessages(setNames(as.numeric(site_data), colnames(site_data)))
# for(n in names(pars_new)) pars[[n]] <- pars_new[[n]]
# site_data[['eq_temp']] <- 276.4545
# pars[['eq_litt']] <- 0.00005
print(GetEquil(pars = pars, site_data = site_data))
