################################################################################
### Settings
################################################################################

opt_eq <- 1
opt_tr <- 0
flag_ads  <- 1  # simulate adsorption to minerals
flag_lea  <- 0  # simulate leakage
diff_fun  <- "hama"  # Options: 'hama', 'cubic'
dec_fun   <- "MM" # One of: 'MM', '2nd', '1st'
upt_fun   <- "1st" # One of: 'MM', '2nd', '1st'
pars.default.file <-  "parsets/pars5test.csv"
pars.optim.file <- "parsets/pars_optim_1.csv"
pars.optim.file <- "parsets/pars_optimeq_1.csv"
spin.years   <- 1    # years for spinup run
flag.cmi     <- 0     # use a constant mean input for spinup
