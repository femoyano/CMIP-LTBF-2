# Plots
library(readr)

# sitenum <- 7
# sdata <- read_csv("site_data.csv")
# sdata <- sdata[sdata$sitenum == sitenum,]
# obs <- read_csv('obs.csv', skip = 1)
# obs <- obs[obs$sitenum == sitenum,]
# site <- sdata$site
# mod <- read_csv(paste0('output_1/output_trans_', site, '.csv'))

mod <- run_fiteq
SOC <- apply(mod[,-1], MARGIN = 1, sum) * 10000 / 1000  # sum and conversion to tons C per ha
# SOC <- (mod$C_P + mod$C_A) * 10000 / 1000 # sum and conversion to tons C per ha
mod$year <- mod$time/360/24
plot(SOC~mod$year, ylim=c(0,max(c(SOC, obs$soc.t.ha))), main = site, type = 'l')
points(obs$soc.t.ha~obs$time, col = 3, pch = 16)
