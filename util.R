# Useful

GetYearly <- function(data, obs, steps) {
  data <- as.data.frame(data)
  data[,-1] <- round(data[,-1] * 10, 2) # conversion to tons C per ha
  data$C <- apply(data[,-1], 1, sum) # out[['C_P']] + out[['C_A']] + out[['C_D']] + out[['C_M']]
  data <- data[0:(nrow(data)-1)%%steps==0, ] # get values every 12 months, i.e. one per year
  data$year <- c(0:(nrow(data)-1))
  # data <- data[data$year %in% obs$time, ]
  return(data)
}
