# MonthlyInput.R

MonthlyInput <- function(x) {

  out <- x %>%
    group_by(year, month) %>%
    select(year, month, hour, temp, moist, litt) %>%
    summarise(
      hour = mean(hour, na.rm = TRUE),
      temp = mean(temp, na.rm = TRUE),
      moist = mean(moist, na.rm = TRUE),
      litt = mean(litt, na.rm = TRUE)
    )
return(out)
}
