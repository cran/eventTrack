exactDatesFromMonths <- function(predicted, nevent){
     
     # from monthly event numbers, compute exact date of targeted number of events
     # via linear interpolation
     
     date.nevent <- approx(x = predicted[, 2], y = predicted[, 1], xout = nevent)$y     
     return(date.nevent)
}





