require(tidyverse)


pt <-read.csv("data/cleaned/plot_transformations.csv")
tt <- read.csv("data/cleaned/trait_transformations.csv")

backtransform <- function(x, var) {
  
  out <- NULL
  
  # P50 log transform
  if(var == "P50_log") {
    out <- -exp(x*tt$scale[1] + tt$center[1])
  
  # P50 square root transform
  } else if(var == "P50_sqrt") {
    out <- -(x*tt$scale[2] + tt$center[2])^2
    
  # rooting depth log transform
  } else if(var == "rd_log") {
    out <- exp(x*tt$scale[3] + tt$center[3])
    
  # aridity log transform  
  } else if(var == "arid_log") {
    out <- exp(x*pt$scale[1] + pt$center[1])
    
  # aridity square root transform  
  } else if(var == "arid_sqrt") {
    out <- (x*pt$scale[2] + pt$center[2])^2
    
  # WTD log transform
  } else if(var == "wtd_log") {
    out <- exp(x*pt$scale[3] + pt$center[3]) - 0.01


# precip seasonality log transform  
  } else if(var == "ps_log") {
    out <- exp(x*pt$scale[4] + pt$center[4])
  }
  return(out)
}


transform <- function(x, var) {
  
  out <- NULL
  
  # P50 log transform
  if(var == "P50_log") {
    out <- (log(abs(x)) - tt$center[1])/tt$scale[1]
    
  # P50 square root transform
  } else if(var == "P50_sqrt") {
    out <- (sqrt(abs(x)) - tt$center[2])/tt$scale[2]
    
  # Rooting depth log transform
  } else if(var == "rd_log") {
    out <- (log(x) - tt$center[3])/tt$scale[3]
    
  # aridity log transform
  } else if(var == "arid_log") {
    out <- (log(x) - pt$center[1])/pt$scale[1]
    
  # aridity square root transform
  } else if(var == "arid_sqrt") {
    out <- (sqrt(x) - pt$center[2])/pt$scale[2]
    
  # WTD log transform  
  } else if(var == "wtd_log") {
    out <- (log(x+0.01) - pt$center[3])/pt$scale[3]
    
  # precip seasonality log transform
  } else if(var == "ps_log") {
    out <- (log(x) - pt$center[4])/pt$scale[4]
  }
  
  return(out)
}
