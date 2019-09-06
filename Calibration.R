library(pipeR)
library(powell)
 
 
#-------------------------------------------------------------------
# Initial settings for EURUSD 1y
#-------------------------------------------------------------------
r_for <- 0.03446
r_dom <- 0.0294
spot <- 1.3465
tau <- 1
pips <- 4
forward <- (spot * exp((r_dom - r_for) * tau)) %>>% {round(., pips)}
deltaType <- "pipspot"
 
vol_atm <- 0.1825
vol_25ms <- 0.0095 #market strangle
vol_25rr <- -0.006 #risk reversal
 
#-------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------
getD <- function(forward, strike, vol, tau, flg){
  #flg = 1 or -1
  return((log(forward/strike) + flg * 0.5 * vol^2 * tau)/(vol * sqrt(tau)))
}
 
getDelta <- function(deltaType, callput, r_for, forward, strike, vol, tau){
  #callput = 1 for call, -1 for put
  if(tolower(deltaType) == "pipspot"){
    d1 <- getD(forward, strike, vol, tau, 1)
    return(callput * exp(-r_for * tau) * pnorm(callput * d1))
  }
}
 
getOptionValue <- function(callput, spot, forward, r_for, r_dom, strike, vol, tau){
  #callput = 1 for call, -1 for put
  d1 <- getD(forward, strike, vol, tau, 1)
  d2 <- getD(forward, strike, vol, tau, -1)
  return(
    callput * spot * exp(-r_for * tau) * pnorm(callput * d1) -
      callput * strike * exp(-r_dom * tau) * pnorm(callput * d2)
  )
}
 
getATMStrike <- function(deltaType, forward, vol, tau){
  if(tolower(deltaType) == "pipspot"){
    return(forward * exp(0.5 * vol^2 * tau))
  }
}
 
#-------------------------------------------------------------------
# Calibration
#-------------------------------------------------------------------
k_atm <- getATMStrike(deltaType, forward, vol_atm, tau) %>>% {round(., pips)}
 
k_25c_ms <- uniroot(function(x){
  return(
    getDelta(deltaType, 1, r_for, forward, x, vol_atm + vol_25ms, tau) - 0.25
  )
}, c(0, 100))$root %>>% {round(., pips)}
 
k_25p_ms <- uniroot(function(x){
  return(
    getDelta(deltaType, -1, r_for, forward, x, vol_atm + vol_25ms, tau) + 0.25
  )
}, c(0, 100))$root %>>% {round(., pips)}
 
v_25_ms <- getOptionValue(-1, spot, forward, r_for, r_dom, k_25p_ms, vol_atm + vol_25ms, tau) +
  getOptionValue(1, spot, forward, r_for, r_dom, k_25c_ms, vol_atm + vol_25ms, tau)