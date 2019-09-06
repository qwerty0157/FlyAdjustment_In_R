sabr <- function(strike, forward, tau, alpha, volvol, rho, beta){
  xi <- function(x){
    numerator <- sqrt(1 - 2 * rho * x + x^2) + x - rho
    denominator <- 1 - rho
    return(log(numerator / denominator))
  }
 
  z <- volvol / alpha * (forward * strike)^((1 - beta)/2) * log(forward / strike)
 
  ret1_denominator <- (forward * strike)^((1 - beta) / 2) *
    (1 + ((1 - beta)^2/24) * (log(forward / strike))^2 +
       ((1 - beta)^4/1920) * (log(forward / strike))^4)
  ret1 <- alpha / ret1_denominator
 
  ret2 <- z / xi(z)
 
  ret3 <- ((1 - beta)^2)/24 * alpha^2 / ((forward * strike)^(1 - beta)) +
    1/4 * (rho * beta * volvol* alpha) / ((forward * strike)^((1 - beta)/2)) +
    (2 - 3 * rho^2)/24 * volvol^2
  ret3 <- 1 + ret3 * tau
 
  return(ret1 * ret2 * ret3)
}
 
#-------------------------------------------------------------------
# SABR
#-------------------------------------------------------------------
itr <- 0
tol <- 1e-10
adj <- 1e-7
vol_25ss <- vol_25ms #initial guess for smile strangle
error <- 1
par_calibration <- c(-1, -1, 0, k_25c_ms, k_25p_ms) #initial parameters
 
while(error > tol){
 
  beta <- 1 #const
  fn <- function(args){
    alpha <- args[1]
    volvol <- args[2]
    rho <- args[3]
    k_25c <- args[4]
    k_25p <- args[5]
   
    alpha <- exp(alpha)
    volvol <- exp(volvol)
    rho <- tanh(rho)
   
    #residuals
    #eq 3.18 - ATM vol
    res1 <- (sabr(k_atm, forward, tau, alpha, volvol, rho, beta) - vol_atm)^2
   
    #eq 3.19 - risk reversal
    res2 <- (sabr(k_25c, forward, tau, alpha, volvol, rho, beta) -
               sabr(k_25p, forward, tau, alpha, volvol, rho, beta) - vol_25rr)^2
   
    #eq 3.20 - smile strangle
    res3 <- (0.5 * (sabr(k_25c, forward, tau, alpha, volvol, rho, beta) +
                      sabr(k_25p, forward, tau, alpha, volvol, rho, beta)) -
               sabr(k_atm, forward, tau, alpha, volvol, rho, beta) - vol_25ss)^2
   
    #eq 3.14 - 25delta
    res4 <- (getDelta(deltaType, -1, r_for, forward, k_25p, sabr(k_25p, forward, tau, alpha, volvol, rho, beta), tau) + 0.25)^2
    res5 <- (getDelta(deltaType, 1, r_for, forward, k_25c, sabr(k_25c, forward, tau, alpha, volvol, rho, beta), tau) - 0.25)^2
   
    return(res1 + res2 + res3 + res4 + res5)
  }
 
  par_calibration <- powell(par_calibration, fn)$par
 
 
  vol_25c_ms <- sabr(k_25c_ms, forward, tau, exp(par_calibration[1]), exp(par_calibration[2]), tanh(par_calibration[3]), beta)
  vol_25p_ms <- sabr(k_25p_ms, forward, tau, exp(par_calibration[1]), exp(par_calibration[2]), tanh(par_calibration[3]), beta)
  v_trial <- getOptionValue(1, spot, forward, r_for, r_dom, k_25c_ms, vol_25c_ms, tau) +
    getOptionValue(-1, spot, forward, r_for, r_dom, k_25p_ms, vol_25p_ms, tau)
 
  error <- v_trial - v_25_ms
 
  if(error > 0){
    vol_25ss <- vol_25ss - adj
  }else{
    vol_25ss <- vol_25ss + adj
  }
 
  cat(itr, error, par_calibration, "\n")
 
  itr <- itr + 1
}
 
#save parameters
par_sabr <- par_calibration
cat(exp(par_sabr[1]), exp(par_sabr[2]), tanh(par_sabr[3]))
 
#check risk reversal
(sabr(par_sabr[4], forward, tau, exp(par_sabr[1]), exp(par_sabr[2]), tanh(par_sabr[3]), 1) -
   sabr(par_sabr[5], forward, tau, exp(par_sabr[1]), exp(par_sabr[2]), tanh(par_sabr[3]), 1)) %>>%
{round(., 4)}
 
#check market strangle
v_25_ms
v_trial
 
#volatility smile using SABR
strike <- seq(from = 1, to = 2, by = 0.01)
mapply(sabr, strike, forward, tau, exp(par_sabr[1]), exp(par_sabr[2]), tanh(par_sabr[3]), 1) %>>%
{cbind(strike, .)} %>>%
{plot(., type="l", xlab="Strike", ylab="Volatility", main="Volatility smile using SABR for EURUSD")}