rm(list = ls())

library(foreach)
library(lmom)

rgv <- function(p) {
  y <- suppressWarnings(-log(-log(p)))
  return(y)
}

##### Parameters #####
gev_par <- c(xi = 3.1496,
             alpha = 0.5045,
             kappa = -0.1627)
N <- 50
R <- 10000
alpha <- 0.1
set.seed(12345)

##### Where are results desired? #####
aep <- c(0.95, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, 0.000005, 0.000002, 0.000001, 0.0000005)
cdf <- 1 - aep
pop_quants <- quagev(cdf, gev_par)

##### Bootstrap Procedure #####
result <- foreach(1:R, .combine = rbind) %do% {
  x <- quagev(runif(N), gev_par)
  x.lmom <- samlmu(x)
  x.gevpar <- pelgev(x.lmom)
  cdfgev(pop_quants, x.gevpar)
}
  
##### Computation of Expected Probability and Percentiles #####
exp_prob <- colMeans(result)
upper_conf <- apply(result, 2, function(x) quantile(x, alpha/2))
lower_conf <- apply(result, 2, function(x) quantile(x, 1-alpha/2))

##### Plotting Logic for EV-1 Plot #####
max_T <- 1/min(aep)
if(!is.na(max_T)) {
  om <- ceiling(log10(max_T))
} else {
  om <- floor(log10(n))+1 #finds the power of ten to plot to
}
oms <- 1:om #creates a few powers of ten to plot
vls <- c(1-10^-oms[oms<3], 3, 0.5, 10^-oms) #places to add "axis" lines
xmax <- rgv(1-10^-om) #pad the plot to the next order of magnitude
xmin <- rgv(10^-min(2, om)) #Mirror

plot(x = rgv(cdf), y = pop_quants,
     axes = F,
     xlab = "", ylab = "",
     xlim = c(xmin, xmax),
     ylim = c(0, max(pop_quants)),
     type = "l", col = "black", lwd = 2)
     
axis(2)
axis(1, at = rgv(1 - vls), labels = vls)
abline(v = rgv(1-vls), lty = 3, col = "gray") # don't mess with this; sets the AEP dotted lines
abline(h = seq(0, 30, 5), lty = 3, col = "gray") # ok to change; sets the quantile dotted lines
box()
title(xlab = "Exceedance Probability", line = 2)
title(ylab = "GEV Quantiles", line = 2.25)
title(main = "Bootstrap Validation")

##### Adding Results to Plot ######
lines(x = rgv(exp_prob), y = pop_quants, col = "red", lwd = 2)
lines(x = rgv(lower_conf), y = pop_quants, col = "red", lwd = 1, lty = 3)
lines(x = rgv(upper_conf), y = pop_quants, col = "red", lwd = 1, lty = 3)

##### Legend #####
legend("topleft", 
       legend = c("Upper Confidence",
                  "Expected Probability",
                  "Computed Frequency Curve",
                  "Lower Confidence"),
       lty = c(3, 1, 1, 3),
       col = c("red", "red", "black", "red"),
       lwd = c(1, 2, 2, 1))

##### Create Plaintext Output #####
output_df <- data.frame("Base_AEP" = aep,
                        "Computed_Curve" = pop_quants,
                        "Expected_Curve" = exp_prob,
                        "Lower_Conf" = lower_conf,
                        "Upper_Conf" = upper_conf)
write.csv(output_df,
          "D:/Statistics/Lmoment_bootstrap_valid.csv",
          row.names = FALSE)

