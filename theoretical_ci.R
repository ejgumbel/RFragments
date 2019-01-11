rm(list = ls())
library(lmom)

est <- function(kap_par, epy = 500, S = 10000) {
  ams <- replicate(S, {
    storms <- lmom::quakap(runif(epy), kap_par)
    max(storms)
  })
  kap <- lmom::pelkap(lmom::samlmu(ams))
  return(as.numeric(kap["k"]))
}


test_moments <- c(1.0, 0.3, 0.2, 0.17)
test_kap_parameters <- pelkap(test_moments)
conv_k <- sign(est(test_kap_parameters))
test_n <- 50

plot_z <- seq(qnorm(1-0.99), qnorm(1-1/10^4), by = 0.1)
plot_p <- pnorm(plot_z)
oms <- 1:4 #creates a few powers of ten to plot
vls <- c(1-10^-oms[oms<3], 0.5, 10^-oms) #places to add "axis" lines
xmax <- qnorm(1-10^-4) #pad the plot to the next order of magnitude
xmin <- qnorm(10^-min(2, 4)) #Mirror
ylabel <- paste0("Kappa Distribution h = ", round(test_kap_parameters["h"], 1))

plot(x = plot_z,
     y = quakap(plot_p, test_kap_parameters),
     type = "l",
     col = "black",
     lty = 1,
     lwd = 2,
     xlab = "Exceedance Probability",
     ylab = ylabel,
     xlim = c(xmin, xmax),
     ylim = c(0, 5),
     axes = F)

axis(2)
axis(1, at = qnorm(1 - vls), labels = vls)
abline(v = qnorm(1-vls), lty = 2, col = "gray")
box()

boots <- 500
ct <- 0
p_test <- numeric()
exc_test_Val <- 2.5


while(ct < boots) {
  x <- quakap(runif(test_n), test_kap_parameters)
  tryCatch ({
    cand <- pelkap(samlmu(x))
    if(sign(cand["h"]) != conv_k) {
    }
    else {
      ct <- ct + 1
      lines(x = plot_z,
            y = quakap(plot_p, cand),
            col = "grey")
      p_test <- c(p_test, cdfkap(exc_test_Val, cand))
    }
  }, error = function(e) {

  })

}

abline(h = exc_test_Val, lty = 3, col = "black")
abline(v = qnorm(mean(p_test)), lty = 3, col = "black")
lines(x = plot_z,
     y = quakap(plot_p, test_kap_parameters),
     type = "l",
     col = "black",
     lty = 1,
     lwd = 2)
