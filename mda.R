#for each pair of t3 and t4, fit a kappa distribution
#then find the maximum of a large number of events
#repeat this a lot
#compute the kappa-h

rm(list = ls())

est <- function(kap_par, epy = 500, S = 10000) {
  ams <- replicate(S, {
    storms <- lmom::quakap(runif(epy), kap_par)
    max(storms)
  })
  kap <- lmom::pelkap(lmom::samlmu(ams))
  return(as.numeric(kap["k"]))
}

plot_4kappa_line <- function(h) {
  test_k <- seq(-0.95, 0.95, 0.001)
  test_t3 <- test_k
  test_t4 <- test_k
  for(i in 1:length(test_k)){
    test_t3[i] <- lmrkap(para = c(0, 1, test_k[i], h))["tau_3"]
    test_t4[i] <- lmrkap(para = c(0, 1, test_k[i], h))["tau_4"]
  }
  lmrdlines(test_t3[complete.cases(test_t3)], test_t4[complete.cases(test_t3)],
            lwd = 4, col = "red")
}

plot_4kappa_line_by_k <- function(k) {
  test_h <- seq(-1.0, 10, 0.1)
  test_t3 <- test_h
  test_t4 <- test_h
  for(i in 1:length(test_h)){
    test_t3[i] <- lmom::lmrkap(para = c(0, 1, k, test_h[i]))["tau_3"]
    test_t4[i] <- lmom::lmrkap(para = c(0, 1, k, test_h[i]))["tau_4"]
  }
  lmom::lmrdlines(test_t3[complete.cases(test_t3)], test_t4[complete.cases(test_t3)],
                  lwd = 2, col = "black", lty = 2)
}

epy <- 100
S <- 5000

l1 <- 2.6
l2 <- 0.1
t3 <- seq(0, 0.5, by = 0.01)
t4 <- seq(0, 0.3, by = 0.01)

lmom::lmrd(xlim = c(0, 0.5), ylim = c(0, 0.3))

for(t3_i in seq_along(t3)) {
  for(t4_i in seq_along(t4)) {
    tryCatch({
      kap <- lmom::pelkap(c(l1, l2, t3[t3_i], t4[t4_i]))
      k <- est(kap, epy, S)
      if(k < 0){
        lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "red", pch = 18)
      } else {
        lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "blue", pch = 18)
      }
      
    }, error = function(e) {
      # lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "grey")
    }, warning = function(w) {
      # lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "grey", pch = 18)
    })
    
  }
}

plot_4kappa_line_by_k(0)
legend("bottomright", 
       legend = c("Fréchet", "Weibull", "Gumbel"), 
       lty = c(NA, NA, 2), 
       pch = c(18, 18, NA), 
       col = c("red", "blue", "black"), 
       lwd = c(NA, NA, 2), 
       cex = 1.0, 
       title = "Maximum Domain of Attraction")
