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

epy <- 100
S <- 5000

l1 <- 2.6
l2 <- 0.1
t3 <- seq(0, 0.5, by = 0.01)
t4 <- seq(0, 0.3, by = 0.01)

lmom::lmrd()

for(t3_i in seq_along(t3)) {
  for(t4_i in seq_along(t4)) {
    tryCatch({
      kap <- lmom::pelkap(c(l1, l2, t3[t3_i], t4[t4_i]))
      k <- est(kap, epy, S)
      if(k < -0.05){
        lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "red", pch = 18)
      } else {
        if(k >  0.05){
          lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "blue", pch = 18)
        } else {
          lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "black", pch = 18)
        }
      }
      
    }, error = function(e) {
      lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "grey")
    }, warning = function(w) {
      lmom::lmrdpoints(t3[t3_i], t4[t4_i], col = "grey", pch = 18)
    })
    
  }
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
