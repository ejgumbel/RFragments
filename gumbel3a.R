rm(list = ls())
library("foreach")
library("doSNOW")
library("parallel")
library("lmom")
closeAllConnections()

colfunc <- colorRampPalette(c("red", "blue"))

min_epy <- 1
min_epy_text = paste0(min_epy, " Events Per Year")
max_epy <- 100
max_epy_text = paste0(max_epy, " Events Per Year")

S <- 1000000
# S <- 1000
epy <- seq(min_epy, max_epy, by = 1)
color_list <- colfunc(length(epy))

ratios <- c(l_1 = 2.6, l_2 = 0.1, t_3 = 0.45, t_4 = 0.13, t_5 = 0.02573354)
kap <- pelkap(ratios)
kap_mom <- lmrkap(kap, nmom = 5)
lmrd(lmrkap(kap))
kap_samp <- quakap(runif(S), kap)
plot(density(kap_samp), col = "black",
     main = "AMS Density",
     xlab = "Variate")
reserve_cores <- 2
cores_to_use <- min(length(epy), detectCores() - reserve_cores)
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("kap",
                  "quakap",
                  "samlmu",
                  "pelkap",
                  "S",
                  "epy",
                  "pelgev")
clusterExport(cl, cluster_vars)
kap_out <- foreach(i = 1:length(epy), .combine = "rbind") %dopar% {
  this_epy <- epy[i]
  ams <- replicate(S, {
    storms <- quakap(runif(this_epy), kap)
    max(storms)
  })
  ams.mom <- samlmu(ams, nmom = 5)
  # kap_par <- pelkap(ams.mom)
  kap_par <- tryCatch({
    pelkap(ams.mom)
  }, error = function(e) {
    c("xi" = NA, "alpha" = NA, "k" = NA, "h" = NA)
  } 
  )
  # gev_par <- pelgev(ams.mom)
  gev_par <- tryCatch({
    pelgev(ams.mom)
  }, error = function(e) {
    c("xi" = NA, "alpha" = NA, "k" = NA)
  }
  )
  # lines(density(ams), col = color_list[i])
  c(this_epy, gev_par["k"], kap_par["h"], ams.mom["t_3"], ams.mom["t_4"], ams.mom["t_5"])
}
stopCluster(cl)
# legend("topright",
#        col = c("black", "red", "blue"),
#        legend = c("Parent", min_epy_text, max_epy_text),
#        lwd = 1,
#        lty = 1,
#        bty = "n")
kap_out <- rbind(c(0, NA, NA, kap_mom["tau_3"], kap_mom["tau_4"], kap_mom["tau_5"]), kap_out)
colnames(kap_out) <- c("n", "k", "h", "t_3", "t_4", "t_5")
plot(kap_out[,"n"], kap_out[,"h"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "Kappa h-Parameter",
     main = "Kappa Convergence")
abline(h = 0, lty = 3)
plot(kap_out[,"n"], kap_out[,"k"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "GEV Kappa Parameter",
     main = "GEV Convergence")
abline(h = 0, lty = 3)

plot(kap_out[,"n"], kap_out[,"t_3"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Skew")
points(x = 0, y = kap_mom["tau_3"])
plot(kap_out[,"n"], kap_out[,"t_4"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Kurtosis")
points(x = 0, y = kap_mom["tau_4"])
plot(kap_out[,"n"], kap_out[,"t_5"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Cinco")
points(x = 0, y = kap_mom["tau_5"])


library(readr)
write_delim(as.data.frame(kap_out), "D:/Statistics/Gumbel Investigations/kap_out3a.csv", delim = ",", col_names = TRUE)


lmrd(xlim = c(-0.5, 0.5), ylim = c(0, 0.2))
# by(wak_out, 1:nrow(wak_out), function(row) lmrdpoints(row[,"t_3"], row[,"t_4"]))
last_t3 <- as.numeric(kap_out[1, "t_3"])
last_t4 <- as.numeric(kap_out[1, "t_4"])
colfunc <- colorRampPalette(c("red", "blue"))
color_list <- colfunc(nrow(kap_out) - 2)
for(i in 2:nrow(kap_out)) {
  if(i > 2) {
    last_t3 <- t3
    last_t4 <- t4
  }
  t3 <- as.numeric(kap_out[i, "t_3"])
  t4 <- as.numeric(kap_out[i, "t_4"])
  segments(x0 = last_t3, y0 = last_t4, x1 = t3, y1 = t4,
           col = color_list[i-1],
           lwd = 2)
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
