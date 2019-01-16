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

S <- 1000
# S <- 1000
epy <- seq(min_epy, max_epy, by = 1)
color_list <- colfunc(length(epy))

gev <- c(1, 0.1, -0.1)
gev_mom <- lmrgev(gev, nmom = 5)
lmrd(gev_mom)
gev_samp <- quagev(runif(S), gev)
plot(density(gev_samp), col = "black",
     main = "AMS Density",
     xlab = "Variate")
reserve_cores <- 2
cores_to_use <- min(length(epy), detectCores() - reserve_cores)
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("gev",
                  "quagev",
                  "samlmu",
                  "pelkap",
                  "S",
                  "epy",
                  "pelgev")
clusterExport(cl, cluster_vars)
gev_out <- foreach(i = 1:length(epy), .combine = "rbind") %dopar% {
  this_epy <- epy[i]
  ams <- replicate(S, {
    storms <- quagev(runif(this_epy), gev)
    max(storms)
  })
  ams.mom <- samlmu(ams, nmom = 5)
  kap_par <- tryCatch({
    pelkap(ams.mom)
  }, error = function(e) {
    c("xi" = NA, "alpha" = NA, "k" = NA, "h" = NA)
  } 
  )
  gev_par <- tryCatch({
    pelgev(ams.mom)
  }, error = function(e) {
    c("xi" = NA, "alpha" = NA, "k" = NA)
  }
  )
  c(this_epy, gev_par["xi"], gev_par["alpha"], gev_par["k"], kap_par["h"], ams.mom["t_3"], ams.mom["t_4"], ams.mom["t_5"])
}
stopCluster(cl)
# legend("topright",
#        col = c("black", "red", "blue"),
#        legend = c("Parent", min_epy_text, max_epy_text),
#        lwd = 1,
#        lty = 1,
#        bty = "n")
gev_out <- rbind(c(0, NA, NA, NA, NA, gev_mom["tau_3"], gev_mom["tau_4"], gev_mom["tau_5"]), gev_out)
colnames(gev_out) <- c("n", "xi", "alpha", "k", "h", "t_3", "t_4", "t_5")
plot(gev_out[,"n"], gev_out[,"h"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "Kappa h-Parameter",
     main = "Kappa Convergence")
abline(h = 0, lty = 3)
plot(gev_out[,"n"], gev_out[,"k"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "GEV Kappa Parameter",
     main = "GEV Convergence")
abline(h = 0, lty = 3)

plot(gev_out[,"n"], gev_out[,"t_3"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Skew")
points(x = 0, y = gev_mom["tau_3"])
plot(gev_out[,"n"], gev_out[,"t_4"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Kurtosis")
points(x = 0, y = gev_mom["tau_4"])
plot(gev_out[,"n"], gev_out[,"t_5"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Cinco")
points(x = 0, y = gev_mom["tau_5"])


library(readr)
# write_delim(as.data.frame(gev_out), "D:/Statistics/Gumbel Investigations/gev_out_maxstable.csv", delim = ",", col_names = TRUE)


lmrd()
last_t3 <- as.numeric(gev_out[1, "t_3"])
last_t4 <- as.numeric(gev_out[1, "t_4"])
colfunc <- colorRampPalette(c("red", "blue"))
color_list <- colfunc(nrow(gev_out) - 2)
for(i in 2:nrow(gev_out)) {
  if(i > 2) {
    last_t3 <- t3
    last_t4 <- t4
  }
  t3 <- as.numeric(gev_out[i, "t_3"])
  t4 <- as.numeric(gev_out[i, "t_4"])
  segments(x0 = last_t3, y0 = last_t4, x1 = t3, y1 = t4,
           col = color_list[i-1],
           lwd = 2)
}

theoretical_scale_parameter <- function(n, gev_par) {
  m <- gev_par[1]
  s <- gev_par[2]
  alpha <- -1/gev_par[3]
  return(s * n ^ (1/alpha))
}

plot(x = theoretical_scale_parameter(epy, gev),
     y = gev_out[min_epy:max_epy, "alpha"],
     xlab = "Theoretical Scale Parameter",
     ylab = "Simulated Scale Parameter",
     main = "Stability Postulate Convergence",
     pch = 18)
abline(a = 0, b = 1)

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
