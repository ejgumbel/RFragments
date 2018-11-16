rm(list = ls())
library("foreach")
library("lmom")
closeAllConnections()

colfunc <- colorRampPalette(c("red", "blue"))

min_epy <- 2
min_epy_text = paste0(min_epy, " Events Per Year")
max_epy <- 100
max_epy_text = paste0(max_epy, " Events Per Year")

S <- 50000
epy <- seq(min_epy, max_epy, by = 1)
color_list <- colfunc(length(epy))

ratios <- c(l_1 = 2.6, l_2 = 0.1, t_3 = 0.25, t_4 = 0.35, t_5 = 0.03)
wak <- pelwak(ratios, verbose = TRUE)
wak_mom <- lmrwak(wak)
lmrd(lmrwak(wak))
wak_samp <- quawak(runif(S), wak)
# plot(density(wak_samp), col = "black",
#      main = "AMS Density",
#      xlab = "Variate")
reserve_cores <- 2
cores_to_use <- min(length(epy), detectCores() - reserve_cores)
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("wak",
                  "quawak",
                  "samlmu",
                  "pelkap",
                  "S",
                  "epy",
                  "pelgev")
clusterExport(cl, cluster_vars)
wak_out <- foreach(i = 1:length(epy), .combine = "rbind") %dopar% {
  this_epy <- epy[i]
  ams <- replicate(S, {
    storms <- quawak(runif(this_epy), wak)
    max(storms)
  })
  ams.mom <- samlmu(ams)
  kap_par <- pelkap(ams.mom)
  gev_par <- pelgev(ams.mom)
  # lines(density(ams), col = color_list[i])
  c(this_epy, gev_par["k"], kap_par["h"], ams.mom["t_3"], ams.mom["t_4"])
}
stopCluster(cl)
# legend("topright",
#        col = c("black", "red", "blue"),
#        legend = c("Parent", min_epy_text, max_epy_text),
#        lwd = 1,
#        lty = 1,
#        bty = "n")
colnames(wak_out) <- c("n", "k", "h", "t_3", "t_4")
plot(wak_out[,"n"], wak_out[,"h"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "Kappa h-Parameter",
     main = "Kappa Convergence")
abline(h = 0, lty = 3)
plot(wak_out[,"n"], wak_out[,"k"], 
     type = "l", lwd = 2, col = "darkblue",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "GEV Kappa Parameter",
     main = "GEV Convergence")
abline(h = 0, lty = 3)

plot(wak_out[,"n"], wak_out[,"t_3"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Skew")

plot(wak_out[,"n"], wak_out[,"t_4"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Kurtosis")