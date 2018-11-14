rm(list = ls())
closeAllConnections()

library("lmom")
source('D:/Greg/Documents/R Scripts/l_moment_extensions.R')

lmrd()
ratios <- c(1, 0.1, 0.1, 0.4, 0)
# kap <- pelkap(ratios) #verify that kappa fails
wak_orig <- pelwak(ratios) #fit a mildly kurtotic wakeby instead (note hyperskewness is 0)

x <- quawak(runif(500), wak_orig)
x.mom <- samlmu(x, 5)
lmrdpoints(x.mom)

wak <- pelwak(x.mom)

R <- 50000
epy <- floor(seq(2, 100, length.out = 50))

start_time <- proc.time()
#Cluster setup
reserve_cores <- 1
cores_to_use <- min(length(epy), detectCores() - reserve_cores)

#Simulation
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("wak",
                  "quawak",
                  "samlmu",
                  "pelkap",
                  "R",
                  "epy",
                  "pelgev")
clusterExport(cl, cluster_vars)
wak_out <- foreach(i = 1:length(epy), .combine = "rbind") %dopar% {
  this_epy <- epy[i]
  R_out <- replicate(R, {
    storm_sample <- quawak(runif(this_epy), wak)
    max(storm_sample)
    
  })
  kap_par <- pelkap(samlmu(R_out))
  gev_par <- pelgev(samlmu(R_out))
  c(this_epy, gev_par["k"], kap_par["h"])
}
stopCluster(cl)

colnames(wak_out) <- c("n", "k", "h")
plot(wak_out[,"n"], wak_out[,"h"], 
     type = "l", lwd = 2, col = "darkgreen",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "Kappa h-Parameter")
abline(h = 0, lty = 3)

plot(wak_out[,"n"], wak_out[,"k"], 
     type = "l", lwd = 2, col = "darkgreen",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "GEV Kappa Parameter")
abline(h = 0, lty = 3)
runtime <- proc.time() - start_time
