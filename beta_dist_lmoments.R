rm(list = ls())
library(lmom)
library(parallel)
library(doSNOW)
library(foreach)
qbeta_4p_affine <- function(p, a, c, alpha, beta) {
  qbeta_ <- qbeta(p, alpha, beta)
  return((c-a)*qbeta_ + a)
}

# test_samp <- qbeta_4p_affine(runif(10000), 0, 5, 0.15, 1)
# summary(test_samp)
# plot(density(test_samp))
# samp_lmom <- samlmu(test_samp)
# lmrd()
# lmrdpoints(samp_lmom)
# kappa_params <- pelkap(samp_lmom)
# 

closeAllConnections()

# library("lmom")
source('D:/Greg/Documents/R Scripts/l_moment_extensions.R')

R <- 100000
epy <- floor(seq(2, 100, length.out = 50))

start_time <- proc.time()
#Cluster setup
reserve_cores <- 1
cores_to_use <- min(length(epy), detectCores() - reserve_cores)

beta_mom <- samlmu(qbeta_4p_affine(runif(100000), 0, 5, 0.15, 1), nmom = 5)

#Simulation
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("qbeta_4p_affine",
                  "samlmu",
                  "pelkap",
                  "R",
                  "epy",
                  "pelgev")
clusterExport(cl, cluster_vars)
beta_out <- foreach(i = 1:length(epy), .combine = "rbind") %dopar% {
  this_epy <- epy[i]
  ams <- replicate(R, {
    storms <- qbeta_4p_affine(runif(this_epy), 0, 5, 0.15, 1)
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
  c(this_epy, gev_par["k"], kap_par["h"], ams.mom["t_3"], ams.mom["t_4"], ams.mom["t_5"])
}
stopCluster(cl)

beta_out <- rbind(c(0, NA, NA, beta_mom["t_3"], beta_mom["t_4"], beta_mom["t_5"]), beta_out)
colnames(beta_out) <- c("n", "k", "h", "t_3", "t_4", "t_5")
plot(beta_out[,"n"], beta_out[,"h"], 
     type = "l", lwd = 2, col = "darkgreen",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "Kappa h-Parameter")
abline(h = 0, lty = 3)

plot(beta_out[,"n"], beta_out[,"k"], 
     type = "l", lwd = 2, col = "darkgreen",
     ylim = c(-1, 1),
     xlab = "Number of Events Per Year",
     ylab = "GEV Kappa Parameter")
abline(h = 0, lty = 3)
runtime <- proc.time() - start_time

lmrd(xlim = c(-0.4, 0.6), ylim = c(-0.1, 0.3))
# by(beta_out, 1:nrow(beta_out), function(row) lmrdpoints(row[,"t_3"], row[,"t_4"]))
last_t3 <- as.numeric(beta_out[1, "t_3"])
last_t4 <- as.numeric(beta_out[1, "t_4"])
colfunc <- colorRampPalette(c("red", "blue"))
color_list <- colfunc(nrow(beta_out) - 2)
for(i in 2:nrow(beta_out)) {
  if(i > 2) {
    last_t3 <- t3
    last_t4 <- t4
  }
  t3 <- as.numeric(beta_out[i, "t_3"])
  t4 <- as.numeric(beta_out[i, "t_4"])
  segments(x0 = last_t3, y0 = last_t4, x1 = t3, y1 = t4,
           col = color_list[i-1],
           lwd = 2)
}
