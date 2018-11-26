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
epy <- seq(min_epy, max_epy, by = 1)
color_list <- colfunc(length(epy))

ratios <- c(l_1 = 2.6, l_2 = 0.1, t_3 = 0.2, t_4 = 0.4, t_5 = 0.2)
wak <- pelwak(ratios, verbose = TRUE)
wak_mom <- lmrwak(wak)
lmrd(lmrwak(wak))
wak_samp <- quawak(runif(S), wak)
plot(density(wak_samp), col = "black",
     main = "AMS Density",
     xlab = "Variate")
reserve_cores <- 4
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
wak_out <- rbind(c(0, NA, NA, wak_mom["tau_3"], wak_mom["tau_4"], wak_mom["tau_5"]), wak_out)
colnames(wak_out) <- c("n", "k", "h", "t_3", "t_4", "t_5")
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
     ylim = c(0, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Skew")
points(x = 0, y = wak_mom["tau_3"])
plot(wak_out[,"n"], wak_out[,"t_4"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Kurtosis")
points(x = 0, y = wak_mom["tau_4"])
plot(wak_out[,"n"], wak_out[,"t_5"], 
     type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1),
     xlab = "Number of Events Per Year",
     ylab = "L-Cinco")
points(x = 0, y = wak_mom["tau_5"])


library(readr)
write_delim(as.data.frame(wak_out), "D:/Statistics/Gumbel Investigations/wak_out.csv", delim = ",", col_names = TRUE)
lmrd(xlim = c(0, 0.75), ylim = c(0, 0.5))
# by(wak_out, 1:nrow(wak_out), function(row) lmrdpoints(row[,"t_3"], row[,"t_4"]))
last_t3 <- wak_out[1, "t_3"]
last_t4 <- wak_out[1, "t_4"]
colfunc <- colorRampPalette(c("red", "blue"))
color_list <- colfunc(nrow(wak_out) - 2)
for(i in 2:nrow(wak_out)) {
  if(i > 2) {
    last_t3 <- t3
    last_t4 <- t4
  }
  t3 <- wak_out[i, "t_3"]
  t4 <- wak_out[i, "t_4"]
  segments(x0 = last_t3, y0 = last_t4, x1 = t3, y1 = t4,
           col = color_list[i-1],
           lwd = 2)
}