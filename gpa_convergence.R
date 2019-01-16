rm(list = ls())
library("foreach")
library("doSNOW")
library("parallel")
library("lmom")
closeAllConnections()

CorrectedLMoments <- function(x, nmom = 5) {
  raw_moments <- samlmu(x, nmom = nmom)
  bias_t3 <- 4 * (length(x) ^ -1) * (0.10 - raw_moments["t_4"])
  bias_t4 <- 4 * (length(x) ^ -1) * (0.15 - raw_moments["t_4"])
  corrected_moments <- raw_moments
  corrected_moments["t_3"] <- corrected_moments["t_3"] - bias_t3
  corrected_moments["t_4"] <- corrected_moments["t_4"] - bias_t4
  return(corrected_moments)
}

colfunc <- colorRampPalette(c("red", "blue"))

S <- 1000000
lower_pct <- 0.1
upper_pct <- 0.99
seq_length <- 100


parent_moments <- lmrglo(c(1,0.1,-0.1), nmom = 5)
parent_dist <- pelglo(parent_moments)
parent_samp <- quaglo(runif(S), parent_dist)

thresholds_to_test <- seq(quantile(parent_samp, lower_pct), quantile(parent_samp, upper_pct), length.out = seq_length)
threshold_percentiles <- seq(lower_pct, upper_pct, length.out = seq_length)
color_list <- colfunc(length(thresholds_to_test))

plot(density(parent_samp), col = "black",
     main = "Parent Density",
     xlab = "Variate")
reserve_cores <- 1
cores_to_use <- min(length(thresholds_to_test), detectCores() - reserve_cores)
cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)
cluster_vars <- c("parent_dist",
                  "quaglo",
                  "samlmu",
                  "pelkap",
                  "S",
                  "thresholds_to_test",
                  "pelgev")
clusterExport(cl, cluster_vars)
gpa_out <- foreach(i = 1:length(thresholds_to_test), .combine = "rbind") %dopar% {
  this_thresh <- thresholds_to_test[i]
  big_storms <- parent_samp[parent_samp > this_thresh]
  big_storm_lmom <- CorrectedLMoments(big_storms, nmom = 5)
  big_storm_kappa <- pelkap(big_storm_lmom)
  c(this_thresh, big_storm_kappa["xi"], big_storm_kappa["alpha"], big_storm_kappa["k"], big_storm_kappa["h"], 
    big_storm_lmom["t_3"], big_storm_lmom["t_4"], big_storm_lmom["t_5"])
}
stopCluster(cl)
gpa_out <- rbind(c(0, NA, NA, NA, NA, parent_moments["tau_3"], parent_moments["tau_4"], parent_moments["tau_5"]), gpa_out)
gpa_out <- cbind(gpa_out, c(NA, threshold_percentiles))
colnames(gpa_out) <- c("threshold", "xi", "alpha", "k", "h", "t_3", "t_4", "t_5", "threshold_percentile")
rownames(gpa_out) <- 1:nrow(gpa_out)
gpa_out <- as.data.frame(gpa_out)

lmrd()
points(x = parent_moments["tau_3"], y = parent_moments["tau_4"])
last_t3 <- as.numeric(gpa_out[2, "t_3"])
last_t4 <- as.numeric(gpa_out[2, "t_4"])
colfunc <- colorRampPalette(c("red", "blue"))
color_list <- colfunc(nrow(gpa_out) - 3)
for(i in 3:nrow(gpa_out)) {
  if(i > 3) {
    last_t3 <- t3
    last_t4 <- t4
  }
  t3 <- as.numeric(gpa_out[i, "t_3"])
  t4 <- as.numeric(gpa_out[i, "t_4"])
  segments(x0 = last_t3, y0 = last_t4, x1 = t3, y1 = t4,
           col = color_list[i-1],
           lwd = 2)
}

with(gpa_out[2:nrow(gpa_out),],
     plot(x = threshold_percentiles,
          y = h,
          type = "l",
          col = "darkgreen",
          lwd = 3,
          xlab = "Threshold Percentile",
          ylab = "Kappa h Parameter",
          main = "Generalized Pareto Convergence",
          xlim = c(0.1, 1),
          ylim = c(-1, 1),
          log = "x"))
abline(h = 1, lty = 2)
points(x = 0.1, y = pelkap(parent_moments)["h"])
