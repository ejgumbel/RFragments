library(flexsurv)
library(readr)
library(stats4)
library(foreach)
library(iterators)
source('D:/Statistics/R Scripts/ev_plots.R')

#Functions
nll_gengamma <- function(x, mu, sigma, Q) {
  -sum(log(dgengamma(x, mu, sigma, Q)))
}

nll_gengamma.orig <- function(x, a, b, k) {
  -sum(log(dgengamma.orig(x, a, b, k)))
}


mle_gengamma <- function(x) {
  a <- mle(function(mu, sigma, Q) nll_gengamma(x, mu, sigma, Q), 
           start = list(mu = mean(x), sigma = sd(x), Q = 1),
           method = "L-BFGS-B", control = list(trace = 6),
           lower = c(0.01, 0.01, 0),
           upper = c(10, 10, 10))
  a@coef
  
}

mle_gengamma.orig <- function(x) {
  a <- mle(function(a, b, k) nll_gengamma.orig(x, a, b, k),
           start = list(a = 1, b = 1, k = 1),
           method = "Nelder-Mead", control = list(trace = 6))
  a@coef
}

mle_gengamma.orig_quiet <- function(x) {
  a <- mle(function(a, b, k) nll_gengamma.orig(x, a, b, k),
           start = list(a = 1, b = 1, k = 1),
           method = "Nelder-Mead")
  a@coef
}

boot_gengamma.orig <- function(x, S) {
  base_fit_param <- mle_gengamma.orig(x)
  foreach(i = 1:S, .combine = rbind) %do% {
    y <- rgengamma.orig(length(x), base_fit_param["a"], base_fit_param["b"], base_fit_param["k"])
    mle_gengamma.orig_quiet(y)
  }
}

#Data input
ba_typed_abvdallas <- read_csv("C:/Statistics/Generalized Gamma Experiments/ba_typed_abvdallas.csv", 
                               col_types = cols(Date = col_date(format = "%m/%d/%Y")))

#Data subsetting
meso <- ba_typed_abvdallas[ba_typed_abvdallas$`Mes/Syn` == "Meso",]
synop <- ba_typed_abvdallas[ba_typed_abvdallas$`Mes/Syn` == "Syn",]
meso_max <- aggregate(meso, list(Year = meso$Year), max)
synop_max <- aggregate(synop, list(Year = synop$Year), max)
all_max <- aggregate(ba_typed_abvdallas, list(Year = ba_typed_abvdallas$Year), max)

#Fitting
greg_evp(all_max$BA, max_T = 200, max_y = 1.2*max(all_max$BA))
fit <- mle_gengamma.orig(all_max$BA)
ggboot <- boot_gengamma.orig(all_max$BA, 100)


#b=iter(ggboot, by='row')
foreach(i = 1:nrow(ggboot), .combine = "c") %do% {
  b = ggboot[i,]
  lines(x = rgv(pp(1:1000, name = "Gringorten")),
        y = qgengamma.orig(pp(1:1000, name = "Gringorten"), b["a"], b["b"], b["k"]),
        col = "lightgrey", lty = 1, lwd = 0.5)
}

lines(x = rgv(pp(1:1000, name = "Gringorten")),
      y = qgengamma.orig(pp(1:1000, name = "Gringorten"), fit["a"], fit["b"], fit["k"]),
      lty = 1, lwd = 3)
points(x = rgv(pp(all_max$BA, name = "Gringorten")),
       y = sort(all_max$BA),
       pch = 18)
