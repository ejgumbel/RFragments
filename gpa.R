dgpa <- function(x, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- (x - xi) / alpha
  } else {
    y <- -k^(-1) * log(1 - k * (x - xi)/ alpha)
  }
  alpha^(-1) * exp(-(1-k)*y)
}

pgpa <- function(q, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- (q - xi) / alpha
  } else {
    y <- -k^(-1) * log(1 - k * (q - xi)/ alpha)
  }
  1 - exp(-y)
}

qgev <- function(p, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- xi - alpha * log(1 - p)
  } else {
    y <- xi + alpha * (1-(1-p)^k)/k
  }
  y
}

rgpa <- function(n, xi = 0, alpha = 1, k = 0) {
  qgpa(runif(n), xi = xi, alpha = alpha, k = k)
}


nllgpa <- function(x, xi, alpha, kappa) {
  R <- suppressWarnings(dgpa(x, xi, alpha, kappa))
  if(any(is.nan(R))) {
    Inf
  } else {
    -sum(log(R))
  }
}

mle_gpa <- function(x) {
  library("stats4")
  par.init <- lmom::pelgpa(lmom::samlmu(x))
  a <- mle(function(xi, alpha, kappa) nllgpa(x, xi, alpha, kappa),
           start = list(xi = par.init[[1]], alpha = par.init[[2]], kappa = par.init[[3]]),
           method = "Nelder-Mead")
  a@coef
}

boot_gpa <- function(x, S) {
  library("foreach")
  base_fit_param <- mle_gpa(x)
  foreach(i = 1:S, .combine = rbind) %do% {
    y <- rgpa(length(x), base_fit_param["xi"], base_fit_param["alpha"], base_fit_param["kappa"])
    mle_gpa(y)
  }
}