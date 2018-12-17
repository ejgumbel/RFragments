dgev <- function(x, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- (x - xi) / alpha
  } else {
    y <- -k^(-1) * log(1 - k * (x - xi)/ alpha)
  }
  alpha^(-1) * exp(-(1-k)*y-exp(-y))
}

pgev <- function(q, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- (q - xi) / alpha
  } else {
    y <- -k^(-1) * log(1 - k * (q - xi)/ alpha)
  }
  exp(-exp(-y))
}

qgev <- function(p, xi = 0, alpha = 1, k = 0) {
  if(abs(k) < 0.001) {
    y <- xi - alpha * log(-log(p))
  } else {
    y <- xi + alpha * (1-(-log(p))^k)/k
  }
  y
}

rgev <- function(n, xi = 0, alpha = 1, k = 0) {
  qgev(runif(n), xi = xi, alpha = alpha, k = k)
}


nllgev <- function(x, xi, alpha, kappa) {
  R <- suppressWarnings(dgev(x, xi, alpha, kappa))
  if(any(is.nan(R))) {
    Inf
  } else {
    -sum(log(R))
  }
}

nllgev_lc <- function(x, xi, alpha, kappa, nl, Cl) {
  #Left censoring of nl observations with censoring threshold Cl
  nCl <- rep(Cl, nl)
  exact <- suppressWarnings(sum(log(dgev(x, xi, alpha, kappa))))
  lc <- suppressWarnings(sum(log(pgev(nCl, xi, alpha, kappa))))
  if(any(is.nan(exact)) | any(is.nan(lc))) {
    Inf
  } else {
    -(exact + lc)
  }
}


mle_gev <- function(x) {
  library("stats4")
  par.init <- lmom::pelgev(lmom::samlmu(x))
  a <- mle(function(xi, alpha, kappa) nllgev(x, xi, alpha, kappa),
           start = list(xi = par.init[[1]], alpha = par.init[[2]], kappa = par.init[[3]]),
           method = "BFGS")
  a@coef
  
}

boot_gev <- function(x, S) {
  library("foreach")
  base_fit_param <- mle_gev(x)
  foreach(i = 1:S, .combine = rbind) %do% {
    y <- rgev(length(x), base_fit_param["xi"], base_fit_param["alpha"], base_fit_param["kappa"])
    mle_gev(y)
  }
}

# mlegev_lc <- function(x, nl, Cl) {
#   #Left censoring of nl observations with censoring threshold Cl
#   
#   par.init <- lmom::pelgev(lmom::samlmu(x))
#   a <- stats4::mle(function(xi, alpha, kappa) nllgev_lc(x, xi, alpha, kappa, nl, Cl), 
#            start = list(xi = par.init[[1]], alpha = par.init[[2]], kappa = par.init[[3]]),
#            method = "L-BFGS-B",
#            lower = c(-Inf, 0, -0.25),
#            upper = c(Inf, Inf, 1))
#   a@coef
# }

mlegev_lc <- function(x, nl, Cl) {
  #Left censoring of nl observations with censoring threshold Cl
  
  par.init <- lmom::pelgev(lmom::samlmu(x))
  a <- stats4::mle(function(xi, alpha, kappa) nllgev_lc(x, xi, alpha, kappa, nl, Cl), 
                   start = list(xi = par.init[[1]], alpha = par.init[[2]], kappa = par.init[[3]]),
                   method = "BFGS")
  a@coef
}

evplot_gevline <- function(x, mle_par, col = "black", lwd = 1, lty = 1) {
  lines(x = rgv(pp(x, "Gringorten")), y = lmom::quagev(pp(x, "Gringorten"), mle_par),
        col = col, lwd = lwd, lty = lty)
}

moments_gev <- function(xi = 0, alpha = 1, k = 0) {
  #eventually should return all the moments
  g <- function(r, k) {
    k <- -1*k
    gamma(1-r*k)
  }
  skew <- function(k) {
    k <- -1*k
    sign(k) * ((g(3,k)-3*g(2,k)*g(1,k)+2*g(1,k)^3)/((g(2,k)-(g(1,k))^2)^(3/2)))
  }
  variance <- function(k, a) {
    a^2*(g(2,k)-g(1,k)^2)/(k^2)
  }
  return(c(variance(k,alpha), -skew(k)))
  
}

# test_x <- seq(0, 10, 0.01)
# plot(x = test_x, y = dgev(test_x, 3, 1, -0.2), type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "f(x)", main = "GEV Densities")
# lines(x = test_x, y = dgev(test_x, 3, 1, 0), col = "black", lwd = 2)
# lines(x = test_x, y = dgev(test_x, 3, 1, 0.2), col = "green", lwd = 2)
# legend("topright", legend = c("Fréchet k = -0.2", "Gumbel k = 0", "Weibull k = 0.2"), 
#        col = c("blue", "black", "green"), lwd = c(2, 2, 2), lty = c(1,1,1), bty = "n")
