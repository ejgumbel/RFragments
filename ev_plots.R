#rm(list = ls())
pp <- function(x, name = "Weibull", a, transform = F) {
  if(!missing(a)) {
    a <- a
  } else {
    a <- switch(name,
                Weibull = {0},
                Median = {0.3175},
                Blom = {0.375},
                Cunnane = {0.40},
                Gringorten = {0.44},
                Hazen = {0.5})
  }
  n <- length(x)
  qi <- (1:n - a)/(n + 1 - 2*a)
  if(name == "Median") {
    qi[1] <- 1-0.5^(1/n)
    qi[n] <- 1-qi[1]
  }
  if(transform == "norm") {
    qi <- qnorm(qi)
  } else if(transform == "gumbel") {
    qi <- rgv(qi)
  }
  return(qi)
}

rgv <- function(p) {
  y <- suppressWarnings(-log(-log(p)))
  return(y)
}

rev_rgv <- function(y) exp(-exp(-y))

qgumbel <- function(p, xi = 0, alpha = 1) {
  xi + alpha * rgv(p)
}

rgumbel <- function(n, xi = 0, alpha = 1) {
  p <- runif(n)
  qgumbel(p, xi = xi, alpha = alpha)
}

pgumbel <- function(q, xi = 0, alpha = 1) {
  exp(-exp(-(q - xi)/alpha))
}

gumbel_mom <- function(x) {
  em <- -digamma(1)
  s <- sd(x)
  u <- mean(x)
  a <- sqrt(6) * s / pi
  x <- u - em * a
  out <- c(xi = x, alpha = a)
  return(out)
}

#extreme value plots using rgv
#measure of tail thickness using k derived from l-moments

greg_evp <- function(x, max_T = NA, max_y = 1, plot_title = "", y_title = NA, color = "black") {
  #Greg's extreme value plot for univariate data
  x <- sort(x) #Sort the data set immediately
  pp <- pp(x, name = "Gringorten") #Comupute plotting positions (Gringorten since we are using RGV)
  y <- rgv(pp) #Convert to RGVs, the standard variable for this is "y" which will be fun later
  n <- length(x)
  if(!is.na(max_T)) {
    om <- ceiling(log10(max_T))
  } else {
    om <- floor(log10(n))+1 #finds the power of ten to plot to
  }
  oms <- 1:om #creates a few powers of ten to plot
  vls <- c(1-10^-oms[oms<3], 3, 0.5, 10^-oms) #places to add "axis" lines
  xmax <- rgv(1-10^-om) #pad the plot to the next order of magnitude
  xmin <- rgv(10^-min(2, om)) #Mirror
  #ymin <- median(x)
  #ymax <- max(x)
  plot(x = y, y = x, 
       axes = F,
       xlab = "", ylab = "",
       xlim = c(xmin, xmax),
       ylim = c(0, max_y),
       pch = 18,
       col = color)
  axis(2)
  axis(1, at = rgv(1 - vls), labels = vls)
  abline(v = rgv(1-vls), lty = 2, col = "gray")
  box()
  title(xlab = "Exceedance Probability", line = 2)
  if(is.na(y_title)) {
    yt = "Quantile"
  } else {
    yt = y_title
  }
  title(ylab = yt, line = 2.25)
  title(main = plot_title)
}

gumbel_rejector <- function(x, alpha = 0.05) {
  suppressWarnings(library("lmom"))
  t <- samlmu(x)
  gp <- pelgev(t)
  n <- length(x)
  Z <- gp["k"]*sqrt(n/0.5633)
  Zc <- qnorm(1-alpha/2)
  detach("package:lmom", unload=TRUE)
  if(abs(Z) > Zc) return(T) else return(F)
}

ramp_em <- function(x) {
  rbPal <- colorRampPalette(c("darkgreen", "darkorchid"))
  rbPal(length(x))[as.numeric(cut(x[[2]], breaks = length(x)))]
}