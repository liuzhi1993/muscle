% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MUSCLE-package.R
\name{MUSCLE}
\alias{MUSCLE}
\title{MUSCLE
Segmentation with MUSCLE.}
\usage{
MUSCLE(
  Y,
  q,
  beta = 0.5,
  dyadic = TRUE,
  split = FALSE,
  m = 2000,
  details = FALSE,
  deconv = FALSE,
  lag = 192
)
}
\arguments{
\item{Y}{observations.}

\item{q}{quantiles, can be simulated by simulQuantile_MUSCLE.}

\item{beta}{target quantile, beta = 0.5 for median.}

\item{dyadic}{boolian variable, indicates test statistics obtained from dyadic subintervals or on full subintervals.}

\item{split}{boolian variable, indicates whether the data will be split by subsets with size m.}

\item{m}{splitting size.}

\item{details}{boolian variable, indicates whether computation details will be displayed.}

\item{deconv}{boolian variable, indicates whether the deconvolution results will be computed.}

\item{lag}{shift lag (the size of dependency), only for deconvolution (deconv = TRUE).}
}
\value{
multiscale quantile segmentation results.
}
\description{
MUSCLE
Segmentation with MUSCLE.
}
\examples{
# example code
library(muscle)
# generate data
n = 2048
alpha = 0.1
beta = 0.5
df = 3
#changepoints of the block signal
blocks <- rep(c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03,
7.68, 15.37, 0), times =c(204, 62, 41, 164, 40, 308, 82, 430, 225, 41, 61,390))
blocks.cpt <- c(205, 267, 308, 472, 512, 820, 902,1332, 1557, 1598, 1659)

signal = blocks
signal.cpt = sort(c(blocks.cpt,390,666,1445))
sd = c(8,0.5,4,1)
set.seed(96)
noise = c(rt(390,df)*sd[1],rt(278,df)*sd[2],rt(779,df)*sd[3],rt(601,df)*sd[4])/sqrt(df/(df-2))
Y = signal + noise
## MUSCLE
# simulate quantiles
q_muscle = simulQuantile_MUSCLE(n,alpha = alpha,beta = beta)

# segmentation with MUSCLE
reg_muscle = MUSCLE(Y,q_muscle,beta)
plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35),xlab = "",
ylab = "", main = "MUSCLE")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle),type = "s",lwd = 2, col = "red")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )
## MUSCLE-S
# segmentation with MUSCLE-S with m = 500
reg_muscle_s = MUSCLE(Y,q_muscle,beta,split = TRUE, m = 500)

# segmentation with MUSCLE-S with m = 500, print computation details
reg_muscle_s = MUSCLE(Y,q_muscle,beta,split = TRUE, m = 500, details = TRUE)

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35),xlab = "",
ylab = "", main = "MUSCLE-S")
lines(1:n, signal, type = "s", lwd =2, col ="black")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )
lines(evalStepFun(reg_muscle_s),type = "s",lwd = 2, col = "blue")

}
