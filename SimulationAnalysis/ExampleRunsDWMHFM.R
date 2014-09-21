require(xtable)
require(ggplot2)

source("(3) MCMC Samplers/Gibbs/DW.R")
source("(3) MCMC Samplers/MH/VMMH.R")
source("(3) MCMC Samplers/Rejection/FM.R")
source('(5) Code for analysis/analysisHelperFunctions.R')
source("(6) Basic codes/DescribeCirc.R")


# Read in data data.
n <- 30
kappa <- .1
J <- 3
i <- 4
readfilename <- paste0(wd, "/(2) Data/Datasets", 
                       "J=", J, "_n=", n, "_kap=", kappa, "/nr", i, ".csv")
th <- read.table(readfilename, sep=",")

Q <- 500

set.seed(1997)
dwres <- DW(th, Q=Q)
mhres <- VMMH(th, Q=Q)
fmres <- FM(th, Q=Q)


scale <- 1


pdf("Examplerun.pdf", width=12*scale, height=8*scale)
par(mfrow=c(2, 3))
par(cex=1, cex.axis=1, cex.lab=1.9, cex.main=1.3, oma=c(0,3,0,0))

ylimmu <- range(setCircularRangeContinuous(mhres$mu[, 1])*(180/pi))
ylimkp <- c(0, 1)

plot(setCircularRangeContinuous(dwres$mu[, 1])*(180/pi), type="l", main="Gibbs"
     , ylab="mu"
     , xlab="", ylim=ylimmu)
plot(setCircularRangeContinuous(mhres$mu[, 1])*(180/pi), type="l", main="MH"
     , ylab=""
     , xlab="", ylim=ylimmu)
plot(setCircularRangeContinuous(fmres$mu[, 1])*(180/pi), type="l", main="Rejection"
     , ylab=""
     , xlab="", ylim=ylimmu)

plot((dwres$kappa), type="l", #main="Gibbs", 
     ylab="kappa"
     , xlab=""
     , ylim=ylimkp
     )
plot((mhres$kappa), type="l", #main="MH",
     ylab=""
     , xlab=""
     , ylim=ylimkp
     )
plot((fmres$kappa), type="l", #main="Rejection",
     ylab=""
     , xlab=""
     , ylim=ylimkp
     )

dev.off()




