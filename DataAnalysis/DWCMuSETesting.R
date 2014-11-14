library(dplyr)

# Load EVERYTHING
source('Simulation/SimulationStudyVM.R')


set.seed(1001)
th <- list(rvmc(10, 2.3, 2), rvmc(20, 2, 2))

Q <- 100000
rdw <- DW(th, Q = Q, lag = 30)
rmh <- VMMH(th, Q = Q)
rfm <- FM(th, Q = Q)


par(mfrow=c(1, 3))
rdw$mu[, 1] %>% density() %>% plot(xlim=c(1.5, 3), ylim=c(0, 3))
rdw$mu[, 2] %>% density() %>% lines()

rmh$mu[, 1] %>% density() %>% plot(xlim=c(1.5, 3), ylim=c(0, 3))
rmh$mu[, 2] %>% density() %>% lines()

rfm$mu[, 1] %>% density() %>% plot(xlim=c(1.5, 3), ylim=c(0, 3))
rfm$mu[, 2] %>% density() %>% lines()




