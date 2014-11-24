# ----------------------------------------------------------
# runSimulationStudy.R
# Script to run simulationStudyVM.R for several different scenario's,
# specifically with different sampling methods, and different number of groups,
# saving the results to .rda R-object files.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the Dutch
# Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------


source('Simulation/SimulationStudyVM.R')
source('Data/GenerateData.R')

nsim <- 1000

# Generate the data.
set.seed(8921747)
Q <- 10000
ns  <- c(5, 30, 100)
kps <- c(0.1, 4, 32)
generateVonMisesData(nsim, ns=ns, kaps=kps, J=c(1, 3))


# The rejection and MH samplers are straightforward to simulate as they do not
# require seperate values for lag for each cell.

# J=1

set.seed(105843)
simulationStudyVM(samplername="FM",   Q=Q, ns=ns, kappas=kps, nsim=nsim, J=1)

set.seed(190362)
simulationStudyVM(samplername="VMMH", Q=Q, ns=ns, kappas=kps, nsim=nsim, J=1)

# J=3
set.seed(215428)
simulationStudyVM(samplername="FM",   Q=Q, ns=ns, kappas=kps, nsim=nsim, J=3, printsim=TRUE)

set.seed(353397)
simulationStudyVM(samplername="VMMH", Q=Q, ns=ns, kappas=kps, nsim=nsim, J=3)



# The Gibbs sampler is given a separate lag value that was tuned to the specific
# datasets.

# J=1
set.seed(367628)
simulationStudyVM(samplername="DW", nsim=nsim, ns=5,   kappas=0.1, Q=Q, J=1, lag=2, printsim=TRUE)
set.seed(941321)
simulationStudyVM(samplername="DW", nsim=nsim, ns=5,   kappas=4,   Q=Q, J=1, lag=25, printsim=TRUE)

set.seed(758945)
simulationStudyVM(samplername="DW", nsim=nsim, ns=30,  kappas=0.1, Q=Q, J=1, lag=3, printsim=TRUE)
set.seed(711397)
simulationStudyVM(samplername="DW", nsim=nsim, ns=30,  kappas=4,   Q=Q, J=1, lag=30, printsim=TRUE)

set.seed(222638)
simulationStudyVM(samplername="DW", nsim=nsim, ns=100, kappas=0.1, Q=Q, J=1, lag=10, printsim=TRUE)
set.seed(117440)
simulationStudyVM(samplername="DW", nsim=nsim, ns=100, kappas=4,   Q=Q, J=1, lag=150, printsim=TRUE)


# J=3
set.seed(891246)
simulationStudyVM(samplername="DW", nsim=nsim, ns=5,   kappas=0.1, Q=Q, J=3, lag=9, printsim=TRUE)
set.seed(419283)
simulationStudyVM(samplername="DW", nsim=nsim, ns=5,   kappas=4,   Q=Q, J=3, lag=60, printsim=TRUE)

set.seed(845006)
simulationStudyVM(samplername="DW", nsim=nsim, ns=30,  kappas=0.1, Q=Q, J=3, lag=10, printsim=TRUE)
set.seed(895087)
simulationStudyVM(samplername="DW", nsim=nsim, ns=30,  kappas=4,   Q=Q, J=3, lag=65, printsim=TRUE)

set.seed(621705)
simulationStudyVM(samplername="DW", nsim=nsim, ns=100, kappas=0.1, Q=Q, J=3, lag=35, printsim=TRUE)
set.seed(800436)
simulationStudyVM(samplername="DW", nsim=nsim, ns=100, kappas=4,   Q=Q, J=3, lag=250, printsim=TRUE)



