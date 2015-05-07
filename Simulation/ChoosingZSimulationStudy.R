source('Simulation/DW_chosen_k_sampler.R')
sourceCpp('Data/rvmc.cpp')

require(xtable)

# Z is the maximum number of N_k that was selected by the algorithm: if we never
# choose a value above some value Z, it is not necessary to sample values above
# that point.
getMaxZ <- function (kp, n, J) {
  th <- lapply(1:J, function(x) rvmc(n, 2, kp))
  res <- DW(th, Q=10000, burn=500, Z=40)$spec$chosenk
  max(res)
}

repeatGetMaxZ <- function (n, kp, J, nrep) {
  cat("\n\n\nCell: [[", paste(n, kp, J, nrep), "]], datasets:\n", sep="")
  sapply(1:nrep, function(x) {
    cat(x, ", ", sep = "")
    getMaxZ(kp=kp, n=n, J=J)
    }
  )
}

runSimStudyK <- function (ns, kps, J, nrep, seed) {
  set.seed(seed)

  # 3-d array to collect output.
  out <- array(dim=c(nrep, length(ns), length(kps)), dimnames=list(1:nrep, ns, kps))

  # Conditions
  cds <- as.data.frame(t(expand.grid(n=ns, kp=kps)))

  # Run the simulation study.
  for (cd in cds) {
    out[ , as.character(cd[1]), as.character(cd[2])] <-
      repeatGetMaxZ(n=cd[1], kp=cd[2], J=J, nrep=nrep)
  }
  out
}

# Set simulation conditions.
nrep <- 100
kps  <- c(0.1, 1, 4, 8, 16, 32)
ns   <- c(10, 30, 100)

resultj1 <- runSimStudyK(ns=ns, kps=kps, J=1, nrep=nrep, seed=25)
save.image(file="choosingZSimStudyResults.RData")
resultj3 <- runSimStudyK(ns=ns, kps=kps, J=3, nrep=nrep, seed=27)
save.image(file="choosingZSimStudyResults.RData")

DWtabj1 <- t(apply(resultj1, 2:3, max))
tabj3 <- t(apply(resultj3, 2:3, max))

combres <- cbind(tabj1, rep("", 4), tabj3)


# Output to LaTeX.
print(xtable(combres, digits=0, align=c("c@{\\hskip 0.7cm}", "c", "c", "c", "c@{\\hskip 0.7cm}", "c", "c", "c"),
       caption=paste0("Maximum $k$ that was picked out as the smallest value ",
                      "after 10000 iterations of the Gibbs sampler applied to 100 datasets for ",
                      "different sample sizes ($n$), concentration ($\\kappa$) and number of groups $(J)$."),
       caption.placement="top",
       label="tab:sufficientk"),
       include.rownames=TRUE, #Don't print rownames
       include.colnames=FALSE, #We create them ourselves
       caption.placement="top",
       hline.after=NULL, #We don't need hline; we use booktabs
       add.to.row=list(pos     = as.list(c(-1, -1, -1, 0, nrow(combres))),
                       command = c("\\toprule \n  & & $J=1$ & & & & $J=3$ & \\\\ \n",
                                   " & \\multicolumn{3}{c}{Sample size} & & \\multicolumn{3}{c}{Sample size} \\\\ \n",
                                   " $\\kappa$ & 10 & 30 & 100 & & 10 & 30 & 100 \\\\ \n",
                                   "\\hline \n ",
                                   "\\bottomrule \n"))
       )
