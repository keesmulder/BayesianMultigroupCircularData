# ----------------------------------------------------------
# GenerateTables.R
# Create tables to be input into LaTeX format.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the Dutch
# Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------

require(xtable)
require(Rcpp)

source('SimulationAnalysis/analysisHelperFunctions.R')
source("DataAnalysis/DescribeCirc.R")
sourceCpp("Data/rvmc.cpp")

# LOAD IN DATA
oldwd <- getwd()

setwd("Simulation/Results/")

# Combine DW data
# Old
# DWj1k.1 <- combineRes(c("[SimResultVM_v6]__[nsim1000]__[DW]__[n5]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2014-04-20 10.13.18].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n30]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2014-04-20 12.40.49].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n100]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2014-04-20 16.33.55].rda"), onDim=3)
# DWj1k4  <- combineRes(c("[SimResultVM_v6]__[nsim1000]__[DW]__[n5]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2014-04-20 12.20.54].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n30]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2014-04-20 15.27.03].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n100]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2014-04-17 15.46.11].rda"), onDim=3)
# DWj3k.1 <- combineRes(c("[SimResultVM_v6]__[nsim1000]__[DW]__[n5]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2014-04-20 10.22.35].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n30]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2014-04-20 16.20.13].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n100]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2014-04-20 23.02.41].rda"), onDim=3)
# DWj3k4  <- combineRes(c("[SimResultVM_v6]__[nsim1000]__[DW]__[n5]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2014-04-20 15.59.14].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n30]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2014-04-20 21.58.57].rda",
#                         "[SimResultVM_v6]__[nsim1000]__[DW]__[n100]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2014-04-18 02.54.39].rda"), onDim=3)


DWj1k.1 <- combineRes(c("[SimResultVM_v6]__[nsim2000]__[DW]__[n10]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2015-02-12 17.52.51].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n30]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2015-02-12 23.52.18].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n100]__[k0.1]__[mudif0.35]__[J1]__[Q10000]__[2015-02-13 08.18.43].rda"), onDim=3)

DWj1k4  <- combineRes(c("[SimResultVM_v6]__[nsim2000]__[DW]__[n10]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2015-02-12 23.14.43].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n30]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2015-02-13 06.17.10].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n100]__[k4]__[mudif0.35]__[J1]__[Q10000]__[2015-02-14 16.20.29].rda"), onDim=3)


DWj3k.1 <- combineRes(c("[SimResultVM_v6]__[nsim2000]__[DW]__[n10]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2015-03-04 13.45.49].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n30]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2015-04-07 22.20.44].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n100]__[k0.1]__[mudif0.35]__[J3]__[Q10000]__[2015-04-08 11.24.59].rda"), onDim=3)

DWj3k4  <- combineRes(c("[SimResultVM_v6]__[nsim2000]__[DW]__[n10]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2015-03-13 13.35.26].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n30]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2015-04-08 06.46.31].rda",
                        "[SimResultVM_v6]__[nsim2000]__[DW]__[n100]__[k4]__[mudif0.35]__[J3]__[Q10000]__[2015-04-08 14.54.44].rda"), onDim=3)

# Make an empty result set to put in instead of k=32
DWNAj1       <- DWj1k4
DWNAj1[,,,,] <- NA
dimnames(DWNAj1)[[4]] <- "32"
DWNAj3       <- DWj3k4
DWNAj3[,,,,] <- NA
dimnames(DWNAj3)[[4]] <- "32"

# Combine these.
DWj1 <- abind(DWj1k.1, DWj1k4, DWNAj1, along=4)
DWj3 <- abind(DWj3k.1, DWj3k4, DWNAj3, along=4)


# Combine DW data with MH and FM data.
rlj1 <- list(DW=DWj1,
             MH=loadssres("[SimResultVM_v6]__[nsim2000]__[VMMH]__[n10,30,100]__[k0.1,4,32]__[mudif0.35]__[J1]__[Q10000]__[2015-02-12 14.17.02].rda"),
             FM=loadssres("[SimResultVM_v6]__[nsim2000]__[FM]__[n10,30,100]__[k0.1,4,32]__[mudif0.35]__[J1]__[Q10000]__[2015-02-12 13.54.02].rda"))
rlj3 <- list(DW=DWj3,
             MH=loadssres("[SimResultVM_v6]__[nsim2000]__[VMMH]__[n10,30,100]__[k0.1,4,32]__[mudif0.35]__[J3]__[Q10000]__[2015-02-12 17.26.15].rda"),
             FM=loadssres("[SimResultVM_v6]__[nsim2000]__[FM]__[n10,30,100]__[k0.1,4,32]__[mudif0.35]__[J3]__[Q10000]__[2015-02-14 22.57.45].rda"))

# Save dimension names and number of methods
dn <- dimnames(rlj1[[1]])
nmeth <- length(rlj1)

str(rlj3)

mean(rlj3$MH[, 39, 1, 1, 1])
mean(rlj3$MH[, 39, 1, 2, 1])
rlj3$MH[, 39, 1, 1, 1]

# CREATE TABLES
# Create comparison tables.
j1res <- comparisontable(rl=rlj1, fnames=c("DW", "MH", "FM"), ks=c("0.1", "4", "32"))
j1resclean <- cleanRes(j1res)

# j3 needs averaged coverages for the mean
j3resbig <- comparisontable(rl=rlj3, J=3, fnames=c("DW", "MH", "FM"), ks=c("0.1", "4", "32"))
j3res <- cbind(j3resbig[, 1:6], round(rowMeans(apply(j3resbig[, 7:9], 2, as.numeric)), 2), j3resbig[, 10:13])

j3resclean <- cleanRes(j3res)



# nsplits denotes the rows where we get a new value for the sample sizes.
# Before, we want a line. The first index already has a line above (thus [-1]), and all
# values will be one too high as we want the line before the actual row (thus -1).
nsplits <- which(j1resclean[, 1]!="")[-1] - 1
spaceLocsFull <- seq(from=nmeth, to=nrow(j3resclean)-nmeth, by=nmeth)-1
spaceLocs <- spaceLocsFull[spaceLocsFull!=nsplits-1]

# Footnotes
kapfoot <- "\\footnote{Posterior $\\kappa$ mode denotes the mode as described in section \\ref{hdimode}. Coverage denotes the proportion of replications in which the true $\\kappa$ fell within the 95 \\% HDI.}"
accfoot <- "\\footnote{Acceptance ratio. For Gibbs sampling, this is always 1. For Metropolis-Hastings, $Q_{acc}/Q$ is given. For the rejection method, this is $Q/Q_{can}$.}"
mctfoot <- "\\footnote{Mean Computation Time of one replication in seconds.}"

# Filenames for the outputfiles
j1file <- "Spread/Tables/j1table.tex"
j3file <- "Spread/Tables/j3table.tex"


setwd(oldwd)

# Place the table
capture.output({

  print(xtable(j1resclean, digits=2, align = rep("c", 10)), floating = FALSE,
        size="footnotesize", #Change size; useful for bigger tables
        include.rownames=FALSE, #Don't print rownames
        include.colnames=FALSE, #We create them ourselves
        caption.placement="top",
        hline.after=NULL, #We don't need hline; we use booktabs
        add.to.row = list(pos = as.list(c(-1, -1,
                                          spaceLocs,
                                          nrow(j1resclean),
                                          nsplits)),
                          command = c(paste("\\toprule \n &&& \\multicolumn{2}{c}{Posterior $\\mu$}",
                                            "& \\multicolumn{2}{c}{Posterior $\\kappa$", kapfoot, "}& ",
                                            "\\\\\n  "), #\\cline{6-9}
                                      paste("n & $\\kappa$ & Method & Mean & Coverage & Mode & Coverage
                                            & Acc.", accfoot, "& MCT", mctfoot,"\\\\\n",
                                            "\\midrule \n"),
                                      rep("\\vspace{0.2cm}", (nrow(j1resclean)/nmeth)-3),
                                      "\\bottomrule \n",
                                      rep("\\cmidrule{1-9}", length(nsplits)))
        )
  )

}, file=j1file)





capture.output({

  print(xtable(j3resclean, digits=2, align = rep("c", 12)), floating=FALSE,
        size="footnotesize", #Change size; useful for bigger tables
        include.rownames=FALSE,
        include.colnames=FALSE, #We create them ourselves
        caption.placement="top",
        hline.after=NULL, #We don't need hline; we use booktabs
        add.to.row = list(pos = as.list(c(-1, -1,
                                          spaceLocs,
                                          nrow(j3resclean),
                                          nsplits)),
                          command = c(paste("\\toprule \n &&& \\multicolumn{4}{c}{Posterior $\\mu$}",
                                            "& \\multicolumn{2}{c}{Posterior $\\kappa$", kapfoot, "}& ",
                                            "\\\\\n  "), #"\\cline{4-7}"),
                                      paste("$n_j$ & $\\kappa$ & Method & $\\mu_1$ & $\\mu_2$ & $\\mu_3$ & Coverage & Mode & Coverage",
                                            "& Acc.", accfoot, "& MCT", mctfoot,"\\\\\n",
                                            "\\midrule \n"),
                                      rep("\\vspace{0.2cm}", (nrow(j3resclean)/nmeth)-3),
                                      "\\bottomrule \n",
                                      rep("\\cmidrule{1-11} \n ", length(nsplits)))
        )
  )

}, file=j3file)




j3resbigmoredigits <- comparisontable(rl=rlj3, J=3, fnames=c("DW", "MH", "FM"), ks=c("0.1", "4", "32"),digits=6)

# Average coverages per group
j3resmoredigits <- cbind(j3resbigmoredigits[, 1:6], round(rowMeans(apply(j3resbigmoredigits[, 7:9], 2, as.numeric)), 2), j3resbigmoredigits[, 10:13])


MHcov <- data.frame(n=rep(c(10, 30, 100), each=9),
                    kappa=factor(paste("kappa =", rep(c(0.1, 4, 32), times=3, each=3)), levels = paste("kappa =", c(0.1, 4, 32))),
                    #paste("kappa =", rep(c(0.1, 4, 32), times=3, each=3)),
                    Method=c("Gibbs", "MH", "Rejection"),
                    Coverage=as.numeric(j3resmoredigits[, c("Kap true in HDI")]), StringsAsFactors=FALSE)

MHcov[is.nan(MHcov[,4]),4] <- NA

require(ggplot2)
require(grid)

scale <- 0.8
pdf("Spread/Figures/Coverages.pdf", width=12*scale, height=5.5*scale)

ggplot(data=MHcov, aes(y=Coverage, x=n, shape=Method, group=Method)) +
  geom_line() + geom_point(size=3) + xlab("n") + ylab("Coverage") + theme_classic() +
  ggtitle("") +
  facet_wrap(~kappa) + geom_line(y=.95, color="black") +
  theme(panel.margin = unit(2, "lines")) +
  scale_x_continuous(breaks=c(10, 30, 100))
dev.off()


pdf("Spread/Figures/CoveragesWide.pdf", width=12*scale, height=3.5*scale)

ggplot(data=MHcov, aes(y=Coverage, x=n, shape=Method, group=Method, col=Method)) +
  geom_line() + geom_point(size=3) + xlab("n") + ylab("Coverage") + theme_classic() +
  ggtitle("") +
  facet_wrap(~kappa) + geom_line(y=.95, color="black") +
  theme(panel.margin = unit(2, "lines")) +
  scale_x_continuous(breaks=c(10, 30, 100)) +
  scale_color_manual(values = c("black", "red", "green"))
dev.off()


scale <- 4
pdf("Spread/Figures/ExampledatasetDeg.pdf", width=scale, height=scale)
plotCircular(rvmc(15, pi/4, 3), outputtype = "degrees")
dev.off()

scale <- 4
pdf("Spread/Figures/ExampledatasetDeg2.pdf", width=scale, height=scale)
plotCircular(rvmc(15, pi/3, 7), outputtype = "degrees")
dev.off()

scale <- 4
pdf("Spread/Figures/Exampledataset1.pdf", width=scale, height=scale)
plotCircular(rvmc(15, 2, 4))
dev.off()


pdf("Spread/Figures/Exampledataset2.pdf", width=scale, height=scale)
plotCircular(rvmc(15, 1.3, 4))
dev.off()

