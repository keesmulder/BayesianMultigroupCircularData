# ----------------------------------------------------------
# GenerateData.R
# Function to create a given number of datasets, on which then simulations can
# be run.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the
# Dutch Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------


library(Rcpp)
sourceCpp("Data/rvmc.cpp")


generateVonMisesData <- function(ndat, ns, kaps, Js=1, meandif=(20*(pi/180))) {
  # FUNCTION generateVonMisesData ---------------------
  #Args:
  # ndat:    Number of datasets to be generated.
  # ns:      Vector of possible sample sizes.
  # kaps:    Vector of possible kappas.
  # Js:      Vector of possible numbers of groups.
  # meandif: Difference between each group.


  for (J in Js) {

    # The true means.
    true_means <- (meandif*1:J)%%(2*pi)

    for (n in ns) {

      for (kap in kaps) {

        # Create a folder for the datasets
        dirname <- paste0(getwd(), "/Data/Datasets/Datasets_", "J=", J, "_n=", n, "_kap=", kap)
        dir.create(dirname, showWarnings=FALSE)

        for (i in 1:ndat){

          output_filename <- paste0(dirname, "/nr", i, ".csv")

          # Generate data from the von Mises distribution
          dat <- sapply(true_means, function(x) rvmc(n, x, kap))

          # Write data
          write.table(dat, file=output_filename, sep=",", row.names=FALSE, col.names=FALSE)

        }
      }
    }
  }
}


