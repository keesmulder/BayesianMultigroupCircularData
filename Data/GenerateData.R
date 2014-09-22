library(Rcpp)
sourceCpp("Data/rvmc.cpp")

# Function to generate datasets
# ndat:    Number of datasets to be generated.
# ns:      Vector of possible sample sizes.
# kaps:    Vector of possible kappas.
# Js:      Vector of possible numbers of groups.
# meandif: Difference between each group.

generateVonMisesData <- function(ndat, ns, kaps, Js=1, meandif=(20*(pi/180))) {


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

# Generate the data.
set.seed(8921747)
generateVonMisesData(10, c(5, 30, 100), c(0.1, 4, 32), J=c(1,3))

