library(Rcpp)
sourceCpp("(6) Basic codes/rvmc.cpp")

# Function to generate datasets
generateData <- function(ndat, ns, kaps, J=1, meandif=(20*(pi/180))) {
  
  # The true means.
  true_means <- (meandif*1:J)%%(2*pi)
  
  for (n in ns) {
    for (kap in kaps) {
      
      # Where to create data
      dirname <- paste0(getwd(), "/(2) Data/Datasets/Datasets_", "J=", J, "_n=", n, "_kap=", kap)
      dir.create(dirname, showWarnings=FALSE)
      
      for (i in 1:ndat){
        outputfile <- paste0(dirname, "/nr", i, ".csv")
        
        # Get data
        dat <- sapply(true_means, function(x) rvmc(n, x, kap))
        
        # Write data
        write.table(dat, file=outputfile, sep=",", row.names=FALSE, col.names=FALSE)
        
      }
    }
  }
}

set.seed(8921747)
generateData(1000, c(5, 30, 100), c(0.1, 4, 32), J=1)
generateData(1000, c(5, 30, 100), c(0.1, 4, 32), J=3)