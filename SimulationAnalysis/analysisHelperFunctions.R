require(abind)

# Load ssres from working directory, and extract the results part without Design. 
loadssres <- function(filename, wd="(4) Code for simulation/Results/") {
  if (wd!="") setwd(wd)
  load(filename)
  ssres$Results
}

# Retrieve the summarized results for each of a list of results from simulation
# studies. nout is the index of the outcome, meanfun is the function to be
# applied for summarizing the data (usually colMeans).
getres <- function(rl, nout, meanfun = colMeans){
  print(dimnames(rl[[1]])$outcomes[nout])
  sapply(rl, function(x) meanfun(x[, nout, , , ]))
}


# Combine several outcomes from their string filenames on a given dimension (ie. kappa or n)
combineRes <- function(rnames, onDim=3, wd="(4) Code for simulation/Results/") {
  setwd(wd)
  ssreslist <- lapply(rnames, function(x) {load(x); ssres$Results})
  abind(ssreslist, along=onDim)
}

# Create a table in which multiple samplers are compared. 
comparisontable <- function(rl, fnames, J=1, digits=2, ns=NULL, ks=NULL, concise=TRUE) {
  if (concise) {
    angInd <- if(J==1) {c(4)} else {8:10}
    linInd <- if(J==1) {c(7, 13, 21, 23, 25)} else {c(17:19, 27, 35, 37, 39)}
  } else {
    angInd <- if(J==1) {c(4)} else {8:10}
    linInd <- if(J==1) {c(7:8, 11:15, 20:21, 23:25)} else {c(17:20, 23:27, 32:33, 35:37)}
  }
  
  if (is.null(ns)) ns <- dimnames(rl[[1]])[[3]]
  if (is.null(ks)) ks <- dimnames(rl[[1]])[[4]]
  
  out <- data.frame(matrix(nc=3+length(angInd)+length(linInd), nr=length(rl)*length(ns)*length(ks)))
  ri <- 1
  
  for (n in ns) {
    dispn <- n
    
    for (k in ks) {
      dispk <- k
      
      for (fname in fnames) {        
        r <- rl[[fname]]
        
        # Find results for the angles, as they need special treatment.
        angleResults  <- round(apply(r[, angInd, n, k, , drop=FALSE], 2, meanDir, na.rm=TRUE)*(180/pi), digits=digits)
        
        linearResults <- round(colMeans(r[,linInd,n,k,], na.rm=TRUE), digits=digits)
        
        out[ri, ] <- c(dispn, dispk, fname, angleResults, linearResults)
        
        # Set for next iteration.
        dispn <- dispk <- ""
        ri <- ri+1
      }
    }
  }
  colnames(out) <- c("$n$", "$\\kappa$", "Method", colnames(r)[angInd], colnames(r)[linInd])
  out
}

# Switch out missingness for dashes, set names to their longer versions,
# addZeros where necessary.
cleanRes <- function (x, addZeros=TRUE) {
  x[is.na(x) | as.matrix(x)=="NaN" | x==0] <- "---"
  x[x=="DW"] <- "Gibbs"
  x[x=="FM"] <- "Rejection"
  if (addZeros) x[, 4:ncol(x)] <- sapply(x[, 4:ncol(x)], addZeros)
  x
}

# Adds forced double digits to a string vector containing strings representing
# numbers as well as non-numeric strings.
addZeros <- function(chvec, digits=2) {
  sapply(chvec, function (ch) {
    dotloc <- gregexpr("[.]", ch)[[1]][1]
    nadd   <- digits - nchar(ch) + dotloc
    paste0(ch, ifelse(nadd>0, rep(0, nadd), ""))
  })
}



# Get datasets for testing errors.
thget <- function(i, J=3, kappa=.1, n=5, wd=paste0(getwd(), "/(2) Data/Datasets")) {
  if (wd!="") setwd(wd)
  readfilename <- paste0("/Datasets_", 
                         "J=", J, "_n=", n, "_kap=", kappa, "/nr", i, ".csv")
  th <- read.table(readfilename, sep=",")
  th
}
