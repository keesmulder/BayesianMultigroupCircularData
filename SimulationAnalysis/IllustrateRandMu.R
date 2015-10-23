# ----------------------------------------------------------
# analysisHelperFunctions.R
# Plot an example of what R and the mean direction represent.
#
# Kees Tim Mulder
# Last updated: November 2014
#
# This work was supported by a Vidi grant awarded to I. Klugkist from the Dutch
# Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------


source("DataAnalysis/DescribeCirc.R")
require(circular)


plotExampleRMu <- function(th, r=0.4, col=FALSE) {
  # FUNCTION plotExampleRMu ------------------------------------------------
  # th: An arbitrary set of data for which to create the visualization.
  # r: The radius, or size of the picture.
  # ------------------------------------------------------------------------
  
  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")
  require(plotrix)
  
  if (col) {
    cols <- c("tomato", "darkolivegreen", "skyblue")
  }
  
  # Draw empty plot
  plot(-1:1, -1:1, type = "n",
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
       xlab = "", ylab = "",
       asp = 1, xaxt = "n", yaxt = 'n', bty = 'n')
  
  # Draw empty circle
  draw.circle(0,0, radius = r, border="gray50", lwd = 1.8)
  
  # Draw 0 degree bar.
  segments(0,0, r, 0, col = "black")
  text(1.1*r, 0, labels = "0", col = "black")
  
  
  # Compute mean direction
  th_bar <- meanDir(th)
  
  # Get sample size
  n <- length(th)
  
  # Start loc at the origin.
  origin <- c(0, 0)
  loc <- origin
  
  # Loop through the angles.
  for (i in n:1) {
    
    # Select colors.
    if (col) {
      col1 <- col2 <- col3 <-  cols[i]
    } else {
      col1 <- "gray"
      col2 <- "gray80"
      col3 <- "black"
    }
    
    # Save the cosine and sine moment of current angle.
    csthi <- r * c(cos(th[i]), sin(th[i]))
    
    # Draw this point as an arrow and a point on the unit circle.
    # points(csthi[1], csthi[2], cex=0.7, pch=16, col=col1)
    arrows(origin[1], origin[2], csthi[1], csthi[2], length = 0.13 * r,
           col=col2, lwd=1.5)
    
    # Draw the arrow that gives R.
    newloc <- loc + csthi
    arrows(loc[1], loc[2], newloc[1], newloc[2],
           length = 0.13 * r, lwd=1.5, col=col3)
    loc <- newloc
  }
  
  # Draw a line to (C, S).
  lines(col="gray1", rbind(origin, newloc), lty="dashed")
  
  # Place the text R.
  Rloc <- 0.7 * newloc + c(-1, 1) * 0.07 * newloc
  text(Rloc[1], Rloc[2], labels="R")
  
  # Draw a part of a circle to the mean direction.
  pts <- seq(from=0, to=th_bar*(180/pi), length.out=100)*(pi/180)
  cdp <- 0.5*r*cos(pts)
  sdp <- 0.5*r*sin(pts)
  lines(cdp, sdp, col="gray50")
  
  # Add a label for the mean direction.
  text(cos(th_bar/2)*0.6*r, sin(th_bar/2)*0.6*r, labels=expression(bar(theta)))
}

#
# th <- c(1.7, -0.28, 1.1)-.2
#
# scale <- 0.6
# pdf("ExampleRMu.pdf", width=12*scale, height=8*scale)
# plotExampleRMu(th, r=0.7, col=TRUE)
# dev.off()
# # plotExampleRMu(rvonmises(50, 3, 4), r=0.3)
# # plotExampleRMu(c( 1.35, -3.28, 1.97), r=0.5)
#
#
#




























plotLinearCircMean <- function(th, r=0.4, col=FALSE) {
  
  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")
  require(plotrix)
  
  if (col) {
    cols <- c("tomato", "darkolivegreen", "skyblue", "tan")
  }
  
  # Draw empty plot
  plot(-1:1, -1:1, type = "n",
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
       xlab = "", ylab = "",
       asp = 1, xaxt = "n", yaxt = 'n', bty = 'n')
  
  # Draw empty circle
  draw.circle(0,0, radius = r, border="gray50", lwd = 1.8)
  
  # Draw 0 degree bar.
  segments(0,0, r, 0, col = "black",  lty = "solid", lwd=2.5)
  text(1.4*r, 0, labels = "Circulair", col = "black")
  
  segments(0,0, -r, 0, col = "grey30", lty = "longdash", lwd=1.7)
  text(-.60*r, .1, labels = "Lineair", col = "black")
  
  
  # Compute mean direction
  th_bar <- meanDir(th)
  
  # Get sample size
  n <- length(th)
  
  # Start loc at the origin.
  origin <- c(0, 0)
  loc <- origin
  
  # Loop through the angles.
  for (i in n:1) {
    
    # Select colors.
    if (col) {
      col1 <- col2 <- col3 <-  cols[i]
    } else {
      col1 <- "gray"
      col2 <- "gray80"
      col3 <- "black"
    }
    
    # Save the cosine and sine of current angle.
    csthi <- r * c(cos(th[i]), sin(th[i]))
    
    # Draw this point as an arrow and a point on the unit circle.
    # points(csthi[1], csthi[2], cex=0.7, pch=16, col=col1)
    arrows(origin[1], origin[2], csthi[1], csthi[2], length = 0.13 * r,
           col=col2, lwd=1.5)
    
    thideg <- th[i]*180/pi
    text(csthi[1]*.94, csthi[2]*1.09, paste0(thideg, "°"), pos = 4, col="gray30")
    
  }
  
  
}


th <- c(1.7, -0.28, 1.1)-.2
th <- c(10, 30, 330, 350)*pi/180

scale <- 0.7
pdf("LinearCircMean.pdf", width=12*scale, height=8*scale)
plotLinearCircMean(th, r=0.7, col=TRUE)
dev.off()
legend(1, 1, 
       legend = c(expression(paste(kappa, " = 3")),
                  expression(paste(kappa, " = 1")),
                  expression(paste(kappa, " = 0"))), 
       col    = c("tomato", 
                  "skyblue", 
                  "darkolivegreen"),
       lty=1, cex=2)

# 
# 
# require(circular)
# 
# plotCircularDensity(function(x) dvonmises(x, 2, 3))
# plotCircularDensity(function(x) dvonmises(x, 1, 10))
# plotCircularDensity(function(x) dvonmises(x, 1, 2))
# 
# 
# 






























