source("(6) Basic codes/DescribeCirc.R")

plotExampleRMu <- function(th, r=0.4) {
  
  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")
  require(plotrix)  
  
  # Draw empty plot
  plot(-1:1, -1:1, type = "n", 
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), 
       xlab = "", ylab = "",
       asp = 1, xaxt = "n", yaxt = 'n', bty = 'n')
  
  # Draw empty circle
  draw.circle(0,0, radius = r)
  
  # Draw 0 degree bar.
  segments(0,0, r, 0, col = "black")
  text(1.15*r, 0, labels = "0", col = "black")
  
  
  # Compute mean direction
  th_bar <- meanDir(th)

  # Get sample size
  n <- length(th)
  
  # Start loc at the origin. 
  origin <- c(0, 0)
  loc <- origin
  
  # Loop through the angles.
  for (i in 1:n) {
    
    # Save the cosine and sine moment of current angle.
    csthi <- r * c(cos(th[i]), sin(th[i]))
    
    # Draw this point as an arrow and a point on the unit circle. 
    points(csthi[1], csthi[2], cex=0.7, pch=16, col="gray")
    arrows(origin[1], origin[2], csthi[1], csthi[2], length = 0.15 * r, col="gray80", lwd=1.8)
    
    # Draw the arrow that gives R.
    newloc <- loc + csthi
    arrows(loc[1], loc[2], newloc[1], newloc[2], length = 0.10 * r, lwd=0.6)
    loc <- newloc
  }
  
  # Draw a line to (C, S).
  lines(col="gray1", rbind(origin, newloc), lty="dotted")
  
  # Place the text R. 
  Rloc <- 0.7 * newloc + c(-1, 1) * 0.05 * newloc
  text(Rloc[1], Rloc[2], labels="R")

  # Draw a part of a circle to the mean direction. 
  pts <- seq(from=0, to=th_bar*(180/pi), length.out=100)*(pi/180)
  cdp <- 0.5*r*cos(pts)
  sdp <- 0.5*r*sin(pts)
  lines(cdp, sdp, col="gray")  
  
  # Add a label for the mean direction.
  text(cos(th_bar/2)*0.65*r, sin(th_bar/2)*0.65*r, labels=expression(bar(theta)))
}


th <- c( 1.35, -0.28, 0.97)

scale <- 0.55
pdf("ExampleRMu.pdf", width=12*scale, height=8*scale)
plotExampleRMu(th, r=0.8)
dev.off()

round(th%%(2*pi)*(180/pi), 0)

