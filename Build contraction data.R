
#This script analyzes fluorescence kymographs to measure 
#and visualize the Contraction or diameter changes over time. 
#It involves:
# 1. Reading and processing TIFF images.
# 2. Computing the diameter of features using autocorrelation.
# 3. Normalizing and plotting the results.
# 4. Saving the results to files for further analysis.

############ Libraries
library(tiff)     #tiff: This package is used to handle TIFF image files, specifically to read kymographs stored in TIFF format.
############ Functions: : This function measures the contraction or diameter of a feature in a kymograph by analyzing the autocorrelation function (ACF) of the data.

measureContraction<-function(f, upL, lowL, ys, Time){   ##Parameters:
                                                        # f: contraction data (one row of the kymograph).
                                                        #upL: Upper limit for contraction measurement.
                                                        #lowL: Lower limit for contraction measurement.
                                                        #ys: Spatial positions in micrometers.
                                                        #Time: Time vector.
  
  f.s<-smooth.spline(ys,f, df = DF)$y                   #Smooth Data: smooth.spline is used to smooth the input data (f) across the y-values (ys) with a specified degrees of freedom (DF).
  a<-acf(f, lag =  (upL/pixelSize)*magnification, plot = F) #Compute ACF: Compute the autocorrelation function (ACF) of the smoothed data to analyze periodicity with specified lags.

  dy<-a$lag*pixelSize/magnification         #The distance (dy) for each lag is calculated based on the pixel size and magnification.
  id<-which(dy>lowL & dy <upL)              #Find Relevant Peaks: Extract lag values and corresponding ACF values within the specified range (lowL to upL). The function checks if this maximum is valid (i.e., it should not be at the boundary or not a local maximum).
  ac<-a$acf[id]
  mx.id<-which.max(ac)
  if(mx.id == 1 | mx.id==length(ac) ){
    return(NA)
    } else if(ac[mx.id]<ac[mx.id-1] | ac[mx.id]<ac[mx.id+1] ){
    return(NA)
    }
  d<-which(a$acf==ac[mx.id])*pixelSize/magnification
  if(d == lowL | d == upL) d <- NA
  return(d)
}                                     #Determine contraction: Find the lag where the ACF is maximum and check if it is within valid bounds. If so, return the corresponding diameter.
computeDiameter<-function(ks, upL, lowL,PLOT=T){        #Computes the diameter for each row in the kymograph (ks) and optionally plots the kymograph with the computed diameter data.
  ################ measure Contraction
  ys<-(1:dim(ks)[2])*pixelSize/magnification
  Time<-(1:dim(ks)[1])*acTime       #Define Axes: ys and Time vectors are created based on the dimensions of the input kymograph.
  d<-NULL
  dim(ks)
  for(i in 1:nrow(ks)){
    d[i]<-measureContraction(ks[i,], upL, lowL, ys, Time)   #Measure contraction: For each row of the kymograph, the measureContraction function is called to compute the diameter.
  }
  if(PLOT){
    image(x = Time, y = ys, z= ks, 
          col = lutRed(255),
          ylab = "diameter (um)",
          xlab = "time (min)",
          main = kName)
    abline(v=3, col = "white", lwd = 2, lty=2)
  }
  return(d)
}
plotdiameter<-function(s, s.norm, acTime, Name){    # Plots the diameter and normalized diameter data.
  Time<-(1:dim(s)[1])*acTime
  colnames(s)<-colnames(s.n)<-l
  par(mfrow = c(1,2))
  plot(s[,1]~Time, main = Name,
       ylab = "diameter (um)",
       xlab = "time (min)",
       pch = 19, col = "black",
       ylim =c(0.5, 1))
  abline(v=3, col = "black", lwd = 2, lty=2)
  for(k in 2:length(l)){
    points(s[,l[k]]~Time,pch=19, col= k)
    points(s[,l[k]]~Time,pch=19, col= k)
  }
  legend("topright", legend=l, col = 1:length(l), cex =0.5, lwd = 3)
  plot(s.n[,1]~Time, main = paste(Name, "norm"),
       ylab = "norm diameter",
       xlab = "time (min)",
       pch = 19, col = "black",
       ylim =c(0.5, 1.4))
  abline(v=3, col = "black", lwd = 2, lty=2)
  for(k in 2:length(l)){
    points(s.n[,l[k]]~Time,pch=19, col= k)
    points(s.n[,l[k]]~Time,pch=19, col= k)
  }
  legend("topright", legend=l, col = 1:length(l), cex =0.5, lwd = 3)
}
############ parameters      Purpose: Defines various constants used in the analysis. These include pixel size, magnification, acquisition time, degrees of freedom for spline smoothing, and limits for measuring stretching.
pixelSize <- 0.117 # um       # Micrometer per pixel
magnification <- 5 # times    # Magnification factor
acTime <- 0.5 #min            # Acquisition time in minutes
DF <- 8                       # Degrees of freedom for smoothing spline
upL = 1                       # Upper limit for contraction measurement
lowL = 0.3                    # Lower limit for contraction measurement
############ LUTs      Defines a color lookup table (LUT) for creating a gradient from black to red to yellow for the plots.
lutRed <- colorRampPalette(c("#000000", "#ff0000", "#ffff00"))
############ Load fm4-64 fluorescence kymograph
wd<-"C:/.."                   # Working directory
setwd(wd)                     # Set working directory
cdir<-list.dirs(wd)           # List all directories in the working directory


############ Compute diameter from fluorescence kymograph 

for(cDir in cdir[-1]){
  print(cDir)
  l<-list.files(cDir, pattern = "um.tif")   # List TIFF files in the directory
s<-s.n<-NULL
for(kName in l){
  print(kName)
  ks<-readTIFF(paste(cDir,kName, sep ="/")) # Read TIFF image
  d<-computeDiameter(ks, upL, lowL, PLOT = T) # Compute diameter
  d.n<-d/mean(d[1:5], na.rm=T)                # Normalize diameter
  s<-cbind(s,d)
  s.n<-cbind(s.n, d.n)
}
colnames(s)<-colnames(s.n)<-l
write.table(s, paste(cDir, "s.txt", sep ="_"))  # Save diameter data
write.table(s.n, paste(cDir, "s.n.txt", sep ="_"))   # Save normalized diameter data
plotdiameter(s, s.norm, acTime, "0.1P")     # Plot the diameter data
}

