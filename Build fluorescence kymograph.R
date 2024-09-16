
##This script processes fluorescence kymographs, normalizes, smooths, and plots them,
##then saves the results.

############ Libraries
library(tiff)       ##Handles TIFF image files.
library(fields)     ##Provides tools for visualization and analysis of spatial and spatio-temporal data, including image.plot.
library(stringr)    ##Simplifies string manipulation.

########### functions
computeKyms<-function(kName, bName, reacTime, channel){
 #Read the TIFF file:readTIFF loads the image stack; K[[channel]] extracts the specified channel.
   K<-readTIFF(kName,all = T, as.is=T)   ##kName: File name for the kymograph TIFF image.
  kf<-K[[channel]]                      ##channel: Channel to use from the TIFF image.
  
 #Define spatial and time axes: xf (spatial) and Time (temporal) based on image dimensions.
   xf<-(1:dim(kf)[2])*pixelSize
  Time<-(1:dim(kf)[1])*acTime
  idx<-which(Time<reacTime)             ## idx defines the previous frames to the induction of acrosome reaction (reacTime)

  #  Load and subtract background:
  bckg<-read.table(bName)$V1           ##Background data is read from a .csv file and used to correct the kymograph.
                                       ##bName: Background data file name.
  B<-read.csv(bName, sep = ";")
  bckg<-B[["fm"]]
  kf.b<-kf.b.n<-NULL                   ##kf.b is the background-subtracted kymograph
  for(i in 1:ncol(kf)){                ##kf.b.n is the normalized kymograph.
    a<-kf[,i]-bckg
    kf.b<-cbind(kf.b, a)
    kf.b.n<-cbind(kf.b.n, a/mean(a[idx]) - 1)  ##Normalize by dividing by the mean of the initial values and then subtract 1.
  }
  ###### build normalized and smoothened kymogram
  kf.ns<-array(0,dim = dim(kf.b.n))
  db<-array(0, dim = c(nrow(kf.b.n)-1, ncol(kf.b.n)))
  for(d in 1:ncol(kf.b.n)){
    fNorm<-kf.b.n[,d]     
    kf.ns[,d]<-smooth.spline(fNorm, df = DF)$y    ##Smooth the normalized kymograph using smooth.spline.
    db[,d]<-diff(kf.ns[,d])                       ##Compute the derivative of the smoothed data (db).
  }
  K<-list(K = K, B= B, f =kf, f.b= kf.b, f.b.n = kf.b.n, f.ns = kf.ns, f.db = db)
}                                                 ##Return a list with all processed data.

plotKyms<-function(k, NAME){        ##k: List of processed kymograph data.
                                    ##NAME: Name for the plot title.
  xf<-(1:dim(k$f)[2])*pixelSize
  Time<-(1:dim(k$f)[1])*acTime
  idx<-which(Time<reacTime)
  
  par(mfrow = c(2,3))             ##Set plotting layout: par(mfrow = c(2, 3)) sets up a 2x3 grid for plotting.
  image.plot(x = Time, y = xf, z= k$f, 
             col = lutRed(255),
             ylab = "Length (µm)",
             main = NAME, xlab = "Time (min)")
  abline(v=3, col = "white", lwd = 2, lty=2)
  image.plot(x = Time, y = xf, z= k$f.b, 
             col = lutRed(255),
             ylab = "Length (µm)",
             main = "Background substracted",xlab = "Time (min)")
  abline(v=3, col = "white", lwd = 2, lty=2)
  image.plot(x = Time, y = xf, z= k$f.b.n, 
             col = lutRed(255),
             ylab = "Length (µm)",
             main = "Normalized",xlab = "Time (min)")
  abline(v=3, col = "white", lwd = 2, lty=2)
  image(x = Time, y = xf, z= k$f.ns, 
        col = lutRed(255),
        ylab = "Length (µm)",
        main = "Smoothed",xlab = "Time (min)")
  abline(v=3, col = "white", lwd = 2, lty=2)
##Plot images: image.plot and image functions plot the original, background-subtracted, 
## normalized, smoothed, and derivative kymographs.
  contour(x = Time, y = xf, z= k$f.ns, col = "white", add = T, nlevels = 2)
##Add contours and lines for visual reference.
  abline(h=c(7,14), col = "grey", lwd = 1, lty=2)
  plot(k$f.ns[nrow(k$f.ns),],xf, pch = 19, main = "F/F0-1 last frame",
       ylab = "Length (µm)", xlab = "F/F0-1")
  abline(v=1, col = "grey", lwd = 2, lty=2)
  abline(h=c(7,14), col = "grey", lwd = 1, lty=2)
  image.plot(x = Time[-1],y = xf, z=k$f.db, main = "(dFs/dt)/F0",
             ylab = "Length (µm)", xlab = "Time (min)")
  abline(v=3, col = "black", lwd = 2, lty=2)
  abline(h=c(7,14), col = "black", lwd = 1, lty=2)
}
############ parameters
pixelSize <- 0.117 # um         #Size of a pixel in micrometers.      
magnification <- 5 # times      #Magnification factor.
acTime <- 0.5 #min              #Time per frame in minutes.
DF <- 8                         #Degrees of freedom for smoothing spline.
reacTime <- 2                   #Time of induction of acrosome reaction
channel <- 2                    #The channel index to process.
############ LUTs (look-up table)
lutRed <- colorRampPalette(c("#000000", "#ff0000", "#ffff00"))    ##Defines a color gradient from black to red to yellow for plotting.
############ Load fm4-64 fluorescence kymograph
wd<-"D:/.../"   #Set working directory: Defines where the data is stored and where results will be saved.
resdir<-paste(wd,"results2",sep="/")
dir.create(resdir)  #Create results directory: dir.create(resdir) ensures the output directory exists.
setwd(wd)
#Process each TIFF file
lk<-list.files(wd, pattern = ".tif")
for(kName in lk){
  print(kName)
  cell<-strsplit(kName,"_")[[1]][1]   #Extract file names.
  bName<-paste(cell,"_Fondo.csv", sep="")
  print(bName)
  k<-computeKyms(kName, bName, reacTime, channel) #Compute kymograph and background correction.
  save(k, file = paste(resdir,"/", cell,".RData", sep=""))  #Save processed data and generate plots.
  pdf(paste(resdir,"/", cell,".pdf", sep=""), height = 6, width = 8.5)
  plotKyms(k, NAME = cell)
  dev.off()
}

############## Read and plot a kymogram
load("D:/.../...RData")   #Load data: Loads a previously saved .RData file.
image.plot(k$f.b)
plotKyms(k, "01")         #Plot background-subtracted kymograph and call plotKyms for a specific plot.
