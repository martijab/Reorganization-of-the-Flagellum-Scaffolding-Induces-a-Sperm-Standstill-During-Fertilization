#This script processes a TIFF image file containing multiple frames. It reads the frames, 
#performs thresholding to create binary images, and calculates the fraction of overlapping 
#pixels between two measures (fm and sa). It then plots these fractions over time and saves 
#the results to a text file. The use of thresholding helps in analyzing features in the images 
#by converting pixel values to binary, facilitating easier comparison and quantification.


library(tiff)   #For reading TIFF image files.
library(fields) #For plotting and spatial data analysis (though image.plot is the main function used here).

#Sets the working directory to the specified path where your TIFF image file (36NR.tif) is located. 
#This allows you to load files from this directory.
wd<-"D:/.."
setwd(wd)
#Reads the TIFF file into a list of matrices. Each element of img represents a separate image frame from the TIFF file.
img<-readTIFF("example.tif",all = T)

length(img)   #Returns the number of frames in the TIFF file.
image.plot(img[[1]])
image.plot(img[[2]])
image.plot(img[[3]])    #Plots the first three frames of the TIFF file for visual inspection.

#Create sequences of indices to select frames from img. idfm selects frames for fm (FM4-64), 
#and idsa selects frames for sa (SiR-actin).
idfm<-seq(2,length(img), by=3)
idsa<-seq(3,length(img), by=3)

fm<-sa<-cf<-NULL    #Initialize empty lists to store these frames.
fm<-sa<-list()
count<-2
for(i in idfm){
  fm[[count]]<-img[[i]]
  sa[[count]]<-img[[i+1]]
  count<-count+1
}               #Loop through the indices idfm to populate the fm and sa lists. 
                #Each fm frame corresponds to a sa frame immediately following it.

#Convert fm and sa frames to binary images based on thresholds th1 (for fm) and th2 (for sa).
#These thresholds act as background subtraction.
th1<-180
  th2<-125
  cffm<-cfsa<-NULL
for(i in 1:length(fm)){
  idx<-which(fm[[i]]<th1)
  fm[[i]][idx]<-0
  idx2<-which(fm[[i]]>=th1)
  fm[[i]][idx2]<-1
  
  idxs<-which(sa[[i]]<th2)
  sa[[i]][idxs]<-0
  idxs2<-which(sa[[i]]>=th2)
  sa[[i]][idxs2]<-1
  
  a<-fm[[i]]+sa[[i]]      # Sum the binary images to create a combined image.
  
  idcor<-which(a==2)    # Indices where both fm and sa are 1 (colocalization).
 
   totfm<-which(fm[[i]]>0)
   totsa<-which(sa[[i]]>0)
   
   
  cffm[i]<-(sum(idcor)/sum(totfm))

  cfsa[i]<-(sum(idcor)/sum(totsa))
  
    }           #Calculate the fraction of overlapping pixels for each frame, manders colocalization coefficients.
                # cffm refers to M1 (Assigned to FM4-64, proportion of F-actin that colocalizes with plasma membrane)
                # cfsa refers to M2 (Assigned to SiR-actin, proportion of plasma membrane that colocolizes with F-actin)

time<- (1:length(cffm))

plot(cffm, x=time, xlab="", ylab="",              #Plot cffm values over time with red points.
     pch=19, col="red", ylim=c(0,1), cex.axis=1.5, cex.lab=1.5)
points(cfsa, x=time, pch=19, col="blue")          #Add cfsa values to the plot with blue points.
abline(v=3, col="black", lty=3, lwd=4)

cffm
cfsa
manders<-cbind(cffm,cfsa) #Combine cffm and cfsa into a matrix.


write.table(manders, file="example.txt")   # Save the matrix to a text file for further analysis.

