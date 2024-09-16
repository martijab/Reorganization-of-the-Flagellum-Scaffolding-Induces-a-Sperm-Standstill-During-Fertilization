#The script processes TIFF images to extract and analyze local maxima, normalizes the distances
#from the center of the flagellum, and visualizes these distances over time. It also performs 
#linear regression to understand trends and saves the results for further analysis.


library(tiff)   #To read TIFF image files.
library(fields) #Provides tools for visualizing spatial data, including image.plot.
wd<-"C:/..."
setwd(wd)       #Sets the working directory to where the TIFF image file is located.
img<-readTIFF("example.tif",all = T)  # Reads the TIFF file into a list of matrices where each matrix
                                      #corresponds to an individual frame from the TIFF file.


idfm<-seq(2,length(img), by=3)
idsa<-seq(3,length(img), by=3)
#Generates sequences to select frames. idfm selects frames starting from the second frame and 
#skipping every 3 frames, while idsa starts from the third frame.

fm<-sa<-cf<-NULL     #Initializes empty lists fm and sa. 
                     #The variable count starts at 2 and is used to index these lists.
fm<-sa<-list()
count<-2

for(i in idfm){     ##Fills fm and sa lists with frames from the TIFF file. 
                    #Each fm frame is paired with the sa frame immediately following it.
  fm[[count]]<-img[[i]]
  sa[[count]]<-img[[i+1]]
  count<-count+1
}


d<-NULL                          #Sums the values in each row of frames in fm and sa, 
                                #creating matrices d and s where each row represents the summed 
                                #values of a corresponding frame.
for(i in 2:length(fm)){
  d<-rbind(d,apply(fm[[i]], 1, sum))
}

s<-NULL
for(i in 2:length(sa)){
  s<-rbind(s,apply(sa[[i]], 1, sum))
}

xs<-(1:ncol(d))*0.024   #xs represents pixel distances multiplied by 0.024 (1 pixel=0.024 micrometers).
time<-1:nrow(d)         #time is a sequence from 1 to the number of rows in d, representing time steps.


find2localmax<-function(f){   #Defines a function to find the two highest local maxima in a vector f. 
                              #It calculates differences between adjacent elements, identifies peaks, 
                              #and selects the two highest peaks.
 
   l<-sign(diff(f))
  posmax<-NULL
  for(i in 1:(length(l)-1)){
    if(l[i]>l[i+1]) posmax<-c(posmax,i)
  }
  m<- sort(f[posmax], decreasing=T)
  
  localmax<-c(which(f==m[1]),
              which(f==m[2]))
  return(localmax)
}


#Applies find2localmax to each row of matrices d and s to find local maxima.
dfm<-apply(d,1,find2localmax)
ssa<-apply(s,1,find2localmax)

#Plots the positions of local maxima found in dfm and ssa over time, 
#with red points representing dfm and blue points representing ssa.


#Calculate and Plot Center and Normalized Distances
plot(dfm[1,]*0.024, ylim=c(0.75,1.75), pch=19,col="red", xlab= "", ylab="")
title(ylab = "Distance (pixel)", line = 2.5, cex.lab = 0.9)
title(xlab = "Time (min)", line = 2.5, cex.lab = 0.9)
points(dfm[2,]*0.024,col="red", pch=19)

points(ssa[1,]*0.024, col="blue", pch=19)
points(ssa[2,]*0.024,col="blue", pch=19)



cent<-apply(dfm,2, mean )   # cent: Calculates the mean position of local maxima across frames.
                            #Normalizes the local maxima by subtracting this mean to center them.
points(cent*0.024, col="green", pch=15)

dfm1.norm<-dfm[1,]-cent   #dfm1.norm and dfm2.norm are distances from the center for the first and second sets of maxima in dfm.
dfm2.norm<-dfm[2,]-cent


ssa1.norm<-ssa[1,]-cent   #ssa1.norm and ssa2.norm are distances from the center for ssa.
ssa2.norm<-ssa[2,]-cent

sa.abs<-NULL            #sa.abs and fm.abs compute the average absolute distances from the center for ssa and dfm, respectively.
for(i in 1:length(ssa1.norm)){
  sa.abs[i] <- mean(abs(ssa1.norm[i]),abs(ssa2.norm[i]))
}
fm.abs<-NULL
for(i in 1:length(dfm1.norm)){
  fm.abs[i] <- mean(abs(dfm1.norm[i]),abs(dfm2.norm[i]))
}

#Plots normalized distances from the center over time, with red and blue points representing dfm and ssa, respectively.
plot(dfm2.norm*0.024, ylim=c(-0.5,0.5),col="red", pch=19, 
     xlab= "", ylab="")
title(ylab = "Distance to the center \n of flagellum (pixel)", line = 2.5, cex.lab = 0.9)
title(xlab = "Time (min)", line = 2.5, cex.lab = 0.9)
points(dfm1.norm*0.024,col="red", pch=19)
points(ssa1.norm*0.024, col="blue", pch=19)
points(ssa2.norm*0.024,col="blue", pch=19)
abline(h=0)



#Combines normalized distances into a table and saves it as a text file named example_maximas.txt.
table<-cbind(dfm1.norm*0.024,dfm2.norm*0.024,ssa1.norm*0.024,ssa2.norm*0.024)
write.table(table, file="example_maximas.txt")


#Linear Regression and Plotting

## Set plotting parameters for margins
par(mar=c(4.5,6,4.5,2))
## Create a plot of absolute deviations of flagellum maxima (FM)
plot(fm.abs [time]*0.024~time,ylim=c(0,0.5),col="red", pch=19,
     xlab="", ylab="")
title(ylab = "Distance to the center \n of flagellum (µm)", line = 2.5, cex.lab = 0.9)
title(xlab = "Time (min)", line = 2.5, cex.lab = 0.9)
# Add points for absolute deviations of flagellum maxima (SA)
points(sa.abs[time]*0.024~time, col="cyan", pch=19)
# Add a horizontal line at y = 0
abline(h=0, lty=2)

# Perform linear regression on the flagellum data (FM)
p<-data.frame(fm.abs*0.024, time)
linFM<-lm(p)
summary(linFM)
# Add regression line to the plot
abline(linFM)
# Add text annotation to the plot with the slope of the regression line
text(10, 0.4, labels= paste(round(linFM$coefficients[2],3),"µm/min"))
# Perform linear regression on the abs deviations of flagellum maxima (SA)
q<-data.frame(sa.abs*0.024,time)
linSA<-lm(q)
summary(linSA)
# Add regression line to the plot
abline(linSA)
# Add text annotation to the plot with the slope of the regression line
text(10, 0.2, labels= paste(round(linSA$coefficients[2],3),"µm/min"))


