library(plotrix)  #This library includes additional plotting functions not available in base R, such as plotCI, which is used for adding error bars to plots.
library(stats)    #This library includes fundamental statistical functions, including lm for fitting linear models, which is essential for regression analysis.

############ Functions   This function extracts and organizes specific data from contraction and florescence kymograph data.
readSingleContVSFnorm<-function(db){
  a<-NULL
  for(i in 1:dim("contraction data")[2]){
    a<-c(a, "contraction data"["data from rows starting from induction time to the end of the data frame."("contraction data"),i])
  }
  #Loop through each column of "contraction data" and concatenate data from rows starting from induction time to the end of the data frame.
  # a will contain the concatenated contraction values from all columns.
    b<-NULL
  for(i in R[1:dim("contraction data")[2]]){
    b<-c(b, "florescence kymograph data"["data from rows starting from induction time to the end of the data frame"("contraction data"),i])
  }
    #Similar to the first step, but uses the indices stored in R to select the appropriate columns from the fluorescence kymographs
    # b will contain the fluorescence values. 
     return(list(dn=a,fn=b))
}
#Returns a list with two elements: dn (data from contraction data) and fn (data from fluorescence data in kymographs).

############ parameters     Defines constants and calculates indices used in data processing.
pixelSize <- 0.117 # um
magnification <- 5 # times
#pixelSize and magnification: Define the physical size of a pixel and the magnification used in imaging. These are crucial for converting pixel measurements to real-world units

r<-(1:200)*pixelSize    #Creates a sequence of distances (in micrometers) from 1 to 200, scaled by pixelSize.
                        #These lines create logical indices for specific ranges of r. For example, r1 includes indices for distances between 2.45 and 2.57 micrometers.
                        #These ranges are used to select subsets of data for further analysis.
r0<-1 
r1<-which(r > 2.45 & r< 2.57)
r2<-which(r > 4.91 & r< 5.03)
r3<-which(r > 7.48 & r< 7.60)
r4<-which(r > 9.94 & r< 10.06)
r5<-which(r > 12.40 & r<12.51)
r6<-which(r > 14.97 & r<15.09)
R<-c(r0,r1,r2,r3,r4,r5,r6)    #Combines all the defined ranges into a single vector for later use in indexing.
############# read sperms databases

wd<-"C:/..."
setwd(wd)
l<-list.files(pattern = "db.RData")
##Loop through each file, construct the full file path, and load the contraction and fluorescence files.
#Use readSingleContVSFnorm to extract and process data from each file.
#Combine the results from each file into d.n and f.n.
d.n<-f.n<-NULL
for(f in 1:length(l)){
 db<-paste(wd,l[f], sep = "/")
 load(db)
 DB<-readSingleContVSFnorm(db)
 d.n<-c(d.n, DB$dn)
 f.n<-c(f.n, DB$fn)
}

#Filter the combined data to keep only the values within a specific range for further analysis.

id1<-which(d.n>=0.5 & d.n<=1.15)    #which(d.n >= 0.5 & d.n <= 1.15) selects indices where d.n falls between 0.5 and 1.15. This range is chosen to focus on a relevant subset of the data.
f.n<-f.n[id1]
d.n<-d.n[id1]     # Apply these indices to both d.n and f.n to filter the data.

#Compute mean and standard error of f.n for bins of d.n values.
f.n.m<-NULL
f.n.se<-NULL      #f.n.m for mean values and f.n.se for standard errors.
xseq<-seq(0.45,1.1, by=0.025);xseq    # xseq creates bins for d.n ranging from 0.45 to 1.1 micrometers, with a step size of 0.025 micrometers.
for(i in xseq){
  id<-which(d.n>=i & d.n<i+0.025)     #For each bin defined by xseq, find indices where d.n falls within the bin range.
  f.n.m <-c(f.n.m, mean(f.n[id],na.rm=T))
  f.n.se <- c(f.n.se, sd(f.n[id],na.rm=T))    #Calculate the mean (f.n.m) and standard deviation (f.n.se) for these bins.
  print(i)    #Print the current bin value.
  }

#Perform linear regression, plot the results, and visualize the relationship between the data.

plot(f.n.m)   #plot(f.n.m) creates a simple plot of the mean values. This is an initial step for visualization.

upL<-f.n.m+f.n.se
lowL<-f.n.m-f.n.se    #upL and lowL define the upper and lower limits of the confidence intervals using the standard error.

df = data.frame(cbind(upL,lowL,f.n.m))    #Combine these into a dataframe df for plotting.


simple.fit<- lm(f.n~d.n)  #lm(f.n ~ d.n) fits a linear model where f.n is the response variable and d.n is the predictor.
summary(simple.fit)   
b<-simple.fit$coefficients[1]
m<-simple.fit$coefficients[2]
expression(y=mx+b)

#Create a scatter plot of f.n vs. d.n.
plot(f.n~d.n, ylab="", xlab="",
     col="grey",cex.main= 2,
     ylim=c(0.1,8),pch=19, cex=0.3, xlim=c(0.5,1.1))
plotCI(xseq, y=df$f.n.m, uiw=df$upL-df$f.n.m, liw=df$f.n.m-df$lowL, err="y",
       pch=20, cex=2, slty=1, lwd=2, scol = "black", add=T)   #Add error bars using plotCI with the confidence intervals from df.
abline(h=0.1)
points(predictedcounts, f.n)
#Overlay the linear regression line with abline(simple.fit).
abline(simple.fit)
Counts=f.n
Time=d.n

linear.model <-lm(Counts ~ Time)
summary(linear.model)

