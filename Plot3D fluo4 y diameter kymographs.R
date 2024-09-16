#This R script uses the plotly package to create two 3D surface plots based on fluorescence and 
#diameter data from a dataset. It sets up the working directory, loads the data, and then subsets 
#it for plotting. The plots visualize how fluorescence and diameter change across the length of a
#flagellum over time. The layout and add_trace functions in plotly are used to configure and display these 3D plots.


library(plotly) #Loads the plotly package, which is used for interactive data visualization. 
                #plotly is particularly good for creating 3D plots and other complex visualizations.

wd<-"D:/.."
setwd(wd)
load("fluorescence kymograph data for Fluo4 and contraction kymograph")

############ parameters
magnification<-0.117 #um
time<- (3:nrow("fluorescence kym"))* "acquisition time"   #the frames after the induction of acrosome reaction are analyzed
ys<-(1:ncol("fluorescence kym")*0.117)    #Creates a vector of length values by multiplying column indices 
                                          #(from 1 to the number of columns in the fluorescence and contraction kymographs) by 0.117



#####fluo4
fluoacot<-"fluo4 kymographs"[3:20,]   #the frames after the induction of acrosome reaction are analyzedSets up the 3D scene with axis titles.

f<- plot_ly()%>%    #Initializes a new plotly plot object.
     layout(scene = list(        # Sets up the 3D scene with axis titles.
    xaxis = list(title = "  Length across <br> flagellum (µm)"),
    yaxis = list(title = "Time (min)"),
    zaxis = list(title = "Fluo4 normalized <br>    fluorescence")
  ))%>% 
  add_trace(z = fluoacot,x=ys,  y=time, type = "surface")%>%  #Adds a surface plot layer using fluoacot data for the z values, ys for x, and time for y.
  layout(scene=list(aspectmode='cube'))   # Ensures the aspect ratio of the plot is equal, so the plot's dimensions are scaled equally. 
f    #Displays the plot.


####diameter     Similar to the fluorescence plot, but this one visualizes the normalized diameter data.
diamacot<-"contraction kymograph"[3:20,]    #the frames after the induction of acrosome reaction are analyzed

d <- plot_ly() %>% 
  layout(scene = list(
    xaxis = list(title = "  Length across <br> flagellum (µm)"),
    yaxis = list(title = "Time (min)"),
    zaxis = list(title = "Normalized <br>  diameter")
  ))%>% 
  add_trace(z = diamacot, x=ys, y=time, type = "surface",
            colorscale = list(c(0, 1),   c("yellow", "#B70F0F")))%>%    #Specifies a color gradient for the surface plot from yellow to a dark red color, based on the values in diamacot.
  layout(scene=list(aspectmode='cube'))

d





