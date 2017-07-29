# Route with minimal dose from point A to B, 8 directions (including diagonal), mask of roads
ptm <- proc.time() # duration of calulation

library(raster); library(rasterVis);library(rgdal)
library(igraph);library(geosphere) ; library(plotrix) #     igraph for work with graphs

basemap<-raster("industrial site.gif") 
grid1<-read.table("GRidded_data_RMO_OTZPVH_test.dat")   # Read grid from interpolation.R script
wsp<-extent(basemap)       
Nrow=length(unique(grid1[,2]))
Ncol=length(unique(grid1[,1]))

grid_wsp<-raster(wsp,nrow=Nrow,ncol=Ncol)    
grid_wsp[]<-abs(grid1[,3])                  # Making raster for grid                 
grid<-flip(grid_wsp,2)                      # flip - first cell is bottom left
plot(basemap,axes=T) 
grid_graph<-grid         

plot(grid_graph,alpha=0.3,add=TRUE,col=rainbow(6),gridded=T) # plot a grid on map

# raster cells = grid cells
wsp2<-as.data.frame(grid_graph)            # map as a vector of values
n1<-1:ncell(grid_graph)
grid_grapn_ncol= ncol(grid_graph)
grid_grapn_nrow= nrow(grid_graph)
wsp3<-cbind(n1,wsp2)                       # dataframe with numbers of cells and dose rate 

# road mask
mask_buildings<-raster("170412good_mask_roads_inv.tif")
extent(mask_buildings)<-wsp
mask_buildings<-resample(mask_buildings, grid_graph, method='bilinear') # the mask in same resolutioan as a grid
nrow(mask_buildings); ncol(mask_buildings)
# Criterion of the action of the mask 1000 - visible, 0- invisible, 150 - color boundary (0-255), giving boundaries of buildings
mask_buildings[which(mask_buildings[]>60)]=1000  
mask_buildings[which(mask_buildings[]<60)]=0
#plot(mask_buildings,alpha=0.7,add=TRUE) 
# end of road mask

#Selection of the start and end points on the radiation situation map
par(plt=c(0,1,0,1))
p7<-click(grid_graph,n=2,cell=TRUE,type="p",plt=c(0,1,0,1))# Select the start and end points on the map
point_1<-p7[1,1];point_2<-p7[2,1]
#Selection end

#Creation of a date-frame describing the weight of the edges of the graph-lattice
i=1
df=data.frame(x=double(),y=double(),dose_rate=double()) # empty dataframe for weights
while(i < ncell(grid_graph)) {
  if (mask_buildings[i]!=0){ # If there is a mask, then do not build an egde 
    if (ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol)){          # Check for the last column
      df=rbind(df,c(i,i+1,wsp3[i,2]))  
    }
    if (ceiling(i/grid_grapn_nrow)<grid_grapn_nrow){          # Check for the last row
      df=rbind(df,c(i,i+grid_grapn_ncol,wsp3[i,2])) 
      if (ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol)){    #1 diagonal "\", Multiply by 1.4 - Longer walk diagonally
        df=rbind(df,c(i,i+1+grid_grapn_ncol,wsp3[i,2]*1.4))  
      }
      if (i>1) { 
        if (ceiling((i-1)/grid_grapn_ncol)!=((i-1)/grid_grapn_ncol)){    #2 diagonal "/", Multiply by 1.4 - Longer walk diagonally
          df=rbind(df,c(i,i-1+grid_grapn_ncol,wsp3[i,2]*1.4))  
        }
      }
    }
    
  }
  
  i=i+1
}
#end of creation of a date-frame

#Creation of a date-frame describing the __distance__weight__ of the edges of the graph-lattice
# time for 1 horizontal edge= 1, for diagonal =1.4
i=1
df_dist=data.frame(x=double(),y=double(),dose_rate=double()) # empty dataframe for weights
while(i < ncell(grid_graph)) {
  if (mask_buildings[i]!=0){ # If there is a mask, then do not build an egde 
    if (ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol)){          # Check for the last column
      df_dist=rbind(df_dist,c(i,i+1,1))  
    }
    if (ceiling(i/grid_grapn_nrow)<grid_grapn_nrow){          # Check for the last row
      df_dist=rbind(df_dist,c(i,i+grid_grapn_ncol,1)) 
      if (ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol)){    #1 diagonal "\", Multiply by 1.4 - Longer walk diagonally
        df_dist=rbind(df_dist,c(i,i+1+grid_grapn_ncol,1.4))  
      }
      if (i>1) { 
        if (ceiling((i-1)/grid_grapn_ncol)!=((i-1)/grid_grapn_ncol)){    #2 diagonal "/", Multiply by 1.4 - Longer walk diagonally
          df_dist=rbind(df_dist,c(i,i-1+grid_grapn_ncol,1.4))  
        }
      }
    }
    
  }
  
  i=i+1
}
el_dist<-as.matrix(df_dist)
g_dist<- add.edges(graph.empty(ncell(grid_graph)), t(el_dist[,1:2]), weight=el_dist[,3])
g_dist<-as.undirected(g_dist)
#end of creation of a distance weight date-frame
#Computation of the shortest path in the graph
el<-as.matrix(df)
g<- add.edges(graph.empty(ncell(grid_graph)), t(el[,1:2]), weight=el[,3])
g<-as.undirected(g)
#plot(g,edge.label=round(E(g)$weight,2)) # plot a graph

#dose by summarizing the edge weights (More precisely, takes into account large doses on the diagonal 1.4)
#dose_dist storage for distance per route
#DF_XY=data.frame(i=integer(),j=integer(),rout=integer()) 
routine<-get.shortest.paths(g,from=point_1,to=point_2,weights = NULL, output = "both") # Find the path with the lowest dose rate
List_epath=unlist(routine$epath)
#Calculation of the dose by the shortest route  
dose<-sum(E(g)[List_epath]$weight)
dose_dist<-sum(E(g_dist)[List_epath]$weight)
rout<-unlist(routine$vpath)
XY<-xyFromCell(grid_graph,rout)        #Obtaining coordinates of waypoints from raster cell numbers
par(plt=c(0,1,0,1))
lines(XY,col='blue',lwd=2)              #Image on the route map with the lowest #dose                               #dose on route

# The dose in the weigths was in mcSv/h, reworked for the size of the map, calculate the size of the map in meters
xdist=distm (c(ymin(wsp), xmin(wsp)), c(ymax(wsp), xmin(wsp)), fun = distHaversine)
ydist=distm (c(ymin(wsp), xmin(wsp)), c(ymin(wsp), xmax(wsp)), fun = distHaversine)
stepx=xdist/Ncol
stepy=ydist/Nrow
# total_dist means sum of all grid cells on the route, where 1 = dist for passing horizontal and 
#vertical cells, 1.4 - for diagonal.
# if the distance within one cell is "stepx" value (for square grid) and walkspeed is
walkspeed=1
# then the distance of transition throw all cell (in meters) =total_dist*stepx
Total_dist_sized = dose_dist*stepx
# then the time for all route is [in seconds]
Total_time=Total_dist_sized/walkspeed
# Dose is =
walktimestep=(stepx/walkspeed)/3600 # in hours
Total_Dose=dose * walktimestep
#Time on the route is conditional, because the diagonal elements are also considered to pass for the "walktimestep"
paste('The dose on the route from the first to the second point is D=',round(Total_Dose,2),'mcSv', 'route length l=',round(Total_dist_sized,2),'m','time on route t=',round(Total_time/60,2),'minutes, speed v=',walkspeed,'m/s') 
# Time of calculation
proc.time() - ptm
