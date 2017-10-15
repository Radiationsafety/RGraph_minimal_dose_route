# Hamiltonian path, route for visiting given points with a minimal dose, 8 directions (including diagonal), mask of roads
ptm <- proc.time() # duration of calulation

library(raster); library(rasterVis);library(rgdal)
library(igraph);library(geosphere);library(plotrix)     #  igraph for work with graphs

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
mask_buildings<-raster("170412good_mask_buildings.tif")
extent(mask_buildings)<-wsp
mask_buildings<-resample(mask_buildings, grid_graph, method='ngb') # the mask in same resolutioan as a grid
nrow(mask_buildings); ncol(mask_buildings)
# Criterion of the action of the mask 1000 - visible, 0- invisible, 150 - color boundary (0-255), giving boundaries of buildings
mask_buildings[which(mask_buildings[]>0.60)]=1000  
mask_buildings[which(mask_buildings[]<0.60)]=0
#plot(mask_buildings,alpha=0.7,add=TRUE,gridded=T) 
# end of road mask

#Selection of the start and end points on the radiation situation map
numberofpoints = 4	# number of points to visit								
par(plt=c(0,1,0,1))
p7<-click(grid_graph,n=numberofpoints,cell=TRUE,type="p",plt=c(0,1,0,1))# (n= number of points to visit)
point_1<-p7[1,1];point_2<-p7[1,1]   # The route starts and ends at the first point
#Selection end

#Creation of a date-frame describing the weight of the edges of the graph-lattice
i=1
df=data.frame(x=double(),y=double(),dose_rate=double()) # empty dataframe for weights
df_dist=data.frame(x=double(),y=double(),dose_rate=double()) # empty dataframe for distance weights
while(i < ncell(grid_graph)) {
  if (mask_buildings[i]!=0){ # If there is a mask, then do not build an egde 
    if ((ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol))&(mask_buildings[i+1]!=0)){          # Check for the last column
      df=rbind(df,c(i,i+1,wsp3[i,2]))  
      df_dist=rbind(df_dist,c(i,i+1,1)) 
    }
    if (ceiling(i/grid_grapn_ncol)<grid_grapn_nrow){  # Check for the last row (строка)
      if (mask_buildings[i+grid_grapn_ncol]!=0){
        df=rbind(df,c(i,i+grid_grapn_ncol,wsp3[i,2])) 
        df_dist=rbind(df_dist,c(i,i+grid_grapn_ncol,1))
      }
      if (ceiling(i/grid_grapn_ncol)!=(i/grid_grapn_ncol)){          #1 diagonal "\", Multiply by 1.4 - Longer walk diagonally
        if (mask_buildings[i+1+grid_grapn_ncol]!=0){
          df=rbind(df,c(i,i+1+grid_grapn_ncol,wsp3[i,2]*1.4))  
          df_dist=rbind(df_dist,c(i,i+1+grid_grapn_ncol,1.4))   
        }
      }
      if ((i!=1)&((i-1)/grid_grapn_ncol!=floor((i-1)/grid_grapn_ncol))){  #2 diagonal "/", Multiply by 1.4 - Longer walk diagonally
        if (mask_buildings[i-1+grid_grapn_ncol]!=0){
          df=rbind(df,c(i,i-1+grid_grapn_ncol,wsp3[i,2]*1.4))  
          df_dist=rbind(df_dist,c(i,i-1+grid_grapn_ncol,1.4))   
        }
      }
      
    }
    
  }
  i=i+1
}
#end of creation of a date-frame

el_dist<-as.matrix(df_dist)
g_dist<- add.edges(graph.empty(ncell(grid_graph)), t(el_dist[,1:2]), weight=el_dist[,3])
g_dist<-as.undirected(g_dist)
#end of creation of a distance weight date-frame

#Computation of the shortest path in the graph
el<-as.matrix(df)
g<- add.edges(graph.empty(ncell(grid_graph)), t(el[,1:2]), weight=el[,3])
g<-as.undirected(g)

# optimal route finder
n_points<-length(p7[,1])    # number of points to visit
dose=matrix(nrow = n_points, ncol = n_points)  #dose by summarizing the edge weights (More precisely, takes into account large doses on the diagonal 1.4)
#dose2=matrix(nrow = n_points, ncol = n_points) #dose by extracting from grid (not right)
dose_dist=matrix(nrow = n_points, ncol = n_points) #storage for distance per route
DF_XY=data.frame(i=integer(),j=integer(),rout=integer()) #storage for routes
for(i in 1:(n_points)){
  for(j in i:n_points){
    if (i!=j) {
      routine<-get.shortest.paths(g,from=p7[i,1],to=p7[j,1],weights = NULL, output = "both") # search for the route with a minimal dose
      List_epath=unlist(routine$epath)
      dose[i,j]<-sum(E(g)[List_epath]$weight)
      dose_dist[i,j]<-sum(E(g_dist)[List_epath]$weight)
      rout<-unlist(routine$vpath)
      DF_XY=rbind(DF_XY,data.frame(p7[i,1],p7[j,1],rout))
    } else {
      dose[i,j]<-0
      dose_dist[i,j]<-0
    }
    print(paste(i,j,dose[i,j],dose_dist[i,j]))
    print(paste(p7[i,1],p7[j,1]))
    print(rout)
    print("------")
    #Calculation of the dose by the shortest route
    #dd<-grid_graph[rout]
    #print (paste("dd=",dd,"sum(dd)=",sum(dd), '| sum(E(g))=',dose[i,j]))
    # obtain the values of dose rate for the optimal trajectory
    #dose2[i,j]<-sum(dd)                       # find the dose along the optimal route by extracting from grid
  }
}
#The calculation of the dose matrix for the shortest routes between the visiting points is completed

#Calculation of doses for all possible Hamiltonian paths between points of visits (the first point is the beginning and the end of the path is always)
dose[lower.tri(dose)] = t(dose)[lower.tri(dose)] # The dose matrix is symmetrized (!!!2017-07-27)
Dose=dose
dose_dist[lower.tri(dose_dist)] = t(dose_dist)[lower.tri(dose_dist)] # The distance matrix is symmetrized
g1 <- graph.ring(n_points)   # Random permutation of a random graph g1
kk<-factorial(n_points)     
D_h<-rep(0,kk)            # A vector is created to store the doses
D_h_dist<-rep(0,kk)       # A vector is created to store the distance on every step in the grid
rr<-matrix(nrow=kk,ncol=n_points,c(0:0)) # A matrix is created for the saving of the order of bypass of selected Points
ss<-matrix(nrow=kk,ncol=n_points,c(0:0)) # A matrix is created for storing the "sample"
for (j in 1:kk){ #A cycle of calculation of doses is started on randomly chosen Hamiltonian paths
  smp<-sample(vcount(g1))
  g2 <- permute.vertices(g1,smp)#A random ring graph is created
  Dosemat<-get.edges(g2,es=1:n_points)       #Created a dataframe with the order of traversing vertices
  el2<-as.matrix(Dosemat)
  dl<-rep(0,n_points)
  el2<-cbind(el2,dl) #The creation of a matrix of weights of one of the possible Hamiltonian graphs
  el2_dist<-el2 #for total time on route calculation
  for (i in 1:n_points) {
    el2[i,3]<-Dose[el2[i,1],el2[i,2]]
    el2_dist[i,3]<-dose_dist[el2[i,1],el2[i,2]]
  }
  D<-colSums(el2)   #The dose value is calculated from the current Hamiltonian cycle
  D_h[j]<-D[3]
  D_h_dist[j]<-colSums(el2_dist)[3] #sum of steps on the route												  
  ss[j,]<-smp
  rr[j,]<-el2[,1] #The dose along the Hamiltonian route and the route itself are stored
}
max(D_h);min(D_h)# Determined min and max dose over all Hamiltonian bypass route
NN<-which(D_h==min(D_h))[1] # first variant with minimum dose
total_dist=D_h_dist[NN]
#Construction of the route of detour of selected visiting points
r_beg<-ss[NN[1],]
r_end<-rep(0,n_points)
linecolors= smoothColors("#FFD000",n_points,"#FF0000") # colors for every i-route
for (i in 1:(n_points-1)){r_end[i]<-r_beg[i+1]}
r_end[n_points]<-r_beg[1]
pairs<-cbind(r_beg,r_end)#A matrix of transitions between the vertices of the Hamiltonian graph with min dose

for(i in 1:n_points){
  print(paste(p7[pairs[i,1],1],"->",p7[pairs[i,2],1]))
  pointS1=p7[pairs[i,1],1]
  pointS2=p7[pairs[i,2],1]
  rout=DF_XY[which((DF_XY[,2]==pointS2) &(DF_XY[,1]==pointS1)),3] 
  # if it is no 1-2, then we mirror it for 2-1.
  if (length(rout)==0) {
    rout=DF_XY[which((DF_XY[,2]==pointS1) &(DF_XY[,1]==pointS2)),3]
  }
  XY<-xyFromCell(grid_graph,rout)        #Obtaining coordinates of waypoints from raster cell numbers
  par(plt=c(0,1,0,1))
  lines(XY,col=linecolors[i],lwd=3)  
}

#dose on route
# The dose in the weigths was in mcSv/h, reworked for the size of the map, calculate the size of the map in meters
xdist=distm (c(ymin(wsp), xmin(wsp)), c(ymax(wsp), xmin(wsp)), fun = distHaversine)
ydist=distm (c(ymin(wsp), xmin(wsp)), c(ymin(wsp), xmax(wsp)), fun = distHaversine)
stepx=xdist/Ncol
#ydist=ymax(wsp) -ymin(wsp) # Step in the horizontal and vertical are equal (specified in interpolation.R)
#stepy=ydist/Nrow
# total_dist means sum of all grid cells on the route, where 1 = dist for passing horizontal and 
#vertical cells, 1.4 - for diagonal.
# if the distance within one cell is "stepx" value (for square grid) and walkspeed is
walkspeed=1
# then the distance of transition throw all cell (in meters) =total_dist*stepx
Total_dist_sized = total_dist*stepx
# then the time for all route is [in seconds]
Total_time=Total_dist_sized/walkspeed
# Dose is =
walktimestep=(stepx/walkspeed)/3600 # in hours
Total_Dose=min(D_h) * walktimestep
#Time on the route is conditional, because the diagonal elements are also considered to pass for the "walktimestep"
paste('The dose on the route after visiting all the given points =',round(Total_Dose,2),'mcSv', 'route length=',round(Total_dist_sized,2),'m','time on route=',round(Total_time/60,2),'minutes, speed=',walkspeed,'m/s') 
# Time of calculation
proc.time() - ptm
