#An algorithm for finding an optimal Chinese postman route is:
#Step 1 List all odd vertices.
#Step 2 List all possible pairings of odd vertices.
#Step 3 For each pairing find the edges that connect the vertices with the minimum weight.
#Step 4 Find the pairings such that the sum of the weights is minimised.
#Step 5 On the original graph add the edges that have been found in Step 4.
#Step 6 The length of an optimal Chinese postman route is the sum of all the edges added to the total found in Step 4.
#Step 7 A route corresponding to this minimum weight can then be easily found.

ptm <- proc.time() # duration of calulation

library(igraph)
#Creating an example of the graph of the transport and pedestrian network of the facility
# random values for test example see the figure "Scheme for Chinese postman problem.jpg"
# edge weight, the dose is in mcSv
el <- matrix(nc=3, byrow=TRUE,
             c(1,2,0.1,
               1,21,0.1,
               1,2,0.1,
               1,3,0.2,
               1,19,0.3,
               2,3,0.1,
               2,4,0.1,
               3,4,0.1,
               3,6,0.1,
               4,6,0.3,
               4,5,0.1,
               5,6,0.6,
               5,7,0.2,
               5,7,1.1,
               7,8,0.1,
               7,10,0.3,
               7,17,0.4,
               8,9,0.2,
               8,22,1.1,
               8,23,0.1,
               9,10,0.4,
               9,15,0.1,
               10,11,1.2,
               11,12,3.3,
               11,13,2.1,
               13,14,0.9,
               13,15,1.2,
               15,16,1.0,
               17,18,0.2,
               17,19,0.2,
               19,20,0.8))
g<- add.edges(graph.empty(max(el)), t(el[,1:2]), weight=el[,3])
g<-as.undirected(g)
plot(g,edge.label=format(E(g)$weight,digits=1))
d<-degree(g)
g1<-g #double of parameter for check

first_v<-1 # The choice of the vertex with which the graph starts to traverse

#the construction of the graph for the vertices of odd degree with minimum weights
di<-which((d/2)>floor(d/2))
sp<-shortest.paths(g1)
nr=length(di) # Number of vertices in a new full graph
ne=(nr*(nr-1))/2 # Number of edges in a new full graph
# creation of a matrix of weights
el1<-matrix(nc=3, nr=nr*nr,rep(0,3*nr*nr))
for (i in 1:nr){
	for (j in 1:nr){
		el1[(i-1)*nr+j,1]<-di[i]
		el1[(i-1)*nr+j,2]<-di[j]
		el1[(i-1)*nr+j,3]<-sp[di[i],di[j]]
	}
}
el2<-el1[el1[,3]>0,]# Remove loops
el2<-el2[el2[,1]<el2[,2],]#Remove repeated weights (2-4 | 4-2)
# The matrix of weights for the new graph is ready
g2<-graph.empty(vcount(g1))#Start creating a new graph
g3<- add.edges(g2, t(el2[,1:2]), weight=el2[,3])
g4<-as.undirected(g3)
########################################################################
#Choose the minimal covering of the graph
np<-nr/2 #Determine the required number of pairs in the graph (the depth of recursion). The number of vertices of nr is always even.
ell<-as.data.frame(el2)
wsp_el<-matrix(nc=ne, nr=ne,rep(0,ne*ne))
#Elements of the matrix wsp_el show which matchings have no common boundaries
#If the value of the element is 10, then this pair has a common vertex,
#Otherwise, this is the total weight of the two edges
for (i in 1:ne){for (j in 1:ne){
	a=ell[i,1];b=ell[i,2];c=ell[j,1];d=ell[j,2];
	wsp_el[i,j]<-!((a==c)|(a==d)|(b==c)|(b==d))
	if(wsp_el[i,j]>0){wsp_el[i,j]<-ell[i,3]+ell[j,3]}
	else wsp_el[i,j]<-10
	}}
ni<-1
ll<-0
repeat{
	l<-which(wsp_el==min(wsp_el),arr.ind = TRUE)
	ll<-c(ll,l[1,])
	wsp_el[,wsp_el[l[1,1],]==10]<-10;wsp_el[,wsp_el[l[1,2],]==10]<-10 # index: two commas
	ni<-ni+1
	if(ni>(np-1))break
}
lll<-unique(ll)
l4<-lll[2:(np+1)]   
l4<-sort(l4) # The numbers of the edges that form the minimal covering of the graph g4
ell_min<-ell[l4,] #Cut out its minimal covering from the graph g4
el3<-as.matrix(ell_min)

#construction of a graph of the minimal covering
g5<-graph.empty(vcount(g1))
g6<- add.edges(g5, t(el3[,1:2]), weight=el3[,3])
g7<-as.undirected(g6)
plot(g7,edge.label=E(g7)$weight) #A minimal covering graph is constructed

#Now you need to draw new (additional) edges of the graph correctly
#Explanation: for example, in fact, the edge "3-8" is not exist, but there is a route on edges "3-2-6-8", which must be restored!
pl_edges<-c(0,0,0)
pl_edges<-as.data.frame(t(pl_edges))
for (i in 1:length(el3[,1])){
	el_list<-get.shortest.paths(g1,from=el3[i,1],to=el3[i,2])
	el_list<-unlist(el_list,use.names = FALSE)
	leng<-length(el_list)
	aa<-bb<-cc<-rep(0,(leng-1))
	for (j in 1:(leng)-1){
		aa[j]<-el_list[j];bb[j]<-el_list[j+1]
		dd<-cbind(aa[j],bb[j],cc[j])
		pl_edges<-rbind(pl_edges,dd)
	}
}
pl_edges<-pl_edges[2:length(pl_edges[,1]),]# dataframe with additional edges

# Present the primary graph as a dataframe
graph.table<-cbind(as.data.frame(get.edgelist(g1)),E(g1)$weight)
el<-graph.table

#restore weights for new edges
for (i in 1:length(pl_edges[,1])){
	aaa<-pl_edges[i,1]
	bbb<-pl_edges[i,2]
	wh<-which((el[,1]==aaa)&(el[,2]==bbb)|(el[,1]==bbb)&(el[,2]==aaa))
	pl_edges[i,3]<-el[wh,3]
}

#build a graph of the pedestrian transport network  with additional edges
g8<- add.edges(g1, t(pl_edges[,1:2]), weight=pl_edges[,3])
g9<-as.undirected(g8)
plot(g9,edge.label=format(E(g9)$weight,digits=1))
#The graph g9 became an Euler graph - all vertices have an even number of edges
degree(g9) #Check whether all vertices of g9 have even degrees
	#Now is needed to find the Euler cycle for the selected vertex.
	#In the absence of vertices of odd degree, there is an Euler cycle, 
	#the algorithm of construction of which allows to leave from any vertex,
	#sequentially remove the traversed edges and do not walk along the edge,
	#which at the moment is an isthmus, after its removal the graph
	#splits into several connectivity components.
el4<-get.edgelist(g9) #conversion
lll<-ecount(g9)
el5<-matrix(nrow=1,ncol=3,c(0,0,0)) #An empty matrix is created where the directed graph of the traversal of the transport network will be stored

g10<-g9 #double the parameter for check
#Beginning of the cycle search cycle of the graph traversal
for (j in (1:lll)){
  el4<-get.edgelist(g10)
  n<-which((el4[,1]==first_v)|(el4[,2]==first_v))
  v2_1<-el4[n,1]
  v2_2<-el4[n,2]
  v2<-c(v2_1,v2_2)
  v2<-v2[!(v2==first_v)]
  i=0
  repeat{
  	i=i+1
	  vr<-v2[i]
	  m<-edge.connectivity(g10,first_v,v2[i])
		if( (m>1)|((m=1)&(i==length(v2))))break
	}
  sequence<-c(first_v,vr,j)
  g10<-g10-path(first_v,vr) #Removing an edge that is not an isthmus
  el5<-rbind(el5,sequence)
  first_v<-vr
} 
el5<-el5[2:length(el5[,1]),] #Remove the zero initial line
g11<- add.edges(graph.empty(vcount(g1)), t(el5[,1:2]),weight=el5[,3]) # graph is built of the bypass of the roads of the facility
autocurve.edges (g11,0.5) # Curving of some edges for better visibility
plot(g11,edge.label=E(g11)$weight,edge.label.cex=0.9,edge.arrow.size=0.7,vertex.size=7)
Dose=sum(E(g11)$weight)
paste('Dose on route=',round(Dose,4),'mcSv')

# Time of calculation
proc.time() - ptm
