# This file contains functions to capture the characteristics of a network.
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  
#### This file is to obtain the statistics of networks #

####Load the following packages first, and set the working directory####
install.packages("igraph")
library(igraph)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


#### obtain property/characteristics from EconLit coauthor network ####
#This function includes functions listed: transitive_triples_clustering(), average_clustering(), get.giant.component() #

get_network_characteristics <- function(gh){

  #### A function to get giant component####
  get.giant.component <- function(g) {
    
    if (!is.igraph(g)) stop("g is not an igraph object")
    
    comps <- clusters(g, mode="weak")
    induced.subgraph(g, which(which.max(comps$csize)==comps$membership))
  }
  
  
  #### The function to calculate average clustering function #######
  # what we do is replace NaN with 0. E.g. transitivity(gh_small,"local")     NaN 0.1666667 1.0000000 1.0000000  NaN.
  # The right calculation is (  0+ 0.1666667+ 1.0000000+ 1.0000000 + 0)/5, according to Newman(2003) page 11
  average_clustering <- function(gh) {
    array_local_clustering <- transitivity(simplify(gh),"local") # simplify the graph, ignore the multiple edges. This is an array with NaN.
    array_local_clustering[is.nan(array_local_clustering)] <- 0
    mean(array_local_clustering)
  }
  
  
  #### The function to calculate the clustering coefficient of transitive triples CTT function############################
  transitive_triples_clustering <- function(gh){
    number_triad <- triad_census(gh)
    number_triad
    number_triplets_triples_para_transitivity <- array(c(0,0,0,0,0,1,1,1,1,3,2,2,2,3,4,6,0,0,0,0,0,0,0,0,1,0,0,2,2,1,3,6), c(16,2))
    Clustering_transitivity_triplets_triples <- apply(number_triad*number_triplets_triples_para_transitivity, c(2),sum) 
    Clustering_transitivity_triplets_triples
    Clustering_transitivity_simu <- Clustering_transitivity_triplets_triples[2]/Clustering_transitivity_triplets_triples[1] #  transitivity clustering by simulation this should be correct  
    Clustering_transitivity_simu
  } 
  
  
  a=length(V(gh)) # 9877 nodes; Econlit 129003 for the 1970-1999 overall network
  cat("Number of nodes:", a,"\n")
  b=gsize(gh) #number of edges, 51971 edges; Econlit 106705
  cat("number of edges (not articles):" , b,"\n")
  c=mean(degree(gh, v=V(gh),mode="all")) # 1.511972 average degree; Econlit 1.654303
  cat("average overall degree:",c,"\n" )
  d=mean(degree(gh, v=V(gh),mode="out")) # average outdegree
  cat("average outdegree:",d ,"\n")
  e=length(V(get.giant.component(gh))) # the number of nodes of the giant component
  cat("The number of nodes of the giant component",e,"\n")
  f=max( clusters(gh)$csize[clusters(gh)$csize!=max(clusters(gh)$csize)] ) #number of nodes of the second largest component
  cat("The number of nodes of the second largest component:",f,"\n")
  g=length(which(clusters(gh)$csize ==1 )) # count isolated nodes
  cat("The number of isolated nodes:",g,"\n")
  var_overall_m = var(degree(gh,mode = "all")) # The variance of overal degree
  cat("The variance of overal degree:",var_overall_m,"\n")
  var_m = var(degree(gh,mode = "out")) #variance of out-degree
  cat("The variance of outdegree:",var_m,"\n")
  od=max(degree(gh,mode = "all")) # get the max value of overall degree
  cat("The max of overal degree:", od,"\n" )
  j=min(degree(gh,mode = "all")) # get the min value of overall degree
  cat("The min of overal degree:",j,"\n")
  h=max(degree(gh,mode = "out")) # get the max value of out degree
  cat("The max of outdegree:",h,"\n")
  k=min(degree(gh,mode = "out")) # get the min value of out degree
  cat("The min of outdegree:",k,"\n")
  i=max(degree(gh,mode = "in")) 
  cat("The max of indegree:",i,"\n")
  l=min(degree(gh,mode = "in")) 
  cat("The min of indegree:",l,"\n")
  #m=length(E(gh)[which(head_of(gh, E(gh))==tail_of(gh, E(gh)))]) #get the number of authors who both work independently and cooperately. The links are already unique. That is, it is impossible to have two loops on the same node.
  #cat("The number of authors who both work independently and cooperately:",m,"\n")
  
  p=transitive_triples_clustering(gh) # get clustering of transitive triples CTT
  cat("clustering of transitive triples:",p,"\n")
  o=transitivity(gh, type = c("global")) # get the total clustering
  cat("total clustering:",o,"\n")
  n=average_clustering(gh) # get average clustering
  cat("average clustering:",n,"\n")
  
  cat("E[m^2]",mean((degree(gh,mode = "out"))^2),"\n")
  cat("E[m^3]",mean((degree(gh,mode = "out"))^3),"\n")
  
  q=diameter(gh) # get the diameter
  cat("diameter:",q,"\n")
  r=assortativity.degree(gh)
  cat("degree assortativity:",r,"\n")
  s=centralization.degree(gh)$centralization #0.006049436 , theoretical_max: 195070752
  cat("centralization.degree",s,"\n")
  
  
}



