# This file is to derive all simplified networks from the original network.
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  

####Load the following packages first, and set the working directory####
install.packages("igraph")
install.packages("foreign")
library(igraph)
library(foreign)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

#### The networks of 1970-1999  ####
gh_kpmul_kploop <- read_graph("digraph_kpmul_kploop.csv") # load the overall network of 1970-1999
gh_kpmul_rmloop <- igraph::simplify(gh_kpmul_kploop,remove.multiple = FALSE, remove.loops = TRUE ) # We create the network where multiple links are kept, while loops are removed.
gh_rmmul_kploop <- igraph::simplify(gh_kpmul_kploop,remove.multiple = TRUE, remove.loops = FALSE ) # We create the network where multiple links are removed, while loops are kept.
gh_rmmul_rmloop <- igraph::simplify(gh_kpmul_kploop,remove.multiple = TRUE, remove.loops = TRUE ) # We create the network where both the multiple links and loops are removed.
# Next we remove all the isolated nodes out of the network of gh_kpmul_kploop 
isolated_nodes_from_kpmulkploop <- V(gh_kpmul_kploop)[which(degree(gh_rmmul_rmloop,mode="out")==0 & degree(gh_rmmul_rmloop,mode="in")==0)]
gh_kpmul_kploop_rmisonodes <- gh_kpmul_kploop - isolated_nodes_from_kpmulkploop
gh_kpmul_rmloop_rmisonodes <- igraph::simplify(gh_kpmul_kploop_rmisonodes, remove.multiple = FALSE, remove.loops = TRUE) # We remove all the isolated nodes out of the network of gh_kpmul_rmloop 
gh_rmmul_kploop_rmisonodes <- igraph::simplify(gh_kpmul_kploop_rmisonodes, remove.multiple = TRUE, remove.loops = FALSE) # We remove all the isolated nodes out of the network of gh_rmmul_kploop 
gh_rmmul_rmloop_rmisonodes <- igraph::simplify(gh_kpmul_kploop_rmisonodes, remove.multiple = TRUE, remove.loops = TRUE) # We remove all the isolated nodes out of the network of gh_rmmul_rmloop 

# generate network for EB model, undirected network - isolated 
gh = get_undirected_network(1970,1999)
isolated_nodes_from_gh <- V(gh)[which(degree(gh,mode="out")==0 & degree(gh,mode="in")==0)]
gh_rmisonodes <- gh - isolated_nodes_from_gh
gh = gh_rmisonodes

#### The networks of 1970-1989  ####
gh_1970_1989_kpmul_kploop <- read_graph("digraph_1970_1989_kpmul_kploop_201909111502_graphml.txt", format = "graphml")
gh_1970_1989_kpmul_rmloop <- igraph::simplify(gh_1970_1989_kpmul_kploop, remove.multiple = FALSE, remove.loops = TRUE) # We create the network where multiple links are kept, while loops are removed.
gh_1970_1989_rmmul_kploop <- igraph::simplify(gh_1970_1989_kpmul_kploop, remove.multiple = TRUE, remove.loops = FALSE) # We create the network where multiple links are removed, while loops are kept.
gh_1970_1989_rmmul_rmloop <- igraph::simplify(gh_1970_1989_kpmul_kploop, remove.multiple = TRUE, remove.loops = TRUE) # We create the network where both the multiple links and loops are removed.

# Next we remove all the isolated nodes out of the network of gh_1970_1989_kpmul_kploop 
isolated_nodes_1970_1989_from_kpmulkploop <- V(gh_1970_1989_kpmul_kploop)[which(degree(gh_1970_1989_rmmul_rmloop,mode="out")==0 & degree(gh_1970_1989_rmmul_rmloop,mode="in")==0)]
gh_1970_1989_kpmul_kploop_rmisonodes <- gh_1970_1989_kpmul_kploop - isolated_nodes_1970_1989_from_kpmulkploop
gh_1970_1989_kpmul_rmloop_rmisonodes <- igraph::simplify(gh_1970_1989_kpmul_kploop_rmisonodes, remove.multiple = FALSE, remove.loops = TRUE) # We remove all the isolated nodes out of the network of gh_1970_1989_kpmul_rmloop 
gh_1970_1989_rmmul_kploop_rmisonodes <- igraph::simplify(gh_1970_1989_kpmul_kploop_rmisonodes, remove.multiple = TRUE, remove.loops = FALSE) # We remove all the isolated nodes out of the network of gh_1970_1989_rmmul_kploop 
gh_1970_1989_rmmul_rmloop_rmisonodes <- igraph::simplify(gh_1970_1989_kpmul_kploop_rmisonodes, remove.multiple = TRUE, remove.loops = TRUE) # We remove all the isolated nodes out of the network of gh_1970_1989_rmmul_rmloop 



#### Next we capture the giant component out of the overall 1970-1999 network ####
get.giant.component <- function(g) {
  
  if (!is.igraph(g)) stop("g is not an igraph object")
  
  comps <- clusters(g, mode="weak")
  induced.subgraph(g, which(which.max(comps$csize)==comps$membership))
} # get giant component ##
gh_giant = get.giant.component(gh_rmmul_rmloop_1970_1989)










