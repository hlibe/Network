# This file corresponds to the method in section "Fitting the Model to Data". We create directed networks based on the seniority of all authors.
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  

####Load the following packages first, and set the working directory####
library(igraph)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


#### We create a 1970-1989 directed network. Indexing from 1 ####
#### First we create a 1970-1989 sorted author namelist based on publication time, the oldest is index 1 
get_network_1970_1989 <- function(){
  data_origin <- read.csv("networkdatac_1970_1989.csv",header=TRUE) 
  data_pubtime <- data_origin[order(data_origin$year, data_origin$journalid, data_origin$articleid ),] #create a table where sort the rows of data in the order of: 1. publication year, 2.journal id, 3. article id
  data_pubtime$auth1 <- data_pubtime$auth1+1 # add 1 to every author's index, such that to avoid index 0 
  data_pubtime$auth2 <- data_pubtime$auth2+1 # add 1 to every author's index, such that to avoid index 0
  data_pubtime$auth3 <- data_pubtime$auth3+1 # add 1 to every author's index, such that to avoid index 0
  write.csv(data_pubtime,"networkdata_1970_1989_sort_inorderof_year_journal_article_time.csv") # store this dataset
  name_list <- vector()  # create a name list which stores authors uniquely in order of time. The senior ones have smaller indices.
  for(i in 1:nrow(data_pubtime)){
    if(data_pubtime$nauthors[i] == 1){ # If this article is single-authored, then add this author into the namelist if the author is not in the namelist
      if( is.element(data_pubtime$auth1[i],name_list) == FALSE  ) { 
        name_list <- c(name_list, data_pubtime$auth1[i])
      } 
    }
    
    if(data_pubtime$nauthors[i] == 2){ # If this artical is double-authored, then check the first and second author sequentially.
      if( is.element(data_pubtime$auth1[i],name_list) == FALSE  ) { # If the first author is not in the namelist, then add the first author into the namelist
        name_list <- c(name_list, data_pubtime$auth1[i])  
      }
      if( is.element(data_pubtime$auth2[i],name_list) == FALSE ) { # If the second author is not in the namelist, then add the second author into the namelist. 
        name_list <- c(name_list, data_pubtime$auth2[i])  
      }
    }
    
    if(data_pubtime$nauthors[i] == 3){ # If this artical is triple-authored, then check the first, second, and third author sequentially.
      if( is.element(data_pubtime$auth1[i],name_list) == FALSE  ) { # If the first author is not in the namelist, then add the first author into the namelist.
        name_list <- c(name_list, data_pubtime$auth1[i]) 
      }
      if( is.element(data_pubtime$auth2[i],name_list) == FALSE ) { # If the second author is not in the namelist, then add the second author into the namelist.
        name_list <- c(name_list, data_pubtime$auth2[i])  
      }
      if( is.element(data_pubtime$auth3[i],name_list) == FALSE) { # If the third author is not in the namelist, then add the third author into the namelist.
        name_list <- c(name_list, data_pubtime$auth3[i]) 
      }
    }   
  }
  ####  According to the namelist we create a directed network of 1970-1989 
  name_list_char <- as.character(name_list) # We change name list from numeric to charactors.
  g1 <- graph(c(name_list_char[1],name_list_char[1]), isolates = c(name_list_char)) ## This creates a network without edges, and the nodes are named following name_list_char. 
  ## We create c(name_list_char[1],name_list_char[1]) as the initial edge because function graph() must have at least one edge. So  This edge is added manipulately, which will have little influence in the following analysis. 
  for(i in 1:30){ # This for-loop is time-consuming, costing 1~2 hours.
    if(data_pubtime$nauthors[i] == 1){ # If the article is single-authored, then create a self-loop.
      g1 <- g1 + edges(as.character(data_pubtime$auth1[i]), as.character(data_pubtime$auth1[i])) 
    }
    if(data_pubtime$nauthors[i] == 2){  # If the article is double-authored, then we compare and know which author is more senior according to the vector of name_list. Consequently, we create a directed edge: junior -> senior.
      if( which(name_list == data_pubtime$auth1[i]) < which(name_list == data_pubtime$auth2[i]) ){ # If the first author is senior, then we create a directed edge 2nd-author -> 1st-author.
        g1 <- g1 + edges(as.character(data_pubtime$auth2[i]), as.character(data_pubtime$auth1[i]))
      } else{ 
        g1 <- g1 + edges(as.character(data_pubtime$auth1[i]), as.character(data_pubtime$auth2[i])) # Otherwise, we create a directed edge: 1st-author -> 2nd-author.
      }
    }
    if(data_pubtime$nauthors[i] == 3){ #If the article is triple-authored, then we compare and know which author is more senior according to the vector of name_list. 
      if( which(name_list == data_pubtime$auth1[i]) < which(name_list == data_pubtime$auth2[i]) ){ #We compare the 1st-author and the 2nd-author. If the 1st-author is more senior than the 2nd-author, then 2nd-author -> 1st-author.
        g1 <- g1 + edges(as.character(data_pubtime$auth2[i]), as.character(data_pubtime$auth1[i]))
      } else{
        g1 <- g1 + edges(as.character(data_pubtime$auth1[i]), as.character(data_pubtime$auth2[i])) # Otherwize, the edge is in the opposite direction.
      }
      if( which(name_list == data_pubtime$auth1[i]) < which(name_list == data_pubtime$auth3[i]) ){ #We compare the 1st-author and the 3rd-author. If the 1st-author is more senior than the 3rd-author, then 3rd-author -> 1st-author.
        g1 <- g1 + edges(as.character(data_pubtime$auth3[i]), as.character(data_pubtime$auth1[i]))
      } else{
        g1 <- g1 + edges(as.character(data_pubtime$auth1[i]), as.character(data_pubtime$auth3[i])) # Otherwize, the edge is in the opposite direction.
      }
      if( which(name_list == data_pubtime$auth2[i]) < which(name_list == data_pubtime$auth3[i]) ){ #We compare the 2nd-author and the 3rd-author. If the 2nd-author is more senior than the 3rd-author, then 3rd-author -> 2nd-author.
        g1 <- g1 + edges(as.character(data_pubtime$auth3[i]), as.character(data_pubtime$auth2[i]))
      } else{
        g1 <- g1 + edges(as.character(data_pubtime$auth2[i]), as.character(data_pubtime$auth3[i])) # Otherwize, the edge is in the opposite direction.
      }
    }
  }
  
  gh_1970_1989_20190911 <- g1
  
  date_time <- format(Sys.time(), format = "%Y%m%d%H%M") # we create a timestamp
  name_with_time <- paste("digraph_1970_1989_kpmul_kploop_",date_time,"_graphml.txt",sep="") # assemble a file name with the time stamp.
  write_graph(g1, format = "graphml", file = name_with_time) # store the digraph of 1970~1989 in the format of "graphml", such that the every node is stored by its charactor name instead of numeric name.
  
  gh_1970_1989_rmmul_kploop <- igraph::simplify(g1, remove.multiple = TRUE, remove.loops = FALSE) # we simplify the network by removing the mul links of the same pair of nodes, but keeping the loops of all nodes.
  name_with_time <- paste("digraph_1970_1989_rmmul_kploop_",date_time,"_graphml.txt",sep="")
  write_graph(gh_1970_1989_rmmul_kploop,format = "graphml", file = name_with_time)
  
  gh_1970_1989_kpmul_rmloop <- igraph::simplify(g1, remove.multiple = FALSE, remove.loops = TRUE) # we simplify the network by removing loops while keeping mul links.
  name_with_time <- paste("digraph_1970_1989_kpmul_rmloop_",date_time,"_graphml.txt",sep="")
  write_graph(gh_1970_1989_kpmul_rmloop,format = "graphml", file = name_with_time)
  
  gh_1970_1989_rmmul_rmloop <- igraph::simplify(g1, remove.multiple = TRUE, remove.loops = TRUE) # we simplify the network by removing both the mul links and loops. 
  name_with_time <- paste("digraph_1970_1989_rmmul_rmloop_",date_time,"_graphml.txt",sep="")
  write_graph(gh_1970_1989_rmmul_rmloop,format = "graphml", file = name_with_time)
}

isolated_nodes_1970_1989_from_rmmulkploop <- V(gh_1970_1989_rmmul_kploop)[which(degree(gh_1970_1989_rmmul_rmloop,mode="out")==0 & degree(gh_1970_1989_rmmul_rmloop,mode="in")==0)] # from the network without loops, we find the index of isolated nodes 
gh_1970_1989_rmmul_kploop_rmisolated <- gh_1970_1989_rmmul_kploop - isolated_nodes_1970_1989_from_rmmulkploop
name_with_time <- paste("digraph_1970_1989_rmmul_kploop_rmisonodes_",date_time,"_graphml.txt",sep="")
write_graph(gh_1970_1989_rmmul_kploop_rmisolated, format = "graphml", file = name_with_time)

isolated_nodes_1970_1989_from_rmmulrmloop <- V(gh_1970_1989_rmmul_rmloop)[which(degree(gh_1970_1989_rmmul_rmloop,mode="out")==0 & degree(gh_1970_1989_rmmul_rmloop,mode="in")==0)] # from the network without loops, we find the index of isolated nodes 
gh_1970_1989_rmmul_rmloop_rmisolated <- gh_1970_1989_rmmul_rmloop - isolated_nodes_1970_1989_from_rmmulrmloop
name_with_time <- paste("digraph_1970_1989_rmmul_rmloop_rmisonodes_",date_time,"_graphml.txt",sep="")
write_graph(gh_1970_1989_rmmul_rmloop_rmisolated, format = "graphml", file = name_with_time)

## The difference between gh_1970_1989_rmmul_kploop_rmisolated and gh_1970_1989_rmmul_rmloop_rmisolated 
## is that the loops of some collabrated nodes are kept in gh_1970_1989_rmmul_kploop_rmisolated

#### We create an undirected network of 1990-1999  ###################
#The function getnetwork() is by Marco van der Leij           

get_undirected_network <- function(startyear,endyear,data=NULL,y=NULL) {
  
  if (is.null(data)) data<-read.csv("networkdatac.csv",header=TRUE)
  if (is.null(y)) y<-subset(data, year>=startyear & year<(endyear+1)) 
  
  n1<-max(data$auth1+1)
  g1<-graph.empty(n1,directed=FALSE)
  
  a<-subset(y, nauthors==2)  
  g1<-add_edges(g1,rbind(a$auth1,a$auth2))
  
  a2<-subset(y, nauthors==3)
  g1<-add_edges(g1,rbind(a2$auth1,a2$auth2))
  g1<-add_edges(g1,rbind(a2$auth1,a2$auth3))
  g1<-add_edges(g1,rbind(a2$auth2,a2$auth3))
  
  igraph::simplify(g1)
}

#### This function is a variant of Marco van der leji's function, creating the network stored in name of networkdatac_1990_1999.csv 
getnetwork_1990_1999 <- function(startyear,endyear,data=NULL,y=NULL) {
  
  if (is.null(data)) data<-read.csv("networkdatac_1990_1999.csv",header=TRUE)
  if (is.null(y)) y<-subset(data, year>=startyear & year<(endyear+1)) 
  
  n1<-max(data$auth1+1)
  g1<-graph.empty(n1,directed=FALSE)
  
  a<-subset(y, nauthors==2)  
  g1<-add_edges(g1,rbind(a$auth1,a$auth2))
  
  a2<-subset(y, nauthors==3)
  g1<-add_edges(g1,rbind(a2$auth1,a2$auth2))
  g1<-add_edges(g1,rbind(a2$auth1,a2$auth3))
  g1<-add_edges(g1,rbind(a2$auth2,a2$auth3))
  
  igraph::simplify(g1)
}
