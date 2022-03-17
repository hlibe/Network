a=length(V(gh)) # 9877 nodes; Econlit 129003 for the 1970-1999 overall network
cat("Number of nodes:", a,"\n")
b=gsize(gh) #number of edges, 51971 edges; Econlit 106705
cat("number of edges (not articles):" , b,"\n")
c=mean(degree(gh, v=V(gh),mode="all")) # 1.511972 average degree; Econlit 1.654303
cat("average overall degree:",c,"\n" )
d=mean(degree(gh, v=V(gh),mode="out")) # average outdegree
cat("average outdegree:",d ,"\n")
g=length(which(clusters(gh)$csize ==1 )) # count isolated nodes
cat("The number of isolated nodes:",g,"\n")

indegree = degree(gh, mode = "in") # list the value of indegree of all nodes
unique_indegree = sort(unique(indegree)) # take the unique value of indegree of all nodes an then sort them min --> max
indegree_values <- 0:max(unique_indegree) # create a vector, which contains all possible value of the indegree. That is 0 --> max
indegree_cumulative_distribution = degree.distribution(gh,mode = "in", cumulative = TRUE) # obtain the indegree_cumulative_distribution, whose value: 1 -> 0
dataframe_directednet <-data.frame(indegree_cumulative_distribution,indegree_values) # this dataframe is based on the authors' namelist which is according to their joining time
write_csv(dataframe_directednet,"gh_kpmul_kploop_rmisonodes_3296.csv")

isolated_nodes_from_gh <- V(gh)[which(degree(gh,mode="out")==0 & degree(gh,mode="in")==0)]
gh_rmisonodes <- gh - isolated_nodes_from_gh

isolated_nodes_from_kpmulkploop <- V(gh_kpmul_kploop)[which(degree(gh_kpmul_kploop,mode="out")==0 & degree(gh_kpmul_kploop,mode="in")==0)]
length(isolated_nodes_from_kpmulkploop)
