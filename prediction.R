# This file correspondes to the first subsection of "Fitting the Model to Data"
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  

####Load the following packages first, and set the working directory####
install.packages("igraph","ggplot2")
library(igraph)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


indegree_rmmul_rmloop = degree(gh_rmmul_rmloop, mode = "in") # create a vector to store the indegree of all nodes: 1~129003
unique_indegree_rmmul_rmloop = sort(unique(indegree_rmmul_rmloop))

outdegree_kpmul_kploop = degree(gh_kpmul_kploop, mode = "out")
outdegree_rmmul_kploop = degree(gh_rmmul_kploop, mode = "out")
outdegree_kpmul_rmloop = degree(gh_kpmul_rmloop, mode = "out")
outdegree_rmmul_rmloop = degree(gh_rmmul_rmloop, mode = "out")

# Obtain the observation of CCDF ####
# For EB model, we observe the indegree CDF in gh_rmmul_rmloop_rmisonodes 
indegree_rmmul_rmloop_rmisonodes = degree(gh_rmmul_rmloop_rmisonodes, mode = "in") # We create a vector to store the indegree for every node in the network.
unique_indegree_rmmul_rmloop_rmisonodes = sort(unique(indegree_rmmul_rmloop_rmisonodes)) # We create a vector to store the indegrees with unique values.
indegree_values_rmmul_rmloop_rmisonodes <- 0:max(unique_indegree_rmmul_rmloop_rmisonodes) # We create a vectore 0,1,2,...,max(indegree).
indegree_cdf_rmmul_rmloop_rmisonodes = degree.distribution(gh_rmmul_rmloop_rmisonodes,mode = "in", cumulative = TRUE) # This is a vector of CCDF of indegree.
dataframe_indegree_rmmul_rmloop_rmisonodes <-data.frame(indegree_cdf_rmmul_rmloop_rmisonodes,indegree_values_rmmul_rmloop_rmisonodes) # We construct a dataframe where the values of CCDF are listed corresponding to different indegrees.
write.csv(dataframe_indegree_rmmul_rmloop_rmisonodes, file ="dataframe indegree cdf prediction EB rmmul_rmloop_rmisonodes.csv") # Store the dataframe in .csv and will be used in the Excel file to plot the figures.

# For JR model, we observe the indegree CDF in gh_rmmul_rmloop 
indegree_rmmul_rmloop = degree(gh_rmmul_rmloop, mode = "in") # We create a vector to store the indegree for every node in the network.
unique_indegree_rmmul_rmloop = sort(unique(indegree_rmmul_rmloop)) # We create a vector to store the indegrees with unique values.
indegree_values_rmmul_rmloop <- 0:max(unique_indegree_rmmul_rmloop) # We create a vectore 0,1,2,...,max(indegree).
indegree_cdf_rmmul_rmloop = degree.distribution(gh_rmmul_rmloop,mode = "in", cumulative = TRUE) # This is a vector of CCDF of indegree.
dataframe_indegree_rmmul_rmloop <-data.frame(indegree_cdf_rmmul_rmloop,indegree_values_rmmul_rmloop) # We construct a dataframe where the values of CCDF are listed corresponding to different indegrees.
write.csv(dataframe_indegree_rmmul_rmloop, file = "dataframe indegree cdf prediction JR rmmul_rmloop.csv") # Store the dataframe in .csv and will be used in the Excel file to plot the figures.



