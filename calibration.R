# This file is to calibrate parameters.
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  

####Load the following packages first, and set the working directory####
install.packages("igraph","minpack.lm")
library(igraph)
library(minpack.lm) # use function nlsLM() nonlinear regression
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )


# Edge-based Model, second tier mu = p+(m-1)*r ####
calibrate_EB_directednet <- function(gh, K, m){
  indegree = degree(gh, mode = "in") # list the value of indegree of all nodes
  unique_indegree = sort(unique(indegree)) # take the unique value of indegree of all nodes an then sort them min --> max
  indegree_values <- 0:max(unique_indegree) # create a vector, which contains all possible value of the indegree. That is 0 --> max
  indegree_cumulative_distribution = degree.distribution(gh,mode = "in", cumulative = TRUE) # obtain the indegree_cumulative_distribution, whose value: 1 -> 0
  dataframe_directednet <-data.frame(indegree_cumulative_distribution,indegree_values) # this dataframe is based on the authors' namelist which is according to their joining time
  
  curve.nlslrc.directednet = nlsLM(indegree_cumulative_distribution ~ ((d_0+m*q_EB/secondtier)/(indegree_values + m*q_EB/secondtier))^(m/(K*secondtier)), #  write down EB model cdf with secondtier = (m-1)*r+p
                                   start = list(q_EB = 0.1,
                                                secondtier = 0.1,
                                                d_0=1),
                                   data = dataframe_directednet)
  parameter<-coef(curve.nlslrc.directednet) 
  #next is to solve p and r, assuming p=r
  
  cat("when K=", K,"\n", "q=", parameter[[1]],  ", second tier=", parameter[[2]], ", d_0=",parameter[[3]] ,"\n", "If p=r, then p=r=",parameter[[2]]/m)
}

##### Integrated Edge-based Model ##### 
calibrate_IE_directednet <- function(gh, K, m,t, E, Pc){
  indegree = degree(gh, mode = "in") # list the value of indegree of all nodes
  unique_indegree = sort(unique(indegree)) # take the unique value of indegree of all nodes an then sort them min --> max
  indegree_values <- 0:max(unique_indegree) # create a vector, which contains all possible value of the indegree. That is 0 --> max
  indegree_cumulative_distribution = degree.distribution(gh,mode = "in", cumulative = TRUE) # obtain the indegree_cumulative_distribution, whose value: 1 -> 0
  dataframe_directednet <-data.frame(indegree_cumulative_distribution,indegree_values) # this dataframe is based on the authors' namelist which is according to their joining time
  curve.nlslrc.directednet = nlsLM(indegree_cumulative_distribution ~ ((d_0+(Pc*m*q/secondtier + (1-Pc)*(q-(m-1)*r)/secondtier))/(indegree_values + (Pc*m*q/secondtier + (1-Pc)*(q-(m-1)*r)/secondtier)))^(E/(K*secondtier*t)), #  write down EB model cdf with secondtier = (m-1)*r+p
                                   start = list(q = 0.1,
                                                secondtier = 0.1,
                                                d_0=1),
                                   data = dataframe_directednet)
  parameter<-coef(curve.nlslrc.directednet) 
  #next is to solve p and r, assuming p=r
  
  cat("when K=", K,"\n", "q=", parameter[[1]],  ", second tier=", parameter[[2]], ", d_0=",parameter[[3]] ,"\n", "If p=r, then p=r=",parameter[[2]]/m)
}

##### Jackson & Rogers Model#####
calibrate_rJR_directed_namelist <- function(gh, m){ #this is for directed network data whose direction is based on the entrance time of all authors
  # construct the data frame
  # prepare data for regression in form of directed graph
  
  indegree = degree(gh, mode = "in") # list the value of indegree of all nodes
  unique_indegree = sort(unique(indegree)) # take the unique value of indegree of all nodes an then sort them min --> max
  indegree_values <- 0:max(unique_indegree) # create a vector, which contains all possible value of the indegree. That is 0 --> max
  indegree_cumulative_distribution = degree.distribution(gh,mode = "in", cumulative = TRUE) # obtain the indegree_cumulative_distribution, whose value: 1 -> 0
  indegree_cumulative_distribution[1] <- 0.99
  dataframe_directednet <-data.frame(indegree_cumulative_distribution,indegree_values) # this dataframe is based on the authors' namelist which is according to their joining time
  dataframe_directednet_1 <- data.frame(indegree_cumulative_distribution[-1],indegree_values[-1])
  #calibrate using nonlinear method
  curve.nlslrc = nlsLM(indegree_cumulative_distribution ~ ((d_0+r_JR*m)/(indegree_values + r_JR*m))^(1+r_JR), # directly write down JR model cdf
                       start = list(
                         d_0 = 1,
                         r_JR = 2),
                       data = dataframe_directednet)
  coef(curve.nlslrc) #r_JR  1.634592;  0.4297748
  
} 

r_JR = coef(curve.nlslrc)
d_0 = 3
CCDF_JR_model <- function(r_JR){
  ((d_0+r_JR*m)/(overall_de_values_but_0/2 + r_JR*m))^(1+r_JR)
}
preds_JR_model <- CCDF_JR_model(r_JR)
preds_JR_model <- preds_JR_model[-1] # remove -Inf
actual <- curvelrc[,1]
actual <- actual[-1]
rss_JR_model <- sum((preds_JR_model - actual) ^ 2) ## residual sum of squares
tss_JR_model <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
rsq_JR_model <- 1 - rss_JR_model/tss_JR_model # R squre = 1 - residual sum-of-squares/total sum-of-squares
rsq_JR_model




#### P_c value ####
## To find the value of P_c (0.5842577 = 53632/129003), we find the born-isolated node ---- whose first article is single-authored. 
data_pubtime <- read.csv("networkdata_sort_inorderof_article_time.csv", header=TRUE) # read the dataset where the articles are ranked in order of time.
born_iso_name_list <- vector() # create an empty vector to store the name list of born-isolated nodes
indept_authors <-vector() #create a vector to store independent authors' names. The independent authors are those who have independent works.
for(i in 1:nrow(data_pubtime)){ #store all these independnt authors
  if(data_pubtime$nauthors[i] == 1 ){ # if the artical is single-authored, then the author is added into the vector of independent authors.
    indept_authors <- c(indept_authors, data_pubtime$auth1[i])
  }
}
unique_indept_authors <- unique(indept_authors) # take the unique names. Totally 75700 names if you operate length(unique_indept_authors).

first_author_vector = as.vector(data_pubtime$auth1) # create a vector storing all first author names, in the order of time, including "NA"
second_author_vector = as.vector(data_pubtime$auth2) # create a vector storing all second author names, in the order of time, including "NA"
third_author_vector = as.vector(data_pubtime$auth3) # create a vector storing all third author names, in the order of time, including "NA"
#### Method 2
born_iso_authors_method2 <- vector() # Create a vector to store the born-isolated authors
for(i in 1:length(unique_indept_authors)){ # We check every item in the vector unique_indept_authors, where every author has published an independent paper.
  index_in_1auth_vector = which(first_author_vector == unique_indept_authors[i])[1]   # In first_author_vector, we find the index of the unique author i. That tells us when author i ever published he(r) first first-authored paper.
  index_in_2auth_vector = which(second_author_vector == unique_indept_authors[i])[1]  # In second_author_vector, we find the index of the unique author i. That tells us when author i ever published he(r) first second-authored paper.
  index_in_3auth_vector = which(third_author_vector == unique_indept_authors[i])[1]   # In third_author_vector, we find the index of the unique author i. That tells us when author i ever published he(r) first third-authored paper.
  if(is.na(index_in_3auth_vector) & is.na(index_in_2auth_vector)){ # If author i is not in index_in_2auth_vector or index_in_3auth_vector, then that means author i never be the second author or third author. Then that means author i never collabrated with others. Therefore, author i is a born isolated node
    born_iso_authors_method2 <- c(born_iso_authors_method2, unique_indept_authors[i] ) # add the author i into the vector born_iso_authors_method2
  }
  if(is.na(index_in_3auth_vector) & is.na(index_in_2auth_vector)==FALSE){ # If author i ever be the second author, but never be the third author
    if(index_in_1auth_vector < index_in_2auth_vector){ # then if the author i became the first author BEFORE becoming the second author, then add author i into born_iso_authors_method2
      born_iso_authors_method2 <- c(born_iso_authors_method2, unique_indept_authors[i] )
    }
  }
  if(is.na(index_in_2auth_vector) & is.na(index_in_3auth_vector)==FALSE){ # if author i ever be the third author, but never be the second author
    if(index_in_1auth_vector < index_in_3auth_vector){ # then if the author i became the first author BEFORE becoming the third author, then add author i into born_iso_authors_method2
      born_iso_authors_method2 <- c(born_iso_authors_method2, unique_indept_authors[i] )
    }
  }
}
length(born_iso_authors_method2) # result: 64475, number of authors whose first paper is independent. 64475/129003=0.4997946
born_iso_and_rminn <- vector() # create an empty vectore to store those authors who have in-neighbors
born_iso_and_rmoutn <- vector() # create an empty vectore to store those authors who have out-neighbors


for (i in 1:length(born_iso_authors_method2)){ # check every item in the vector born_iso_authors_method2 length(born_iso_authors_method2)
  indegree_temp <- neighbors(gh_directed_simplify_rmmul_rmloop, born_iso_authors_method2[i], mode = c("in")) # store the in-neighbours of born_iso_authors_method2[i] from the network gh_rmmul_rmloop
  if(length(indegree_temp) == 0 ){ #If there is no in-neighbours for node born_iso_authors_method2[i], then we add node born_iso_authors_method2[i] into born_iso_and_rminn
    born_iso_and_rminn <- c(born_iso_and_rminn, born_iso_authors_method2[i])
  }
  outdegree_temp <- neighbors(gh_directed_simplify_rmmul_rmloop, born_iso_authors_method2[i], mode = c("out")) # store the out-neighbours of born_iso_authors_method2[i] from the network gh_rmmul_rmloop
  if(length(outdegree_temp) == 0 ){ #If there is no out-neighbours for node born_iso_authors_method2[i], then we add node born_iso_authors_method2[i] into born_iso_and_rmoutn
    born_iso_and_rmoutn <- c(born_iso_and_rmoutn, born_iso_authors_method2[i])
  }
}
length(born_iso_and_rmoutn)#result: 53632, 53632/129003 = 0.4157423. P_c = 1-0.4157423=0.5842577
cat("The probability of born-collaborated nodes P_c =", 1-length(born_iso_and_rmoutn))


