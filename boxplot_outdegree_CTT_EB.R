# This file is to create directed networks following the models in the paper. Afterwards, we check the indegree distributions using QQ plot.
# This file is writen by Dr. Hao LI, 2018-2019, in PhD study at Department of Industrial Engineering and Decision Analytics,
# Hong Kong University of Science and Technology  

####Load the following packages first, and set the working directory####
#install.packages("igraph")
library(igraph)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

set.seed(10)
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

nNets = 100
mat_outdegree_EB = matrix(0 , nrow = nNets, ncol = 4 )
mat_CTT_EB = matrix(0 , nrow = nNets, ncol = 4 )

for (net in 1:nNets) {
  ###########################################################
  ## Simulation: create a network of the Edge-based model ###
  ###########################################################
  nRounds <- 1 # the number of networks we are going to create
  nNodes <- 1000 # in total we will have nodes number nNodes
  initial_nActors <- 10  #initially number of nodes
  indegree <- matrix(0 , nrow = nNodes, ncol = nRounds) # create a matrix of all 0, to store indegrees
  K <- 1#10 #number of searching process of a new node
  p_1<- runif(1)#random generate
  q_1<- runif(1)
  r_1 <- runif(1)
  trans_triple<- make_graph(c(1, 2, 2, 3, 1,3), directed = TRUE) # create transitive triples, will be used in graph.get.subisomorphisms.vf2 
  trans_triplets<- make_graph(c(1, 2, 2, 3), directed = TRUE) # create triplets, will be used in graph.get.subisomorphisms.vf2 
  
  indegree_arr <- vector() 
  ptm <- proc.time() 
  for(rounds in 1:nRounds){ # 1st loop
    mat <- matrix(0, initial_nActors, initial_nActors)  #create a matrix with all 0
    mat[lower.tri(mat)] <- 1            #change the all 0 matrix into a matrix with lower triple 1, hence to make a network "2->1, 3->1, 3->2, 4->1, 4->2, 4->3, 5->1, 5->2, 5->3, 5->4,..."
    gh <- graph_from_adjacency_matrix(mat)  #create the initial network
    V(gh)$name <- 1:(initial_nActors) #name the initial nodes
    E(gh)$name <- 1:(initial_nActors*(initial_nActors-1)/2) #name the edges: 2->1 3->1 3->2 4->1 4->2 4->3 5->1 5->2 5->3 5->4 6->1 6->2 6->3 6->4 6->5 7->2 7->6 7->3 7->6 7->3
    ini_gh <- gh  #store initial graph independently
    
    ntriple_fa_mo <- 0 # 0:K*nNodes
    ntriple_fa_sm <- 0 # 0:K*(m-1)
    ntriple_sm_sm <- 0 # 0:K*(m(m-1)/2)
    ntriple_mo_sm <- 0 # 0:K*(m-1)
    ntriple_total <- length(graph.get.subisomorphisms.vf2(ini_gh,trans_triple)) # firstly add the number of triple of initial network
    ntriple_fa_mo_overlap <-0
    
    for(nActors in (initial_nActors+1):nNodes){    #2nd loop: add the 7th 8th ... 100th nodes (if initially 6 nodes, totally 100 nodes)
      gh <- add_vertices(gh, 1) 
      V(gh)$name[nActors] <- nActors # name each new node
      redge <- sample(E(gh), size = K, replace = TRUE)
      pt_mo <- numeric(K) # create a vector to store the nodes found and maybe connected by K searching processes
      pt_fa <- numeric(K) # create a vector to store the nodes found and maybe connected by K searching processes
      pt_sm <- list() #create a vector of length 0
      
      for(k_find in 1:K){                                           #3rd loop: count the number of three tpyes of nodes( mother, father, stepmother), in K multiple searching processes allNodes[(1=mother, 2=father, 3=stepmother ),k,nActors,nRounds]
        pt_mo[k_find] <- head_of(gh, redge[k_find]) # fill the vector with  potential mother nodes' name
        pt_fa[k_find] <- tail_of(gh, redge[k_find]) # fill the vector with  potential father nodes' name
        fa_outnb_dif_mo <- setdiff( neighbors(gh, pt_fa[k_find], mode = c("out")), pt_mo[k_find] ) # temperaly store one father node's out neighbour nodes except the mother node
        pt_sm <- c(pt_sm, list( fa_outnb_dif_mo)) # combine lists. Each list is a vector of potential stepmothers found in one searching process of one new node
      }
      
      mo_new_K_sp<-numeric() # 0:K items. For a new node, create a vector to store the name of nodes that are connected as mother node of a new node
      fa_new_K_sp<-numeric() # 0:K items. For a new node, create a vector to store the name of nodes that are connected as father node of a new node
      sm_new_K_sp<-numeric() # 0:(K*(potential stepmother-1))several items. For a new node, create a vector to store the name of nodes that are connected as stepmother node of a new node
      ntriple_nAc_fa_mo_K_sp <- 0 # 0:K. Count the number of triples in all searching processes of one node. After the for loop below, this item is the number of triples in all searching process of one node
      ntriple_nAc_fa_sm_K_sp <- 0 # 0:K*(m-1)
      ntriple_nAc_mo_sm_K_sp <- 0 # 0:K*(m-1)
      ntriple_nAc_sm_sm_K_sp <- 0 # 0:K*(m(m-1)/2)
      
      for(k in 1:K){
        
        mo_new_one_sp <- numeric(0) #one node or 0. Only in one searching process, there may be 1 or 0 mother node. At beginning we set it as none
        fa_new_one_sp <- numeric(0) #one node or 0. same as above line
        sm_new_one_sp <- numeric(0) #(potential stepmother-1) nodes to 0.  only in one searching process, there may be multiple or 0 stepmother node.
        sm_new_one_sp_sub <- numeric(0) # one node or 0, whether a node found as potential stepmother is connected as true stepmother
        
        p <- rbinom(1,1,p_1)                                  #generate a number, 0 or 1, prob(1) = 0.1, attached to mother node
        if(p == 1){ gh <- add_edges(gh, c(nActors, pt_mo[k])) #add edge "New ->mother"
        mo_new_one_sp <- pt_mo[k]                           # if p==1, then new node connects a node as mother node
        mo_new_K_sp <- c( mo_new_K_sp, mo_new_one_sp ) }    #combine mother node that are connected by one new node into one vector
        
        q <- rbinom(1,1,q_1)                                   #generate a number, 0 or 1, prob(1) = 0.5, attached to father node
        if(q == 1){ gh <- add_edges(gh, c(nActors, pt_fa[k])) #add edge "New ->father"
        fa_new_one_sp <- pt_fa[k]                           # if q==1, then new node connects a node as father node
        fa_new_K_sp <- c( fa_new_K_sp, fa_new_one_sp ) }     #combine father node that are connected by one new node into one vector
        
        sm_new_one_sp_sub <- numeric(0) #0:length(pt_sm[[k]]
        ntriple_nAc_fa_sm_one_sp <- 0 #0:length(pt_sm[[k]]. For one searching process, how many triples "new->father->stepmother" are formed
        ntriple_nAc_mo_sm_one_sp <- 0 #0:length(pt_sm[[k]]. For one searching process, how many triples "new->mother->stepmother" are formed
        
        if(length(pt_sm[[k]]) == 0){gh <- gh} else for( lc in 1:length(pt_sm[[k]]) ) { 
          r <- rbinom(1,1,r_1) #generate a number, 0 or 1, prob(1) = r_1, the probability of successfully connect to ONE stepmother node
          if(r == 1){ 
            gh <- add_edges(gh, c(nActors, pt_sm[[k]][lc]))
            sm_new_one_sp_sub <- pt_sm[[k]][lc] #temprarily store the stepmother node
            subgh_nAc_fa_sm_one_subsp <- induced_subgraph(gh, c(nActors, fa_new_one_sp, sm_new_one_sp_sub) )
            subgh_nAc_mo_sm_one_subsp <- induced_subgraph(gh, c(nActors, mo_new_one_sp, sm_new_one_sp_sub) )
            nfs_triple_0_1_temp <- length( graph.get.subisomorphisms.vf2(subgh_nAc_fa_sm_one_subsp,trans_triplets) )# 0 or 1, if there is a triple "New->father->stepmother"
            nms_triple_0_1_temp <- length( graph.get.subisomorphisms.vf2(subgh_nAc_mo_sm_one_subsp,trans_triplets) ) # 0 or 1, if there is a triple "New->mother->stepmother"
            ntriple_nAc_fa_sm_one_sp <- ntriple_nAc_fa_sm_one_sp + nfs_triple_0_1_temp # count all "New->father->stepmother" triples in one searching process
            ntriple_nAc_mo_sm_one_sp <- ntriple_nAc_mo_sm_one_sp + nms_triple_0_1_temp
            sm_new_one_sp<- c(sm_new_one_sp, sm_new_one_sp_sub) #if r==1, we store one more node into the stepmother nodes new node found in one searching process
          } #if(r == 1) 
        } # for( lc in 1:length(pt_sm[[k]]) )
        subgh_nAc_fa_mo_one_sp <- induced_subgraph(gh, c(nActors ,mo_new_one_sp, fa_new_one_sp )) # catch a subgraph, which is consisted of only new nodes and father & mother nodes found by a new node in only one searching process. 
        ntriple_nAc_fa_mo_K_sp <- ntriple_nAc_fa_mo_K_sp + length(graph.get.subisomorphisms.vf2(subgh_nAc_fa_mo_one_sp,trans_triple)) # if the subgraph captured from the sentence above contains a triple, then count +1
        ntriple_nAc_fa_sm_K_sp <- ntriple_nAc_fa_sm_K_sp + ntriple_nAc_fa_sm_one_sp #after each searching process, how many "New->father->stepmother" triples are made
        ntriple_nAc_mo_sm_K_sp <- ntriple_nAc_mo_sm_K_sp + ntriple_nAc_mo_sm_one_sp #after each searching process, how many "New->mother->stepmother" triples are made
        subgh_nAc_sm_sm_one_sp <- induced_subgraph(gh, c(nActors, sm_new_one_sp) ) #after one searching process, the subgraph made by New node and stepmother nodes of this searching process
        triples_nAc_sm_sm_one_sp <- graph.get.subisomorphisms.vf2(subgh_nAc_sm_sm_one_sp ,trans_triplets)       
        subgh_only_sm_sm_one_sp <- induced_subgraph(gh, c( sm_new_one_sp) )
        triples_only_sm_sm_one_sp <- graph.get.subisomorphisms.vf2(subgh_only_sm_sm_one_sp ,trans_triplets)  
        ntriple_nAc_sm_sm_K_sp <- ntriple_nAc_sm_sm_K_sp + length( setdiff(triples_nAc_sm_sm_one_sp, triples_only_sm_sm_one_sp) )  
        
      } # Searching process
      ntriple_fa_mo <- ntriple_fa_mo + ntriple_nAc_fa_mo_K_sp # sum "New->father->mother" triples number to total number
      ntriple_fa_sm <- ntriple_fa_sm + ntriple_nAc_fa_sm_K_sp
      ntriple_mo_sm <- ntriple_mo_sm + ntriple_nAc_mo_sm_K_sp
      ntriple_sm_sm <- ntriple_sm_sm + ntriple_nAc_sm_sm_K_sp
      
    } # FOR loop of nActors
  
    indegree_arr <- cbind(indegree_arr, degree(gh, mode = "in")) # store indegree for average indegree calculation
  } ##FOR loop of nRounds
  proc.time() - ptm
  #file_name <- paste0("EBmodel","K",K ,"p",substr(p_1, 3,5),"q",substr(q_1, 3,5),"r",substr(r_1, 3,5), ".txt")
  #write_graph(gh, file_name)
  
  var_m = var(degree(gh,mode = "out")) #variance of out-degree
  #var_m = gh_chars$variance_of_outdegree
  #gh_chars = get_network_characteristics(gh)
  P_i = sum(which_loop(gh))/nNodes # total number of born-isolated nodes = total number of loops 
  P_c = 1-P_i 
  mu_m = mean(degree(gh, mode = c("out"))) # average outdegree
  mu_novar = (K*(p_1+q_1-r_1))/(2*(1-K*r_1))
  mu_EB = mu_novar + sqrt((mu_novar)^2 + K*r_1/(1-K*r_1)*var_m) 
  mu_IE = (P_c-1+K+K*P_c*(p_1+q_1-r_1)+sqrt((P_c-1+K*P_c*(p_1+q_1-r_1))^2 + 4*K*P_c*(1-K*r_1)*((p_1+q_1)*(1-P_c)+P_c*r_1*sqrt(var_m))))/(2*P_c*(1-K*r_1)) 
  mat_outdegree_EB[net,1] = mu_m
  mat_outdegree_EB[net,2] = mu_novar
  mat_outdegree_EB[net,3] = mu_EB
  mat_outdegree_EB[net,4] = mu_IE  

  nEdges = gsize(gh)
  h_m =  (nNodes*P_c*q_1)/(nEdges-nNodes*P_c*(p_1+(mu_m-1)*r_1)) 
  CTT_sim = transitive_triples_clustering(gh)
  CTT_novar = (K*q_1*mu_m*(p_1-r_1)+K*q_1*r_1*(mu_m^2+0))/(mu_m^3-2*K*p_1*r_1*(mu_m^2+0)-K*r_1^2*(mu_m^3+3*mu_m*0-2*mu_m*0-2*mu_m^2-2*0))
  CTT_EB = (K*q_1*mu_m*(p_1-r_1)+K*q_1*r_1*(mu_m^2+var_m))/(mu_m^3-2*K*p_1*r_1*(mu_m^2+var_m)-K*r_1^2*(mu_m^3+3*mu_m*var_m-2*mu_m^2-2*var_m))
  CTT_IE = (q_1*mu_m*(p_1-r_1)+q_1*r_1*(mu_m^2+var_m))/((K*P_c*nNodes)/(mu_m^2*h_m*nEdges)-r_1^2*(r_1*(mu_m^3+3*mu_m*var_m)+(p_1+q_1-3*r_1)*(mu_m^2+var_m)+2*mu_m*(r_1-p_1-q_1))-2*p_1*r_1*(r_1*(mu_m^2+var_m)+(p_1+q_1-r_1)*mu_m))
  mat_CTT_EB[net,1] = CTT_sim
  mat_CTT_EB[net,2] = CTT_novar
  mat_CTT_EB[net,3] = CTT_EB
  mat_CTT_EB[net,4] = CTT_IE  
}
df_outdegree_EB = as.data.frame(mat_outdegree_EB)
colnames(df_outdegree_EB) = c('Simulation','Mean field','EB model','IE model')
#write.csv(df_outdegree_EB,file="/Users/HaoLI/Dropbox/Network/Network2019/QJ/data/EB_random_sim_outdegree.csv",row.names=TRUE) # drops the rownames
df_CTT_EB = as.data.frame(mat_CTT_EB)
colnames(df_CTT_EB) = c('Simulation','Mean field','EB model','IE model')
#write.csv(df_CTT_EB,file="/Users/HaoLI/Dropbox/Network/Network2019/QJ/data/EB_random_sim_CTT.csv",row.names=TRUE) # drops the rownames
boxplot(df_outdegree_EB[,c(1,2,3)], outline = F, names = c('Simulation','Mean field','EB model'), 
         ylab = "Expected out-degree")#main = "Comparison of out-degree with and without variance (EB model)",
res_outdegree_EB <- boxplot(df_outdegree_EB[,c(1,2,3)], outline = F)
res_outdegree_EB$stats
boxplot(df_CTT_EB[,c(1,2,3)], outline = F, names = c('Simulation','Mean field','EB model'), 
        ylab = "Clustering coefficient of transitive triples") # main = "Comparison of clustering with and without variance (EB model)",
res_CTT_EB <- boxplot(df_CTT_EB[,c(1,2,3)], outline = F)
median_CTT = res_CTT_EB$stats
text(x= 1, y= 0.3, labels= round(median_CTT[3,1],3))
text(x= 2, y= 0.19, labels= round(median_CTT[3,2],3))
text(x= 3, y= 0.3, labels= round(median_CTT[3,3],3))

