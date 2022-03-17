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
set.seed(100)
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
# a function to calculate h
cal_h <- function(gh, network_size){ 
  h_ratio_list = c()
  for(i in 1:network_size){ # check the node's out-neighbours
    # get born-isolated nodes
    loop_index = which_loop(gh)
    loop = E(gh)[loop_index]
    born_iso_node = head_of(gh,loop) # born-isolated nodes = nodes who have loops
    
    nb = neighbors(gh, i, mode = c("out")) # get out-neighbours of node i
    cnt_borniso_i_outnb = 0 # count born-isolated nodes in i's out-neighbours
    for(j in nb){ # check the node's out-neighbours' out-neighbours
      if(is.element(j, born_iso_node)){
        cnt_borniso_i_outnb=cnt_borniso_i_outnb+1
      }  # TRUE means this out-neighbour of i is born-isolated
    }
    cnt_borncol_i_outnb = length(nb) - cnt_borniso_i_outnb # calculate the number of born-collaborated nodes in i's out-neighbours
    h_ratio_list[i] = cnt_borncol_i_outnb/length(nb) #the fraction of born-collaborated nodes in one’s out-neighbors
  }
  h_ratio_list = na.omit(h_ratio_list) # remove na from h_ratio_list
  h_ratio = sum(h_ratio_list)/network_size
}
nNets = 100
mat_outdegree = matrix(0 , nrow = nNets, ncol = 4 )
mat_CTT = matrix(0 , nrow = nNets, ncol = 9 )

for (net in 1:nNets){
  ############################################################################
  ## Simulation: create a network following the Integrated Edge-based model ##
  ############################################################################
  nRounds <- 1 # the number of networks we are going to create
  nNodes <- 1000 # in total we will have nodes number nNodes
  initial_nActors <- 10  #initially number of nodes
  indegree <- matrix(0 , nrow = nNodes, ncol = nRounds) # create a matrix of all 0, to store indegrees
  K <- 1#10 #number of searching process of a new node
  p_1<- runif(1)#random generate
  q_1<- runif(1)#random generate
  r_1 <- runif(1)#random generate
  trans_triple<- make_graph(c(1, 2, 2, 3, 1,3), directed = TRUE) # create transitive triples, will be used in graph.get.subisomorphisms.vf2 
  trans_triplets<- make_graph(c(1, 2, 2, 3), directed = TRUE) # create triplets, will be used in graph.get.subisomorphisms.vf2 
  
  indegree_arr <- vector() 
  
  for(rounds in 1:nRounds){ # 1st level loop networks
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
    ntriple_total <- length(graph.get.subisomorphisms.vf2(ini_gh,trans_triple)) # firstly add the number of triples of initial network
    ntriple_fa_mo_overlap <-0
    
    for(nActors in (initial_nActors+1):nNodes){    #2nd level loop new node: add the 11st 12nd ... 100th nodes, if initially 10 nodes; totally 100 nodes)
      gh <- add_vertices(gh, 1) # add the new node to gh, not yet create any connection
      V(gh)$name[nActors] <- nActors # name each new node
      redge <- sample(E(gh), size = K, replace = TRUE) # randomly search K edges
      pt_mo <- numeric(K) # create a vector to store the potential mother nodes found and maybe connected by K searching processes
      pt_fa <- numeric(K) # create a vector to store the potential father nodes found and maybe connected by K searching processes
      pt_sm <- list() #create a list to store the potential stepmother nodes found and maybe connected by K searching processes
      
      for(k_find in 1:K){             #3rd level loop of search processes: count the number of three types of nodes( mother, father, stepmother), in K multiple searching processes allNodes[(1=mother, 2=father, 3=stepmother ),k,nActors,nRounds]
        pt_mo[k_find] <- head_of(gh, redge[k_find]) # fill the vector with  potential mother nodes' name
        pt_fa[k_find] <- tail_of(gh, redge[k_find]) # fill the vector with  potential father nodes' name
        fa_outnb_dif_mo <- setdiff( neighbors(gh, pt_fa[k_find], mode = c("out")), pt_mo[k_find] ) # temporarily store father node's outneighbour nodes except the mother node
        pt_sm <- c(pt_sm, list( fa_outnb_dif_mo)) # combine lists. Each list is a vector of potential stepmothers found in one searching process of one new node
      }
      
      mo_new_K_sp<-numeric() # 0 to K items. store the names of nodes that are connected as mother node of this new node
      fa_new_K_sp<-numeric() # 0 to K items. store the names of nodes that are connected as father node of this new node
      sm_new_K_sp<-numeric() # 0 to (K*(potential stepmother-1))several items. store the name of nodes that are connected as stepmother node of this new node
      ntriple_nAc_fa_mo_K_sp <- 0 # 0:K. number of "new-father-mother" triples. Count the number of triples in all searching processes of one node. After the for loop below, this item is the number of triples in all searching process of one node
      ntriple_nAc_fa_sm_K_sp <- 0 # 0:K*(m-1). number of "new-father-stepmother" triples
      ntriple_nAc_mo_sm_K_sp <- 0 # 0:K*(m-1). number of "new-mother-stepmother" triples
      ntriple_nAc_sm_sm_K_sp <- 0 # 0:K*(m(m-1)/2). number of "new-stepmother-stepmother" triples
      
      for(k in 1:K){ # parallel to 3rd loop level of K search processes. for each randomly found edge, connect three kinds of nodes
        
        mo_new_one_sp <- numeric(0) #The name of mother node. 0 or 1.  Only in 1 searching process, there may create 1 or 0 mother node. At beginning we set it as none
        fa_new_one_sp <- numeric(0) #The name of father node. 0 or 1. same as above line
        sm_new_one_sp <- numeric(0) #The names of stepmother nodes. 0 to (potential stepmother-1) nodes.  only in one searching process, there may be multiple stepmother node.
        # sm_new_one_sp_sub <- numeric(0) # 0 or 1 node, whether a node found as potential stepmother is connected as true stepmother
        
        # connect potential mother node
        p <- rbinom(1,1,p_1)                                  #generate a number, 0 or 1, prob(gen 1) = p_1, attached to mother node
        if(p == 1){                                           # with probability p_1, add edge "New ->mother"
        gh <- add_edges(gh, c(nActors, pt_mo[k]))           # add edge "New ->mother"
        mo_new_one_sp <- pt_mo[k]                           # if p==1, then new node connects a node as mother node
        mo_new_K_sp <- c( mo_new_K_sp, mo_new_one_sp )     #combine mother node that are connected by one new node into one vector
        }    
        # connect potential father node
        q <- rbinom(1,1,q_1)                                   #generate a number, 0 or 1, prob(gen 1) = q_1, attached to father node
        if(q == 1){    #if q==1, add edge "New ->father"
        gh <- add_edges(gh, c(nActors, pt_fa[k]))           # add edge "New ->father"
        fa_new_one_sp <- pt_fa[k]                           # if q==1, then new node connects a node as father node
        fa_new_K_sp <- c( fa_new_K_sp, fa_new_one_sp )      #combine father node that are connected by one new node into one vector
        }     
        
        sm_new_one_sp_sub <- numeric(0) #0 to length(pt_sm[[k]]), whether a node found as potential stepmother is connected as true stepmother
        ntriple_nAc_fa_sm_one_sp <- 0 #The number of "New->father->stepmother". 0 to length(pt_sm[[k]]. For one searching process, how many triples "new->father->stepmother" are formed
        ntriple_nAc_mo_sm_one_sp <- 0 #The number of "New->mother->stepmother". 0 to length(pt_sm[[k]]. For one searching process, how many triples "new->mother->stepmother" are formed
        
        if(length(pt_sm[[k]]) == 0){gh <- gh} else for( lc in 1:length(pt_sm[[k]]) ) { # if no potential stepmother nodes, do nothing. Otherwise, new node tries to connect to them. 
          # lc indexes K potential stepmother nodes
          r <- rbinom(1,1,r_1) #generate a number, 0 or 1, prob(gen 1) = r_1, the probability of successfully connect to ONE stepmother node
          if(r == 1){ 
            gh <- add_edges(gh, c(nActors, pt_sm[[k]][lc]))
            sm_new_one_sp_sub <- pt_sm[[k]][lc] #temporarily store the stepmother node
            subgh_nAc_fa_sm_one_subsp <- induced_subgraph(gh, c(nActors, fa_new_one_sp, sm_new_one_sp_sub) ) # for counting triples of "New->father->stepmother"
            subgh_nAc_mo_sm_one_subsp <- induced_subgraph(gh, c(nActors, mo_new_one_sp, sm_new_one_sp_sub) ) # for counting triples of "New->mother->stepmother"
            nfs_triple_0_1_temp <- length( graph.get.subisomorphisms.vf2(subgh_nAc_fa_sm_one_subsp,trans_triplets) )# 0 or 1, if there is a triple "New->father->stepmother"
            nms_triple_0_1_temp <- length( graph.get.subisomorphisms.vf2(subgh_nAc_mo_sm_one_subsp,trans_triplets) ) # 0 or 1, if there is a triple "New->mother->stepmother"
            ntriple_nAc_fa_sm_one_sp <- ntriple_nAc_fa_sm_one_sp + nfs_triple_0_1_temp # "New->father->stepmother" triples +1
            ntriple_nAc_mo_sm_one_sp <- ntriple_nAc_mo_sm_one_sp + nms_triple_0_1_temp # "New->mother->stepmother" triples +1
            sm_new_one_sp<- c(sm_new_one_sp, sm_new_one_sp_sub) #stepmother nodes +1 (name)
          } #if(r == 1) 
        } # for( lc in 1:length(pt_sm[[k]]) )
        subgh_nAc_fa_mo_one_sp <- induced_subgraph(gh, c(nActors ,mo_new_one_sp, fa_new_one_sp )) # catch a subgraph, which is consisted of new, father & mother nodes found by a new node in this searching process. 
        ntriple_nAc_fa_mo_K_sp <- ntriple_nAc_fa_mo_K_sp + length(graph.get.subisomorphisms.vf2(subgh_nAc_fa_mo_one_sp,trans_triple)) # "New->father->mother" +1. 
        ntriple_nAc_fa_sm_K_sp <- ntriple_nAc_fa_sm_K_sp + ntriple_nAc_fa_sm_one_sp #number of "New->father->stepmother" triples in this search
        ntriple_nAc_mo_sm_K_sp <- ntriple_nAc_mo_sm_K_sp + ntriple_nAc_mo_sm_one_sp #number of "New->mother->stepmother" triples in this search
        subgh_nAc_sm_sm_one_sp <- induced_subgraph(gh, c(nActors, sm_new_one_sp) ) #catch a subgraph of New and stepmother nodes of this searching process
        triples_nAc_sm_sm_one_sp <- graph.get.subisomorphisms.vf2(subgh_nAc_sm_sm_one_sp ,trans_triplets)  # triples of new and stepmother     
        subgh_only_sm_sm_one_sp <- induced_subgraph(gh, c( sm_new_one_sp) ) #catch a subgraph of stepmother nodes
        triples_only_sm_sm_one_sp <- graph.get.subisomorphisms.vf2(subgh_only_sm_sm_one_sp ,trans_triplets)  # name of "stepmother-stepmother-stepmother" triples
        ntriple_nAc_sm_sm_K_sp <- ntriple_nAc_sm_sm_K_sp + length( setdiff(triples_nAc_sm_sm_one_sp, triples_only_sm_sm_one_sp) )  # obtain number of triples of "new-stepmother-stepmother", by subtract "stepmother-stepmother-stepmother" triples from triples_nAc_sm_sm_one_sp
        
        sm_new_K_sp <- c( sm_new_K_sp, sm_new_one_sp ) # names of stepmother nodes
      } # K search processes
      ntriple_fa_mo <- ntriple_fa_mo + ntriple_nAc_fa_mo_K_sp # count "New->father->mother" triples number to total number
      ntriple_fa_sm <- ntriple_fa_sm + ntriple_nAc_fa_sm_K_sp # count "New->father->stepmother" triples
      ntriple_mo_sm <- ntriple_mo_sm + ntriple_nAc_mo_sm_K_sp # count "New->mother->stepmother" triples
      ntriple_sm_sm <- ntriple_sm_sm + ntriple_nAc_sm_sm_K_sp # count "New->stepmother->stepmother" triples
      
      if (length(fa_new_K_sp)+length(mo_new_K_sp)+length(sm_new_K_sp)==0) { #if this node is born isolated, i.e., connects no one, then create a self loop
        gh <- add_edges(gh, c(nActors, nActors)) 
        }
      
    } # FOR loop of new nodes
  
    indegree_arr <- cbind(indegree_arr, degree(gh, mode = "in")) # store indegree for average indegree calculation
  } ##FOR loop of nRounds of networks
  
  #file_name <- paste0("IEmodel","K",K ,"p",substr(p_1, 3,5),"q",substr(q_1, 3,5),"r",substr(r_1, 3,5), ".txt")
  #write_graph(gh, file_name)
  
  mu_m = mean(degree(gh, mode = c("out"))) # average outdegree
  P_i = sum(which_loop(gh))/nNodes # total number of born-isolated nodes = total number of loops 
  P_c = 1-P_i
  var_m = var(degree(gh,mode = "out")) #variance of out-degree
  mu_novar = (K*(p_1+q_1-r_1))/(2*(1-K*r_1))
  mu_EB = mu_novar + sqrt((mu_novar)^2 + K*r_1/(1-K*r_1)*var_m) 
  mu_IE = (P_c-1+K*P_c*(p_1+q_1-r_1)+sqrt((P_c-1+K*P_c*(p_1+q_1-r_1))^2 + 4*K*P_c*(1-K*r_1)*((p_1+q_1)*(1-P_c)+P_c*r_1*var_m)))/(2*P_c*(1-K*r_1)) 
  mat_outdegree[net,1] = mu_m
  mat_outdegree[net,2] = mu_novar
  mat_outdegree[net,3] = mu_EB
  mat_outdegree[net,4] = mu_IE
  
  nEdges = gsize(gh)
  h_m = cal_h(gh,network_size = nNodes)
  Em3 = mean((degree(gh,mode = "out"))^3)
  Em2 = mean((degree(gh,mode = "out"))^2)
  CTT_sim = transitive_triples_clustering(gh)
  CTT_novar = (K*q_1*mu_m*(p_1-r_1)+K*q_1*r_1*(mu_m^2+0))/(mu_m^3-2*K*p_1*r_1*(mu_m^2+0)-K*r_1^2*(mu_m^3+3*mu_m*0-2*mu_m*0-2*mu_m^2-2*0))
  CTT_EB = (K*q_1*mu_m*(p_1-r_1)+K*q_1*r_1*(mu_m^2+var_m))/(mu_m^3-2*K*p_1*r_1*(mu_m^2+var_m)-K*r_1^2*(mu_m^3+3*mu_m*var_m-2*mu_m^2-2*var_m))
  CTT_IE = (K*P_c*q_1*mu_m*(p_1-r_1)+K*P_c*q_1*r_1*(mu_m^2+var_m))/(mu_m^3*h_m-2*K*P_c*p_1*r_1*(mu_m^2+var_m)-K*P_c*r_1^2*((mu_m^3+3*mu_m*var_m)-2*(mu_m^2+var_m)))
  eta = p_1+(mu_m-1)*r_1
  mat_CTT[net,1] = CTT_sim # 1st columns store simulation result
  mat_CTT[net,2] = CTT_novar # 2nd columns store calculations with no variance
  mat_CTT[net,3] = CTT_EB # 3rd columns store calculations of EB model
  mat_CTT[net,4] = CTT_IE # 4th columns store calculations of IE model
  mat_CTT[net,5] = eta
  mat_CTT[net,6] = p_1
  mat_CTT[net,7] = r_1
  mat_CTT[net,8] = mu_m
  mat_CTT[net,9] = var_m
}

df_outdegree = as.data.frame(mat_outdegree)
colnames(df_outdegree) = c('Simulation','Meanfield','EB model','IE model')
#write.csv(df_outdegree,file="/Users/HaoLI/Dropbox/Network/Network2019/QJ/data/IE_random_sim_outdegree.csv",row.names=TRUE) # drops the rownames
df_CTT = as.data.frame(mat_CTT)
colnames(df_CTT) = c('Simulation','Mean field','EB model','IE model','eta', 'p','r','mu','var_m')
#write.csv(df_CTT,file="/Users/HaoLI/Dropbox/Network/Network2019/QJ/data/IE_random_sim_CTT.csv",row.names=TRUE) # drops the rownames
# boxplot, 2 figures. First is for outdegree, second is for clustering
boxplot(df_outdegree[,c(1,2,4)], outline = F, names = c('Simulation','Mean field','IE model'), 
             ylab = "Expected out-degree") #main = "Comparison of out-degree with and without variance (IE model)",
res_outdegree_IE <- boxplot(df_outdegree[,c(1,2,4)], outline = F)
res_outdegree_IE$stats
boxplot(df_CTT[,c(1,2,4)], outline = F, names = c('Simulation','Mean field','IE model'),
             ylab = "Clustering coefficient of transitive triples") #main = "Comparison of clustering with and without variance (IE model)",
res_CTT_IE <- boxplot(df_CTT[,c(1,2,4)], outline = F)
res_CTT_IE$stats


