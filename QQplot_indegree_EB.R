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


###################################################################
## Simulation: create a network following the Edge-based model ##
###################################################################
getwd() # get the working directory 
nRounds <- 20 # the number of networks we are going to create
nNodes <- 2000 # in total we will have nodes number nNodes
initial_nActors <- 10  #initially number of nodes
indegree <- matrix(0 , nrow = nNodes, ncol = nRounds) # create a matrix of all 0, to store indegrees
K <- 1#10 #number of searching process of a new node
p_1<- 0.8#0.028
q_1<- 0.9#0.044
r_1 <- 0.7 #0.027
trans_triple<- make_graph(c(1, 2, 2, 3, 1,3), directed = TRUE) # create transitive triples, will be used in graph.get.subisomorphisms.vf2 
trans_triplets<- make_graph(c(1, 2, 2, 3), directed = TRUE) # create triplets, will be used in graph.get.subisomorphisms.vf2 

indegree_arr <- vector() 

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

##### store the network, write graph 
#write_graph( gh, "gh_EBmodel.txt")


#### Q-Q plot. Here we use the indegree matrix indegree_arr to draw the Q-Q plot ####################
indegree_mean <- rowMeans(indegree_arr ) # average in-degree for comparison between theoretic result. the indegree_arr is generated from the above network formation process.
#write.csv(indegree_mean,file = "meanindegree_2Knodes20rounds.csv")
DATA <- indegree_mean[5:2000]
N         <- length(DATA);
PERCS     <- ((1:N)-0.5)/N; 

mu <- mean(degree(gh, mode="out"))
sub_exp <- mu/((mu-1)*r_1+p_1)
QUANTILES <- q_1*sub_exp/(1-PERCS)^(K/sub_exp)-q_1*sub_exp; # use inverse function of CDF, calculate Quantile
PLOTDATA <- data.frame(Sample = sort(DATA),
                       Theoretical = QUANTILES);

#Generate custom QQ plot
library(ggplot2);
theme_update(plot.title    = element_text(size = 15, hjust = 0.5),
             plot.subtitle = element_text(size = 10, hjust = 0.5),
             axis.title.x  = element_text(size = 20, hjust = 0.5),
             axis.title.y  = element_text(size = 20, vjust = 0.5),
             axis.text.x = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
             axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
             plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm"));

QQPLOT <- ggplot(data = PLOTDATA, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 2, colour = 'black') + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  #ggtitle('Quantile-Quantile Plot') + 
  #labs(subtitle = '(Comparison to exponential distribution with unit scale)') + 
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles')+ xlim(0,40) + ylim(0,40);

QQPLOT; # a few outliers are the initial small network
ggsave("/Users/HaoLI/Dropbox/Network/Network2019/QJ/images/QQplot_indegree_theory_simulation_EB.pdf", width = 10,   height = 8)

