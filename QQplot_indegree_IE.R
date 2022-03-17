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
set.seed(1)


###################################################################
## Simulation: create a network following the Edge-based model ##
###################################################################
getwd() # get the working directory 
nRounds <- 20 # the number of networks we are going to create
nNodes <- 2000 # in total we will have nodes number nNodes
initial_nActors <- 10  #initially number of nodes
indegree <- matrix(0 , nrow = nNodes, ncol = nRounds) # create a matrix of all 0, to store indegrees
K <- 1#10 #number of searching process of a new node
p_1<- 0.5#0.028
q_1<- 0.6#0.044
r_1 <- 0.2 #0.027
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

##### store the network, write graph 
#write_graph( gh, "20210226_IEB_10ininodes10KNodeK1p05q07r03_A.txt")



mean(degree(gh, mode = c("out"))) # average outdegree
sum(which_loop(gh)) # total number of born-isolated nodes = total number of loops 
var(degree(gh, v = V(gh), mode = c("out"),loops = FALSE)) # calculate variance of out-degree 1.198598





#### Q-Q plot. Here we use the indegree matrix indegree_arr to draw the Q-Q plot ####################
indegree_mean <- rowMeans(indegree_arr ) # average indegree for comparison between theoretic result. the indegree_arr is generated from the above network formation process.
#write.csv(indegree_mean,file = "meanindegree_2Knodes20rounds.csv")
DATA <- indegree_mean[9:2000]
N         <- length(DATA);
PERCS     <- ((1:N)-0.5)/N; 
mu <- mean(degree(gh, mode="out")) # average outdegree
Pc <- (nNodes - sum(which_loop(gh)))/nNodes # proportion of born-collaborated nodes
Edge_cnt <- gsize(gh) # total number of edges

exp_IEB <- K*(p_1+(mu-1)*r_1)*nNodes/Edge_cnt # the exponential term
D <- Pc*mu*q_1/(p_1+(mu-1)*r_1) + (1-Pc)*(q-(mu-1)*r)/(p_1+(mu-1)*r_1) # construct the complex term in the theory
QUANTILES <- D/(1-PERCS)^exp_IEB - D; # use inverse function of CDF, calculate Quantile. The 

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
  geom_point(size = 2.5, colour = 'black') + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  #ggtitle('Quantile-Quantile Plot') + 
  #labs(subtitle = '(Comparison to exponential distribution with unit scale)') + 
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles')+ xlim(0,40) + ylim(0,40);

QQPLOT; # a few outliers are the initial small network
ggsave("/Users/HaoLI/Dropbox/Network/Network2019/QJ/images/QQplot_indegree_theory_simulation_IE.pdf", width = 10,   height = 8)
