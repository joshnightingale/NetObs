#### Script to test for successful detections of top sites based on sampling individual bird movements

### Contents:
## 1. Setup
# load packages etc
## 2. Make network structure
# simulate a scale-free network with preferential attachment
## 3. Simulate "true" network
# movements of a large number of individuals through the network
# is our latent reality we are trying to sample
# NOTE! Runs quite slowly with large sample size (~30 mins) 
## 4. Effect of sample size with perfect observations (~tracking)
# Build a network from a small subsample of the simulated individuals
# How many of the top sites (using betweenness) do we detect?
## 5. Effect of sample size with imperfect observations (~ringing)
# As above, build a network from a small subsample of the simulated individuals
# BUT each site has observation probability (fixed per site within each rep)
# which determines whether that obs of that individual is included.
# How many of the top sites (using betweenness) do we detect?


#### 1. Setup ####
library(igraph)
library(magrittr)

set.seed(69)


##### 2. Make network structure #####
Nsite <- 1000

# sample random network
g <- sample_pa(Nsite, power = 1, 
               out.seq = rpois(Nsite, 1), # number of links to add per node
               zero.appeal = 0.1, # relative attractiveness of sites with no links
               directed = F)



#### 3. Simulate true network ####

### individual movements 
Nbird <- 100000 # number of birds for "true" network: a big number
Nstep <- 10 # (max) number of movements each bird makes

# # (commented out because slow)
# movements <- matrix(NA, nrow=Nbird, ncol=Nstep)
# ## fill starting point for each bird
# # this code allows starting at any connected node with equal probability
# movements[,1] <- sample(1:Nsite, Nbird, replace=T, prob=( degree(g)>1)  )
# # the code below has starting probability relative to node degree (more realistic?)
# #, prob=( degree(g)/sum( degree(g)) ) )


## loop through birds 
# ii <- 1
# for(ii in 1:Nbird) {
#   ij <- 2
#   for(ij in 2: Nstep){
#     # find adjacent edges
#     ae <- neighbors(g, movements[ii, ij-1] )
#     if(length(ae) > 0) {
#       movements[ii, ij] <- sample(ae, 1, replace = F)
#     } else { 
#       movements[ii, ij] <- movements[ii, ij-1]} # close if(length(ae))
#     
#   } # close for(ij)
# } # close for(ii)

# load previous movement simulation
movements <- read.csv("movements.csv") %>% as.matrix


### function to define network from movements
obs_net <- function(Nsite, movements) {
  pb <- txtProgressBar(style=3)
  
  edge_list <- vector("list", nrow(movements))
  
  # add each bird's sequential movements to net
  ik <- 1
  for(ik in 1:nrow(movements)){
    ee <- na.omit(movements[ik,][!duplicated(movements[ik,])])
    edge_list[[ik]] <- cbind(ee[1:(length(ee)-1)], ee[2:(length(ee))])
   
    setTxtProgressBar(pb, ik/nrow(movements))
  } # close for(ik)
  
  edge_all <- do.call("rbind", edge_list)
  og <- make_empty_graph(n=Nsite)
  og <- add_edges(og, edge_all)
  
  return(og)
  close(pb)
}
x <- obs_net(Nsite, movements[1:10,]) # test

### make true graph
# tg <- obs_net(Nsite, movements)
saveRDS(tg, "true_graph.rds")
tg <- readRDS("true_graph.rds") # load previous

dtg <- degree(tg)
btg <- betweenness(tg)




#### 4. Effect of sample size - perfect observations ####
Nsamp <- c(10, 100, 323, 1000)
Nrep <- 100

perfect <- matrix(NA, nrow=length(Nsamp), ncol=Nrep)

ii <- 1
for (ii in 1:length(Nsamp)) {
  res <- rep(NA, Nrep)
  
  ir <- 1
  for(ir in 1:Nrep){
    samp <- movements[sample(1:Nbird, Nsamp[ii], replace = F),]
    
    sg <- obs_net(Nsite, samp)
    
    bsg <- betweenness(sg)
    dsg <- degree(sg)
    
    top50tru <- which(rank(btg) %in% 951:1000)
    top50obs <- which(rank(bsg) %in% 951:1000)
    
    perfect[ii,ir] <- sum(top50tru %in% top50obs) / 50 * 100
  }
  
  print(round(ii/length(Nsamp)*100,1)) # progress bar
}
rowMeans(perfect)
write.csv(perfect, "perfect.csv", row.names = F)

#### 5. Sample sizes with imperfect observations ####
Nsamp <- c(100, 1000, 10000)
Nrep <- 100

imperfect <- matrix(NA, nrow=length(Nsamp), ncol=Nrep)

ii <- 1
for (ii in 1:length(Nsamp)) {
  res <- rep(NA, Nrep)
  
  ir <- 1
  for(ir in 1:Nrep){
    # random observation probabilities at each site
    obs_probs <- runif(Nsite, min=0.1, max=0.9)
    
    # sample some individuals
    samp <- movements[sample(1:Nbird, Nsamp[ii], replace = F),]
    
    # simulate observation process
    samp_obs <- matrix(NA, ncol=ncol(samp), nrow=nrow(samp))
    iis <- 1
    for (iis in 1:length(samp_obs)) {
      # samp is observed depending on site probability
      samp_obs[iis] <- samp[iis] * rbinom(1, 1, prob=obs_probs[samp[iis]])
    }
    samp_obs[samp_obs==0] <- NA
    
    # only keep birds observed in at least 2 sites
    samp_obs <- samp_obs[which(apply(samp_obs, 1, \(x) {length(unique(na.omit(x)))}) >= 2),]
    
    sg <- obs_net(Nsite, samp_obs)
    
    bsg <- betweenness(sg)
    dsg <- degree(sg)
    
    top50tru <- which(rank(btg) %in% 951:1000)
    top50obs <- which(rank(bsg) %in% 951:1000)
    
    imperfect[ii,ir] <- sum(top50tru %in% top50obs) / 50 * 100
  }
  
  print(ii/length(Nsamp))
}
rowMeans(imperfect)
write.csv(imperfect, "imperfect.csv", row.names = F)