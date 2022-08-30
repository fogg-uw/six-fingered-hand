## usage:
## Rscript 1_sim_networks.R [seed] [nnet] ... [d_0]
## see below for all required arguments.

#cat("1_sim_networks.R\n")
args = commandArgs(trailingOnly=TRUE)
#print(args)
args = as.numeric(args)

seed   = args[1] # random seed
nnet   = args[2] # number of networks to simulate
ntaxa  = args[3] # for each simulation, stop when there are this many extant taxa
lambda = args[4] # speciation rate, in CUs
mu     = args[5] # extinction rate, in CUs
nu     = args[6] # hybridization rate, in CUs
M      = args[7] # proportion of hybridizations that are lineage generative 
Y      = args[8] # proportion of hybridizations that are lineage degenerative
d_0    = args[9] # lineages cannot hybridize if further apart than d_0

###

#library(tictoc)
library(SiPhyNetwork)
library(ape)

###

hybrid_proportions <-c(  M,      ##Lineage Generative
                           Y,    ##Lineage Degenerative
                        -M-Y+1 ) ##Lineage Neutral
 
hybrid_success_prob <- make.stepwise(
  probs     = c(  1,   0),
  distances = c(d_0, Inf)
)

inheritance.fxn <- make.beta.draw(1,1) # beta(1,1) distribution (i.e. unif(0,1))

# sometimes networks go extinct, so i do a while-loop to keep running
# SiPhyNetwork until i get the number of extant networks i want.

#tic("starting ssa_nets")
ssa_nets_full = list()
while(length(ssa_nets_full) < nnet) {
  numbsim = nnet - length(ssa_nets_full)
  ssa_nets <- sim.bdh.taxa.ssa(
    n             = ntaxa,
    numbsim       = numbsim,
    lambda        = lambda,
    mu            = mu,
    nu            = nu,
    hybprops      = hybrid_proportions,
    hyb.inher.fxn = inheritance.fxn,
    complete      = FALSE # do not return extinct taxa
  )
  isExtinct = sapply(X=ssa_nets, FUN = function(x) identical(x,0))
  extantNets = ssa_nets[!isExtinct]
  ssa_nets_full = c(ssa_nets_full, extantNets)
}
#toc()

outputdir = "SiPhyNetwork_output"
unlink(outputdir, recursive=TRUE) # delete all existing output!
dir.create(outputdir)
NL = nchar(as.character(nnet))
for(i in 1:nnet) {
  j = as.character(i)
  while(nchar(j) < NL) {
    j = paste0(0, j)
  }
  file = paste0('sim', j, '.tree')
  file = file.path(outputdir, file)
  try(write.net(ssa_nets_full[[i]], file = file), TRUE)
}
