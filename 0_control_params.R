# usage:
# change the 9 sets of parameters below and set the local julia path (sorry).
# then run from console with:

#    Rscript 0_control_params.R

# output is "results.csv".

seed   = 1440               # i like to do current time
nnet   =  800               # number of networks per scenario
ntaxa  =  c(5, 7)           # number of taxa per network
lambda =  c(0.1, 0.3, 1, 3) # speciation rate, in CUs
mu     =  c(0.1, 0.9)       # extinction rate, as a % of lambda
nu     =    0.5             # hybridization rate, as as % of lambda
M      =  c(0.5, 0.25)      # % lineage generative hybridizations.  Y always 0.25, H picks up the slack
d_0    =    0.5             # forbid hybridizations between lineages more than this % of lambda away
ngt    =  800               # number of gene trees per quartet
julia = "/home/john/julia-1.7.3/bin/julia"

# R will expand.grid the 9 parameter sets and ask SiPhyNetworks to simulate
# under every combination thereof ("scenario"), then pass the simulated networks to other
# julia + R scripts for further analysis.  the output is "results.csv": it has
# one row for every scenario, and summary statistics about each one.

# for example,
# as of 2022-08-19 there are 4 values for lambda, 2 for mu, 2
# for M, and 1 for all other parameters, so we have to test
# 4*2*2 = 16 scenarios.  results.csv will have 16 rows.

###

library(tictoc)

###

scenarios = expand.grid(seed=seed,
                        nnet=nnet,
                        ntaxa=ntaxa,
                        lambda=lambda,
                        mu=mu,
                        nu=nu,
                        M=M, 
                        d_0=d_0,
                        ngt=ngt)

scenarios$mu  = scenarios$mu  * scenarios$lambda # convert to CUs
scenarios$nu  = scenarios$nu  * scenarios$lambda
scenarios$d_0 = scenarios$d_0 * scenarios$lambda

scenarios$seed = scenarios$seed + 1:nrow(scenarios)

scenarios$Y = 0.25

scenarios = scenarios[,c("seed", "nnet", "ntaxa", "lambda", "mu", "nu", "M", "Y", "d_0", "ngt")]

###

try(system("rm results.csv"), TRUE)

for(i in 1:nrow(scenarios)) {
  
  command1 = "Rscript 1_sim_networks.R"
  command2 = paste(julia, "2_extract_quartet_subnetworks.jl", sep=" ")
  command3 = "Rscript 3_summarize_findings.R"
  
  # parameters for siphynetwork
  params1 = scenarios[i, 1:9]
  params1 = paste(unlist(params1), collapse=" ")

  # parameters for phylocoalsimulations
  params2 = scenarios[i, 10]
  params2 = paste(unlist(params2), collapse=" ")
  
  # parameters for final summary
  params3 = paste(params1, params2, sep=" ")
  
  print(scenarios[i,])

  command1 = paste(command1, params1, sep=" ")
  command2 = paste(command2, params2, sep=" ")
  command3 = paste(command3, params3, sep=" ")
  
  tic("1_sim_networks.R")
    system(command1)
  toc()
  
  tic("2_extract_quartet_subnetworks.R")
    system(command2)
  toc()
  
  tic("3_summarize_findings.R")
    system(command3)
  toc()
  
}
