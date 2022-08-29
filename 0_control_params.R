# usage:
# change the 9 sets of parameters below and set the local julia and R paths (sorry).
# then run from console with:

#    [Rscript path] 0_control_params.R

# output is "results.csv".

seed   = 1440               # i like to do current time
#nnet   =  800               # number of networks per scenario
nnet    = 2
ntaxa  =  c(5, 7)           # number of taxa per network
lambda =  c(0.1, 0.3, 1, 3) # speciation rate, in CUs
mu     =  c(0.1, 0.9)       # extinction rate, as a % of lambda
nu     =  c(0.2, 0.5)       # hybridization rate, as as % of lambda
M      =  c(0.5, 0.25)      # % lineage generative hybridizations.  Y always 0.25, H picks up the slack
d_0    =  c(0.1, 0.3, 0.6)  # forbid hybridizations between lineages more than this % of 1/lambda away
#ngt    =  800               # number of gene trees per quartet
ngt    = 2

julia  = "/u/f/o/fogg/julia-1.8.0/bin/julia"
R      = "Rscript"

#on john's machine: julia = "/home/john/julia-1.7.3/bin/julia"

# R will expand.grid the 9 parameter sets and ask SiPhyNetworks to simulate
# under every combination thereof ("scenario"), then pass the simulated networks to other
# julia + R scripts for further analysis.  the output is "results.csv": it has
# one row for every scenario, and summary statistics about each one.

# for example,
# as of 2022-08-19 there are 4 values for lambda, 2 for mu, 2
# for M, and 1 for all other parameters, so we have to test
# 4*2*2 = 16 scenarios.  results.csv will have 16 rows.

###

library(parallel)
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
scenarios$d_0 = scenarios$d_0 / scenarios$lambda

scenarios$seed = scenarios$seed + 1:nrow(scenarios)

scenarios$Y = 0.25

scenarios = scenarios[,c("seed", "nnet", "ntaxa", "lambda", "mu", "nu", "M", "Y", "d_0", "ngt")]

###

startingdir = getwd()
try(system("rm results.csv"), TRUE)

script1 = file.path(startingdir, "1_sim_networks.R")
script2 = file.path(startingdir, "2_extract_quartet_subnetworks.jl")
script3 = file.path(startingdir, "3_summarize_findings.R")

parallel_job = function(i) {
  
  jobdir = paste0("job", i)
  unlink(jobdir)
  dir.create(jobdir)
  setwd(jobdir)
  
  command1 = paste(R,     script1, sep=" ")
  command2 = paste(julia, script2, sep=" ")
  
  # parameters for siphynetwork
  params1 = scenarios[i, 1:9]
  params1 = paste(unlist(params1), collapse=" ")
  
  # parameters for phylocoalsimulations
  params2 = scenarios[i, 10]
  params2 = paste(unlist(params2), collapse=" ")
  
  command1 = paste(command1, params1, sep=" ")
  command2 = paste(command2, params2, sep=" ")
  
  system(command1)
  system(command2)
  
  setwd(startingdir)
  return(i)
}

serial_job = function(i) {
  
  source(script3)
  
  setwd(startingdir)
  jobdir = paste0("job", i)
  setwd(jobdir)
  
  #command3 = paste(R,     "3_summarize_findings.R",           sep=" ")
  
  # parameters for final summary
  params3 = scenarios[i, 1:10]
  #params3 = paste(unlist(params2), collapse=" ")
  params3 = as.numeric(params3)
  
  #command3 = paste(command3, params3, sep=" ")
  #system(command3)
  table_to_write = summarize_findings(params3)
  
  setwd(startingdir)
  unlink(jobdir, recursive=TRUE)
  return(table_to_write)
  
}

numCores = detectCores()-1 #let's not be greedy
results_parallel = mclapply(X=1:nrow(scenarios), FUN=parallel_job, mc.cores = numCores)
results_serial   =   lapply(X=1:nrow(scenarios), FUN=  serial_job                     )

results = results_serial[[1]]
for(i in 2:length(results_serial)) {
  results = rbind(results, results_serial[[i]])
}
