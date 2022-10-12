# usage:
# change the 10 sets of parameters below and set the local julia and R paths (sorry).
# then run from console with:

#    [Rscript path] 0_control_params.R [r XOR j XOR nothing]

# arg = r if you just want it to sim networks in R, arg = j if you just want to
# analyze the nets in julia, arg = blank if you want both / don't know what
# i'm talking about.

# output is "results.csv".

seed   =  9321
nnet   =  200                # number of networks per scenario
ntaxa  =  c(4, 6, 8)         # number of taxa per network
lambda =  c(0.1, 0.3, 1, 3)  # speciation rate, in CUs
mu     =  c(0.1, 0.9)        # extinction rate, as a % of lambda
nu     =    0.5              # hybridization rate, as as % of lambda
MHY    =  c("M", "H")        # type of hybridization that is dominant
d_0    =  c(0.1, 0.3, 0.6)*2 # forbid hybridizations between lineages more than this % of 1/lambda away
model  =  1                  # ssa = 0, gsa = 1.  see hartmann wong stadler 2010
ngt    =  200                # number of gene trees per quartet

ncores_julia = 8
ncores_R     = 8 # just for looping over scenarios.  not for looping over nets

simulate_nets_in_R = FALSE
analyze_nets_in_julia = FALSE

julia  = "/u/f/o/fogg/julia-1.8.0/bin/julia"
R      = "Rscript"

#on john's machine: julia = "/home/john/julia-1.7.3/bin/julia"
#on franklin00: julia = "/u/f/o/fogg/julia-1.8.0/bin/julia"

timeout = "4h" # i don't know if i've ever seen it work in a time >13min but <20min
delete1 = FALSE # whether to simulate up to N+1 taxa, then delete 1 later

arg = commandArgs(trailingOnly=TRUE)
if(length(arg) > 1) error("too many args")
if(is.na(arg[1])) arg == 'rj'
if(arg=='r') simulate_nets_in_R = TRUE # if false, you need to provide the nets yourself
if(arg=='j') analyze_nets_in_julia = TRUE
if(arg=='rj') {
  simulate_nets_in_R = TRUE
  analyze_nets_in_julia = TRUE
}

# R will expand.grid the 10 parameter sets and ask SiPhyNetworks to simulate
# under every combination thereof ("scenario"), then pass the simulated networks to other
# julia + R scripts for further analysis.  the output is "results.csv": it has
# one row for every scenario, and summary statistics about each one.

# for example,
# as of 2022-08-19 there are 4 values for lambda, 2 for mu, 2
# for M, and 1 for all other parameters, so we have to test
# 4*2*2 = 16 scenarios.  results.csv will have 16 rows.

###

library(parallel)
library(tictoc, lib.loc=.libPaths()[1])

###

tic("main job")
scenarios = expand.grid(seed=seed,
                        nnet=nnet,
                        ntaxa=ntaxa,
                        lambda=lambda,
                        mu=mu,
                        nu=nu,
                        MHY=MHY, 
                        d_0=d_0,
                        model=model,
                        ngt=ngt)

scenarios$mu  = scenarios$mu  * scenarios$lambda # convert to CUs
scenarios$nu  = scenarios$nu  * scenarios$lambda
scenarios$d_0 = scenarios$d_0 / scenarios$lambda

scenarios$seed = min(scenarios$seed) - 1 + 1:nrow(scenarios)

scenarios$M = 1/4 + (1/4)*as.numeric(scenarios$MHY=="M")
scenarios$Y = 1/4 + (1/4)*as.numeric(scenarios$MHY=="Y")

scenarios =
  scenarios[
    ,
    c("seed",
      "nnet",
      "ntaxa",
      "lambda",
      "mu",
      "nu",
      "M",
      "Y",
      "d_0",
      "model",
      "ngt"
      )
    ]

###

startingdir = getwd()
print(startingdir)
try(system("rm results.csv"), TRUE)

script1 = paste0('\"', file.path(startingdir, "1_sim_networks.R"),                 '\"')
script2 = paste0('\"', file.path(startingdir, "2_extract_quartet_subnetworks.jl"), '\"')
script3 = paste0('\"', file.path(startingdir, "3_summarize_findings.R"),           '\"')

parallel_job = function(i) {
  
  jobdir = paste0("job", i)
  tic(jobdir)
  cat(paste('start', jobdir, '\n'))
  if(simulate_nets_in_R) { # then we need to clean the directory.  o/w leave alone
    unlink(jobdir, recursive = TRUE)
    dir.create(jobdir)
  }
  setwd(jobdir)
  
  # parameters for siphynetwork
  params1 = scenarios[i, 1:10]
  if(delete1) params1["ntaxa"] = params1["ntaxa"] + 1
  params1 = paste(unlist(params1), collapse=" ")
  
  # parameters for phylocoalsimulations
  params2 = scenarios[i, c(1,11)]
  params2 = paste(unlist(params2), collapse=" ")
  params2 = paste(params2, as.numeric(delete1), collapse=" ")
  
  command1 = paste("timeout", timeout, R,                                script1, params1, sep=" ")
  command2 = paste("timeout", timeout, julia, "--threads", ncores_julia, script2, params2, sep=" ")
  
  if(simulate_nets_in_R) {
    cat(paste(jobdir, command1, "\n"))
    system(command1)
  }
  if(analyze_nets_in_julia) {
    cat(paste(jobdir, command2, "\n"))
    system(command2)
  }
  
  setwd(startingdir)
  cat('finish', jobdir, '\n')
  time = toc()
  return(time$toc - time$tic)
}

serial_job = function(i) {
  
  source(noquote(script3))
  
  setwd(startingdir)
  jobdir = paste0("job", i)
  setwd(jobdir)
  
  # parameters for final summary
  params3 = scenarios[i, 1:11]
  params3 = as.numeric(params3)
  

  table_to_write = summarize_findings(params3)
  
  setwd(startingdir)
  #unlink(jobdir, recursive=TRUE) # i want to leave these for followup
  return(table_to_write)
  
}


if(TRUE)                  results_parallel = mclapply(X=1:nrow(scenarios), FUN=parallel_job, mc.cores = ncores_R)
if(analyze_nets_in_julia) results_serial   =   lapply(X=1:nrow(scenarios), FUN=  serial_job                     )


results = results_serial[[1]]
for(i in 2:length(results_serial)) {
  results = rbind(results, results_serial[[i]])
}

results$time = unlist(results_parallel)

write.csv(results, "results.csv")
toc()
