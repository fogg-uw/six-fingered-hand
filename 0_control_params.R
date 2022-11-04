# usage:
# change the sets of parameters below, set the local julia and R paths
# then run from console with:

#    [Rscript path] 0_control_params.R [1 XOR 2 XOR 3 XOR nothing]

# arg = 1 if you just want it to sim networks in R, arg = 2 if you just want to
# analyze the nets in julia, arg = 3 if you just want to summarize the results
# in R, arg = blank if you want all three / don't know what
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

ncores_julia = 1
ncores_R     = 64 # for looping over scenarios (not over networks)
# 64 dual cores on franklinxx machines: ncores_julia * ncores_R = 64 would make sense

julia  = "julia" # assuming it's in your path
R      = "Rscript"

# on john's machine: julia = "/home/john/julia-1.7.3/bin/julia"
# on franklin00: julia = "/u/f/o/fogg/julia-1.8.0/bin/julia"
#        or shared path:  /s/julia-1.8.1/bin/julia

timeout = "4h" # max seen: 38 min (after increasing # gene trees to remove ambiguous tests)
delete1 = FALSE # whether to simulate up to N+1 taxa, then delete 1 later

# parse arguments to decide which steps of the analysis to run
# if you don't run step 1, you need to manually put the inputs for step 2
# where the program expects them; similarly if you don't run step 3
run_step = c(FALSE, FALSE, FALSE)
arg = commandArgs(trailingOnly=TRUE)
if (length(arg) > 1) stop("too many arguments. Use no arguments to run all steps")
if (length(arg) == 0){
  run_step = c(TRUE, TRUE, TRUE)
} else {
  if(arg!='1' & arg!='2' & arg!='3') stop(paste("unacceptable arg:",arg))
  if(is.na(arg) | arg=='1') run_step[1] = TRUE
  if(is.na(arg) | arg=='2') run_step[2] = TRUE
  if(is.na(arg) | arg=='3') run_step[3] = TRUE
}

# R will expand.grid the parameter sets and ask SiPhyNetworks to simulate
# under every combination thereof ("scenario"), then pass the simulated networks to other
# julia + R scripts for further analysis.  the output is "results.csv": it has
# one row for every scenario, and summary statistics about each one.

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
num_parameters = ncol(scenarios) # 10, to avoid hard-coding 10 or 11 below

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
if(run_step[3]) try(system("rm results.csv"), TRUE)

script1 = paste0('\"', file.path(startingdir, "1_sim_networks.R"),                 '\"')
script2 = paste0('\"', file.path(startingdir, "2_extract_quartet_subnetworks.jl"), '\"')
script3 = paste0(      file.path(startingdir, "3_summarize_findings.R")                )


parallel_job = function(i) {
  
  jobdir = paste0("job", i)
  tic(jobdir)
  cat(paste('start', jobdir, '\n'))
  if(run_step[1]) { # then we need to clean the directory.  o/w leave alone
    unlink(jobdir, recursive = TRUE)
    dir.create(jobdir)
  }
  setwd(jobdir)
  
  # parameters for siphynetwork
  params1 = scenarios[i, 1:num_parameters]
  if(delete1) params1["ntaxa"] = params1["ntaxa"] + 1
  params1 = paste(unlist(params1), collapse=" ")
  
  # parameters for phylocoalsimulations
  params2 = scenarios[i, c(1,num_parameters+1)]
  params2 = paste(unlist(params2), collapse=" ")
  params2 = paste(params2, as.numeric(delete1), collapse=" ")
  
  command1 = paste("timeout", timeout, R,                                script1, params1, sep=" ")
  command2 = paste("timeout", timeout, julia, "--project --threads", ncores_julia, script2, params2, sep=" ")
  
  if(run_step[1]) {
    cat(paste(jobdir, command1, "\n"))
    system(command1)
  }
  if(run_step[2]) {
    cat(paste(jobdir, command2, "\n"))
    system(command2)
  }
  
  setwd(startingdir)
  cat('finish', jobdir, '\n')
  time = toc()
  return(time$toc - time$tic)
}

serial_job = function(i) {
  
  source(script3)
  
  setwd(startingdir)
  jobdir = paste0("job", i)
  setwd(jobdir)
  
  # parameters for final summary
  params3 = scenarios[i, 1:(num_parameters+1)]
  params3 = as.numeric(params3)


  table_to_write = summarize_findings(params3)
  
  setwd(startingdir)
  #unlink(jobdir, recursive=TRUE) # i want to leave these for followup
  return(table_to_write)
  
}


if(run_step[1] | run_step[2]) results_parallel = mclapply(X=1:nrow(scenarios), FUN=parallel_job, mc.cores = ncores_R)
if(run_step[3]              ) results_serial   =   lapply(X=1:nrow(scenarios), FUN=  serial_job                     )

if(run_step[3]) {
  results = results_serial[[1]]
  for(i in 2:length(results_serial)) {
    results = rbind(results, results_serial[[i]])
  }
  if(run_step[2]) results$time = unlist(results_parallel)
  write.csv(results, "results.csv", quote=F)
}

toc()
