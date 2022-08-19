# usage:
# change the 8 sets of parameters below.
# R will expand.grid them and ask SiPhyNetworks to simulate
# under every combination thereof.

# for example,
# as of 2022-08-19 there are 4 values for lambda, 2 for mu, 2
# for M, and 1 for all other parameters, so we have to test
# 4*2*2 = 16 scenarios.

seed   = 1440               # current time
nsim   =  100
ntaxa  =    5
lambda =  c(0.1, 0.3, 1, 3)
mu     =  c(0.1, 0.9)       # multiply by lambda
nu     =    0.5             # multiply by lambda
M      =  c(0.5, 0.25)      # Y always 0.25, H picks up the slack
d_0    =    0.5             # multiply by lambda

###

scenarios = expand.grid(seed=seed,
                        nsim=nsim,
                        ntaxa=ntaxa,
                        lambda=lambda,
                        mu=mu,
                        nu=nu,
                        M=M, 
                        d_0=d_0)

scenarios$mu  = scenarios$mu  * scenarios$lambda
scenarios$nu  = scenarios$nu  * scenarios$lambda
scenarios$d_0 = scenarios$d_0 * scenarios$lambda

scenarios$seed = scenarios$seed + 1:nrow(scenarios)

scenarios$Y = 0.25

scenarios = scenarios[,c("seed", "nsim", "ntaxa", "lambda", "mu", "nu", "M", "Y", "d_0")]

###

try(system("rm results.csv"), TRUE)

for(i in 1:nrow(scenarios)) {
  
  command1 = "Rscript 1_sim_networks.R"
  command2 = "julia 2_extract_quartet_subnetworks.jl"
  command3 = "Rscript 3_summarize_findings.R"
  
  params = scenarios[i,]
  params = paste(unlist(params), collapse=" ")

  command1 = paste(command1, params, sep=" ")
  #command2 = paste(command2,     "", sep=" ")
  command3 = paste(command3, params, sep=" ")
  
  system(command1)
  system(command2)
  system(command3)
  
}
