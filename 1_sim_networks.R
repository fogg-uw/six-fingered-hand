library(SiPhyNetwork)
library(ape)

###

set.seed(1347) # current time
#numbsim = 100
numbsim = 400
#age = 2 # number of time units to simulate to... more than 2 can cause trouble
n = 7 # number of taxa to simulate to

###

inheritance.fxn <- make.beta.draw(1,1) # beta(1,1) distribution (i.e. unif(0,1))
hybrid_proportions <-c(0.5,  ##Lineage Generative
                       0.25, ##Lineage Degenerative
                       0.25) ##Lineage Neutral
ssa_nets <- sim.bdh.taxa.ssa(
  #age = age,
  n = 7,
  numbsim = numbsim,
  lambda=1, # speciation rate
  mu=0.2, # extinction rate
  nu=0.25, # hybridization rate,
  hybprops = hybrid_proportions,
  hyb.inher.fxn = inheritance.fxn,
  complete=FALSE
)
#plot(ssa_nets[[9]]) # oh goodness
#for numbsim=1, age=2, age=4 happen "fast" (before i can hit stopwatch)
#age=8 takes more than 60sec
#took more than 60s for age=20, numbsim=1

outputdir = "SiPhyNetwork_output"
unlink(outputdir, recursive=TRUE)
dir.create(outputdir)
for(i in 1:numbsim) {
  file = paste0('sim', i, '.tree')
  file = file.path(outputdir, file)
  try(write.net(ssa_nets[[i]],file = file))
  #try(ape::write.evonet(ssa_nets[[i]],file = file))
}
