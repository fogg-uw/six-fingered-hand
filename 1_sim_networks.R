library(tictoc)
library(SiPhyNetwork)
library(ape)

###

set.seed(1426) # current time
#numbsim = 100
numbsim = 1600
#age = 2 # number of time units to simulate to... more than 2 can cause trouble
n = 7 # number of taxa to simulate to

###

inheritance.fxn <- make.beta.draw(1,1) # beta(1,1) distribution (i.e. unif(0,1))
hybrid_proportions <-c(0,  ##Lineage Generative
                       1, ##Lineage Degenerative
                       0) ##Lineage Neutral
#f3<-make.stepwise(probs = c(1,0),distances = c(0.03125/2,Inf))

tic("starting ssa_nets")
ssa_nets <- sim.bdh.taxa.ssa(
  #age = age,
  n = n,
  numbsim = numbsim,
  lambda=0.03125, # speciation rate
  mu=0, # extinction rate
  nu=0.00625, # hybridization rate,
  hybprops = hybrid_proportions,
  hyb.inher.fxn = inheritance.fxn,
  complete=FALSE
)
toc("finished ssa_nets")
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
