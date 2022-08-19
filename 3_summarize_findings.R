args = commandArgs(trailingOnly=TRUE)
print(args)
args = as.numeric(args)

summarize_findings = function(args) {
#setwd("/media/john/Phylo/research/2022-05-18 six-fingered hand")
quartets = read.csv("quartets.csv")
nets = data.frame(sim_num = unique(quartets$sim_num))
nrow(nets)

# starting with the quartet data, narrow in successively on different categories

sum(quartets$num3blob_col  > 0 & quartets$num4blob_col >  0)
sum(quartets$num3blob_col == 0 & quartets$num4blob_col >  0)
sum(quartets$num3blob_col  > 0 & quartets$num4blob_col == 0)
sum(quartets$num3blob_col == 0 & quartets$num4blob_col == 0)

# A: quartet has 3blob and no 4blob.  analysis is desirable
quartets$A = quartets$num3blob_col  > 0 & quartets$num4blob_col == 0

sum(quartets$A & quartets$split1 > -1)
sum(quartets$A & quartets$split1 <= -1)

# B: A is true and we have qCF data.  analysis is possible
quartets$B = quartets$A & quartets$split1 > -1

quartets$n1 = as.integer(quartets$qCF_n * quartets$split1)
quartets$n2 = as.integer(quartets$qCF_n * quartets$split2)
quartets$n3 = as.integer(quartets$qCF_n * quartets$split3)

#View(quartets[quartets$B, c("qCF_n", "n1", "n2", "n3")])

sum(quartets$B & quartets$split1 >= 1/3)
sum(quartets$B & quartets$split1 < 1/3 )

# C: B is true and observed major qCF is less than 1/3.  analysis finds
# "anomaly", in a loose sense
quartets$C = quartets$B & quartets$split1 < 1/3

quartets$test1_p = as.numeric(NA)
quartets$test2_p = as.numeric(NA)
for (i in 1:nrow(quartets)) {
  if(quartets$C[i]) {
    # are the two minor concordance factors equal?
    quartets$test1_p[i] = binom.test(x=quartets$n2[i],
                                     n=quartets$n2[i]+quartets$n3[i]
                                     )[['p.value']]
    # is the major concordance factor bigger than 1/3?
    quartets$test2_p[i] = binom.test(x=quartets$n1[i],
                                     n=quartets$qCF_n[i],
                                     p=1/3,
                                     alternative='less'
                                     )[['p.value']]
  }
}

alpha = 0.05
numtests = sum(quartets$C)
# D: C is true and population major CF tested to be < 1/3.  analysis finds
# "anomaly" in a stricter sense
quartets$D = quartets$C & quartets$test2_p < alpha / numtests

table_to_write = data.frame(
  seed   = args[1], # random seed
  nsim   = args[2], # number of networks to simulate
  ntaxa  = args[3], # for each simulation, stop when there are this many extant taxa
  lambda = args[4], # speciation rate, in CUs
  mu     = args[5], # extinction rate, in CUs
  nu     = args[6], # hybridization rate, in CUs
  M      = args[7], # proportion of hybridizations that are lineage generative 
  Y      = args[8], # proportion of hybridizations that are lineage degenerative
  d_0    = args[9], # in [0,1]: lineages cannot hybridize if further apart than d_0 * lambda 
  q = nrow(quartets),
  A = sum(quartets$A),
  B = sum(quartets$B),
  C = sum(quartets$C),
  D = sum(quartets$D)
)
write.col.names.YN = !("results.csv" %in% dir())
write.table(
  x=table_to_write,
  file="results.csv",
  append=TRUE,
  sep=",",
  row.names=F,
  col.names=write.col.names.YN)
}

summarize_findings(args)
