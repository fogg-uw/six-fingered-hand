
summarize_findings = function(args) {
table_to_write = data.frame(
  seed   = args[01], # random seed
  nsim   = args[02], # number of networks to simulate
  ntaxa  = args[03], # for each simulation, stop when there are this many extant taxa
  lambda = args[04], # speciation rate, in CUs
  mu     = args[05], # extinction rate, in CUs
  nu     = args[06], # hybridization rate, in CUs
  M      = args[07], # proportion of hybridizations that are lineage generative 
  Y      = args[08], # proportion of hybridizations that are lineage degenerative
  d_0    = args[09], # in [0,1]: lineages cannot hybridize if further apart than d_0 * lambda 
  ngt    = args[10], # number of gene trees for phylocoalsims
  q = NA,
  A = NA,
  B = NA,
  C = NA,
  D = NA
)
if(!(quartets.csv %in% dir())) return(table_to_write)
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

table_to_write$q = nrow(quartets)
table_to_write$A = sum(quartets$A)
table_to_write$B = sum(quartets$B)
table_to_write$C = sum(quartets$C)
table_to_write$D = sum(quartets$D)

return(table_to_write)

}
