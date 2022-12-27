
summarize_findings = function(args) {
# initialize summary table: 1 row per scenario (=job)
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
  model  = args[10], # 0 = ssa, 1 = gsa
  ngt    = args[11], # baseline number of gene trees for phylocoalsims
  ngenes = NA, # actual number of simulated gene trees: may be higher than baseline
  rho = args[12],
  q   = NA, # total # of 4-taxon sets: ntaxa choose 4 x number of simulated nets
  q10 = NA, # number of tree-like 4-taxon sets
  q11 = NA, # number of 4-taxon sets with a 3_1 blob
  q12 = NA, # number of 4-taxon sets with a 3_2 blob
  q2  = NA, # number of 4-taxon sets displaying 2 splits
  q3  = NA, # number of 4-taxon sets displaying 3 splits
  A   = NA, # q12 + q2: when anomalies are theoretically possible
  A1  = NA, # q12: 1 split but 3_2 blob
  A2  = NA, # q2 : 2 splits
  B   = NA, # subset of A for which simulation was successful
  C   = NA,
  C1  = NA,
  C2  = NA,
  D0good = NA, # tree or 3_1 blob, and estimated as NOT anomalous
  D0anom = NA,
  D0ambi = NA,
  D1good = NA, # not anomalous
  D1anom = NA, # anomalous
  D1ambi = NA, # ambiguous: CI has anomaly threshold, even with more gene trees
  D2good = NA,
  D2anom = NA,
  D2ambi = NA,
  D0g_ex = NA, # same, but from exact calculation of expected qCFs
  D0a_ex = NA, # anomalous. No ambiguity with exact calculation
  D1g_ex = NA, # not anomalous
  D1a_ex = NA, # anomalous
  D2g_ex = NA,
  D2a_ex = NA
)
if(!("quartets.csv" %in% dir())) return(table_to_write)
quartets = read.csv("quartets.csv")
quartets$is32blob   = as.logical(quartets$is32blob  )
quartets$flag_class = as.logical(quartets$flag_class)
#nets = data.frame(sim_num = unique(quartets$sim_num))
#nrow(nets)

# first, break down quartets according to 1) # displayed trees and
# 2) if only 1 displayed tree, whether there is no 3 blob, or a 3-1 blob,
# or a 3-2 blob.

quartets$q10  = quartets$nsplit == 1 & quartets$num3blob_col == 0
quartets$q11  = quartets$nsplit == 1 & quartets$num3blob_col  > 0 & !quartets$is32blob
quartets$q12  = quartets$nsplit == 1 & quartets$num3blob_col  > 0 &  quartets$is32blob
quartets$q2   = quartets$nsplit == 2
quartets$q3   = quartets$nsplit == 3

# A: quartet has 2 displayed trees, or 1 displayed tree and a 3-2 blob.  analysis is desirable

quartets$A1   = quartets$q12
quartets$A2   = quartets$q2

quartets$A    = quartets$A1 | quartets$A2

# B: A is true and we have qCF data.  analysis is possible

quartets$B = quartets$A & quartets$split1 > -1

# C: B is true and (I or II)
# I. 1 displayed split and observed major qCF is less than 1/3 
# II. 2 displayed splits and observed minor qCF is not the smallest.
# i.e., analysis finds anomaly

quartets$smallestSplit                                                                        = 1
quartets$smallestSplit[quartets$split2 < quartets$split1]                                     = 2
quartets$smallestSplit[quartets$split3 < quartets$split2 & quartets$split3 < quartets$split1] = 3

quartets$C1 = quartets$A1 & quartets$B & quartets$split1         < 1/3
quartets$C2 = quartets$A2 & quartets$B & quartets$smallestSplit !=   3

quartets$C  = quartets$C1 | quartets$C2

# D: value of is_ambiguous (good, ambiguous, anomalous) separated by group.
# D0: class-1 networks theoretically NOT anomalous (tree or 3_1 blobs only)
# D1: class-1 networks with a 3_2 blob
# D2: class-2 network (with a 4-blob but only 2 displayed trees)

quartets$D0good = with(quartets, (q10 | q11) & (split1 > -1) & (is_anomalous == "good"))
quartets$D0anom = with(quartets, (q10 | q11) & (split1 > -1) & (is_anomalous == "anomalous"))
quartets$D0ambi = with(quartets, (q10 | q11) & (split1 > -1) & (is_anomalous == "ambiguous"))
quartets$D1good = with(quartets, A1 & B & (is_anomalous == "good"))
quartets$D1anom = with(quartets, A1 & B & (is_anomalous == "anomalous"))
quartets$D1ambi = with(quartets, A1 & B & (is_anomalous == "ambiguous"))
quartets$D2good = with(quartets, A2 & B & (is_anomalous == "good"))
quartets$D2anom = with(quartets, A2 & B & (is_anomalous == "anomalous"))
quartets$D2ambi = with(quartets, A2 & B & (is_anomalous == "ambiguous"))

quartets$D0g_ex = with(quartets, (q10 | q11) & (s1_exact >= 1/3))
quartets$D0a_ex = with(quartets, (q10 | q11) & (s1_exact <  1/3))
quartets$D1g_ex = with(quartets, A1 & (s1_exact >= 1/3))
quartets$D1a_ex = with(quartets, A1 & (s1_exact <  1/3))
quartets$D2g_ex = with(quartets, A2 & (s3_exact <= s1_exact) & (s3_exact <= s2_exact))
quartets$D2a_ex = with(quartets, A2 & ((s3_exact > s1_exact) | (s3_exact > s2_exact)))

quartets$n1 = as.integer(quartets$qCF_n * quartets$split1)
quartets$n2 = as.integer(quartets$qCF_n * quartets$split2)
quartets$n3 = as.integer(quartets$qCF_n * quartets$split3)

table_to_write$ngenes = mean(quartets$ngenes)
table_to_write$q    = nrow(quartets    )
table_to_write$q10  = sum( quartets$q10)
table_to_write$q11  = sum( quartets$q11)
table_to_write$q12  = sum( quartets$q12)
table_to_write$q2   = sum( quartets$q2 )
table_to_write$q3   = sum( quartets$q3 )
table_to_write$A    = sum( quartets$A  )
table_to_write$A1   = sum( quartets$A1 )
table_to_write$A2   = sum( quartets$A2 )
table_to_write$B    = sum( quartets$B  )
table_to_write$C    = sum( quartets$C  )
table_to_write$C1   = sum( quartets$C1 )
table_to_write$C2   = sum( quartets$C2 )
table_to_write$D0good = sum(quartets$D0good)
table_to_write$D0anom = sum(quartets$D0anom)
table_to_write$D0ambi = sum(quartets$D0ambi)
table_to_write$D1good = sum(quartets$D1good)
table_to_write$D1anom = sum(quartets$D1anom)
table_to_write$D1ambi = sum(quartets$D1ambi)
table_to_write$D2good = sum(quartets$D2good)
table_to_write$D2anom = sum(quartets$D2anom)
table_to_write$D2ambi = sum(quartets$D2ambi)
table_to_write$D0g_ex = sum(quartets$D0g_ex)
table_to_write$D0a_ex = sum(quartets$D0a_ex)
table_to_write$D1g_ex = sum(quartets$D1g_ex)
table_to_write$D1a_ex = sum(quartets$D1a_ex)
table_to_write$D2g_ex = sum(quartets$D2g_ex)
table_to_write$D2a_ex = sum(quartets$D2a_ex)

return(table_to_write)

}
