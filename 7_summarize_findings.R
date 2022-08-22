setwd("/media/john/Phylo/research/2022-05-18 six-fingered hand/six-fingered-hand")
quartets = read.csv("quartets.csv")
source("/home/john/Documents/projects/code/checkUnits3.R")
nets = data.frame(sim_num = unique(quartets$sim_num))
nrow(nets)

# starting with the quartet data, narrow in successively on different categories

sum(quartets$num3blob_col  > 0 & quartets$num4blob_col >  0)
sum(quartets$num3blob_col == 0 & quartets$num4blob_col >  0)
sum(quartets$num3blob_col  > 0 & quartets$num4blob_col == 0)
sum(quartets$num3blob_col == 0 & quartets$num4blob_col == 0)
sum(quartets$num2cycle_col > 0)
sum(quartets$num5blob_col > 0)

# A: quartet has 3blob and no 4blob.  analysis is desirable
quartets$A = quartets$num3blob_col  > 0 & quartets$num4blob_col == 0

sum(quartets$A & quartets$split1 > -1)
sum(quartets$A & quartets$split1 <= -1)

# B: A is true and we have qCF data.  analysis is possible
quartets$B = quartets$A & quartets$split1 > -1

quartets$n1 = as.integer(quartets$qCF_n * quartets$split1)
quartets$n2 = as.integer(quartets$qCF_n * quartets$split2)
quartets$n3 = as.integer(quartets$qCF_n * quartets$split3)

View(quartets[quartets$B, c("qCF_n", "n1", "n2", "n3")])

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


nets = smoosh(
  fromData = quartets,
  fromColumn = "A",
  byColumns = "sim_num",
  toData = nets,
  aggFun = function(x) max(as.numeric(x)), # or
  newName = 'anyquartet3blob'
)
sum(nets$anyquartet3blob)

nets = smoosh(
  fromData = quartets,
  fromColumn = "B",
  byColumns = "sim_num",
  toData = nets,
  aggFun = function(x) max(as.numeric(x)), # or
  newName = 'noHLproblem'
)
sum(nets$noHLproblem)

nets = smoosh(
  fromData = quartets,
  fromColumn = "C",
  byColumns = "sim_num",
  toData = nets,
  aggFun = function(x) max(as.numeric(x)), # or
  newName = 'anypathologicalquartet'
)
sum(nets$anypathologicalquartet)

nets = smoosh(
  fromData = quartets,
  fromColumn = "D",
  byColumns = "sim_num",
  toData = nets,
  aggFun = function(x) max(as.numeric(x)), # or
  newName = 'anytrulypathologicalquartet'
)
sum(nets$anytrulypathologicalquartet)

#transform = matrix(c(1/sqrt(3),0,0,-1/sqrt(3),0,0,0,1,0), ncol=3) # "project" simplex in R^3 into R^2
transform = matrix(c(0,1,0,-1/sqrt(3),0,0,1/sqrt(3),0,0), ncol=3) # "project" simplex in R^3 into R^2

good_splits_for_simplex = t(transform %*% t(quartets[quartets$B & !quartets$C,              c('split1', 'split2', 'split3')]))[,c(1,2)]
ok_splits_for_simplex   = t(transform %*% t(quartets[quartets$B &  quartets$C,              c('split1', 'split2', 'split3')]))[,c(1,2)]
bad_splits_for_simplex  = t(transform %*% t(quartets[quartets$B &  quartets$C & quartets$D, c('split1', 'split2', 'split3')]))[,c(1,2)]


shoulder1 = cbind(x=seq(from=0, to=1, by=0.01), y=seq(from=1, to=0, by=-0.01), z=0)
shoulder2 = cbind(x=seq(from=1, to=0, by=-0.01), y=0, z=seq(from=0, to=1, by=0.01))
shoulder3 = cbind(x=0, y=seq(from=0, to=1, by=0.01), z=seq(from=1, to=0, by=-0.01))
shoulders = rbind(shoulder1, shoulder2, shoulder3)
shoulders_2 = t(transform %*% t(shoulders))
shoulders_2 = shoulders_2[,c(1,2)]

arm_major = seq(from=1,    to=0.33, by=-0.01)
hand_major  = seq(from=0.33, to=0,    by=-0.01)

arm_minor = seq(from=0,     to=0.335, by=0.005)
hand_minor  = seq(from=0.335, to=0.5,   by=0.005)

arm1 = cbind(x=arm_major, y=arm_minor, z=arm_minor)
arm2 = cbind(x=arm_minor, y=arm_major, z=arm_minor)
arm3 = cbind(x=arm_minor, y=arm_minor, z=arm_major)
arms = rbind(arm1, arm2, arm3)
arms_2 = t(transform %*% t(arms))[,c(1,2)]

hand1 = cbind(x=hand_major, y=hand_minor, z=hand_minor)
hand2 = cbind(x=hand_minor, y=hand_major, z=hand_minor)
hand3 = cbind(x=hand_minor, y=hand_minor, z=hand_major)
hands = rbind(hand1, hand2, hand3)
hands_2 = t(transform %*% t(hands))[,c(1,2)]

data_to_plot = rbind(shoulders_2, arms_2, hands_2, good_splits_for_simplex, ok_splits_for_simplex, bad_splits_for_simplex)
data_to_plot = data.frame(data_to_plot)
data_to_plot$type = 
  c(rep("1. shoulder", nrow(shoulders)),
    rep("2. arm", nrow(arms)),
    rep("3. hand", nrow(hands)),
    rep("4. major CF > 1/3", nrow(good_splits_for_simplex)),
    rep("5. major CF < 1/3 but retain null", nrow(ok_splits_for_simplex)),
    rep("6. major CF < 1/3 and reject null", nrow(bad_splits_for_simplex))
    )

library(ggplot2)
myplot =
  ggplot(
    data=data_to_plot,
    aes(x=X1, y=X2, color=type)
    ) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(
    values=c('black', "grey50", "grey75", "green", "yellow", "red")
    ) +
  theme(legend.position='bottom') +
  theme(legend.title = element_blank()) +
  xlab("") +
  ylab("")
myplot

quartets$split1_s = as.numeric()
quartets$split2_s = as.numeric()
for(i in 1:nrow(quartets)) {
  split = c(quartets$split1[i], quartets$split2[i], quartets$split3[i])
  split_s = transform %*% split
  split_s = split_s[c(1,2)]
}

a = c(1/sqrt(3),0)
b = c(-1/sqrt(3),0)
c = c(0,1)

x = sapply(list(a,b,c, split_s), FUN=function(x) x[1], simplify=TRUE)
y = sapply(list(a,b,c, split_s), FUN=function(x) x[2], simplify=TRUE)

#quartets$B = quartets$A & quartets$split1 > -1

set.seed(1545 + 6) # current time + millisecond on stopwatch
good_sample = sample(x=(1:nrow(quartets))[              quartets$B], size=3, replace=F)
bad_sample  = sample(x=(1:nrow(quartets))[quartets$A & !quartets$B], size=3, replace=F)

# 2022-08-16: trying to investigate a handful of outliers,
# which are off both the major arm AND its extension

quartets$arm_error = (quartets$split2 - quartets$split3)^2
summary(quartets$arm_error[quartets$B])
summary(quartets$arm_error[quartets$B & quartets$sim_num == 1585]) # contains outliers
test=binom.test(x=1060, n=1060+686, p=0.5, alternative="two.sided", conf.level=0.95) # from sim 1585 quartet 4,5,6,7
# the p-value of this test is way lower than i expect, even for so many quartets--something else must be going on
