results = read.csv("results.csv")

results$d = results$d_0
results$d_0 = results$d * results$lambda
results$mu_lambda = results$mu / results$lambda
#results$l_m_m = results$lambda - results$mu
#results$l_p_m = results$lambda + results$mu

results$A_q = results$A / results$q # possibly anomalous / total
results$C_B = results$C / results$B # anomalous / possibly anomalous and analyzable
results$C_q = results$C / results$q # anomalous / total
results$A1_q = results$A1 / results$q
results$A2_q = results$A2 / results$q
results$C1_A1 = results$C1 / results$A1
results$C2_A2 = results$C2 / results$A2

softLogit = function(x) {
  # logit(p) = log(p/(1-p))
  x = boot::logit(x)
  y = x[!is.na(x)]
  s = sd(y[y!=-Inf & y!=Inf])
  m = min(y[y!=-Inf])
  M = max(y[y!=+Inf])
  # replace -Inf with smallest non-inf value - 1 sd
  y[y==-Inf] = m - s 
  # similar
  y[y== Inf] = M + s
  x[!is.na(x)] = y
  return(x)
}

results$logitAq = softLogit(results$A_q) 
results$logitCB = softLogit(results$C_B)
results$logitCq = softLogit(results$C_q)
results$logitA1q = softLogit(results$A1_q)
results$logitA2q = softLogit(results$A2_q)
results$logitC1A1 = softLogit(results$C1_A1)
results$logitC2A2 = softLogit(results$C2_A2)

fit1 = lm(logitCq ~ ntaxa + lambda + mu_lambda + M + d_0,
          data=results)
summary(fit1)

fit00 = lm(logitA1q ~ ntaxa + lambda + mu_lambda + M + d_0,
          data=results)
fit01 = lm(logitA2q ~ ntaxa + lambda + mu_lambda + M + d_0,
           data=results)
fit10 = lm(logitC1A1 ~ ntaxa + lambda + mu_lambda + M + d_0,
           data=results)
fit11 = lm(logitC2A2 ~ ntaxa + lambda + mu_lambda + M + d_0,
           data=results)

jitter1 = rnorm(n=nrow(results), sd=0.2)
jitter2 = rnorm(n=nrow(results), sd=0.1)

summary(fit00); plot(y=results$logitA1q, x=results$d_0)
summary(fit01); plot(y=results$logitA2q, x=results$d_0)
summary(fit10); plot(y=results$logitC1A1 + jitter1, x=log(results$lambda) + jitter2)
summary(fit11); plot(y=results$logitC2A2 + jitter1, x=log(results$lambda) + jitter2)

