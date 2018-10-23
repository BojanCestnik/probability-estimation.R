### Experimental framework for probability estimation methods

### Setup of experimental data and functions - BEGIN

pbase_global <- seq(0.0, 1.0, length = 11)
pbase_global[1] <- 0.01
pbase_global[11] <- 0.99
pbase_global

# generate n (1000) bernoulli binomial instances (outcomes of 0 and 1)
# with a given probablility of obtaining 1
set.seed(123)
n <- 1000
inst_0_01 <- rbinom(n, 1, 0.01)
inst_0_10 <- rbinom(n, 1, 0.10)
inst_0_20 <- rbinom(n, 1, 0.20)
inst_0_30 <- rbinom(n, 1, 0.30)
inst_0_40 <- rbinom(n, 1, 0.40)
inst_0_50 <- rbinom(n, 1, 0.50)
inst_0_60 <- rbinom(n, 1, 0.60)
inst_0_70 <- rbinom(n, 1, 0.70)
inst_0_80 <- rbinom(n, 1, 0.80)
inst_0_90 <- rbinom(n, 1, 0.90)
inst_0_99 <- rbinom(n, 1, 0.99)
