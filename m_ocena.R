### Experimental framework for probability estimation methods 
### relfr, laplace, piegat, cestnik's m-estimate, ... 

# parameters for testing: compare the estimated probability to the known prior
# 1. prior probability of 1's in 0/1 a sample: 01, 05, 10, 15, ..., 90, 95, 99
# 2. number of instances from which the probability is estimated
# 3. method for frobability estimation
# 4. the role of prior probability pa in m-estimate
# 5. the role of parameter m in m-estimate

### Setup of experimental data and functions - BEGIN

pbase_global <- seq(0.0, 1.0, length = 21)
pbase_global[1] <- 0.01
pbase_global[21] <- 0.99
pbase_global

# generate n (100500) bernoulli binomial instances (outcomes of 0 and 1) 
# with a given probablility of obtaining 1
set.seed(123)
n <- 100500
inst_0_01 <- rbinom(n, 1, 0.01)
inst_0_05 <- rbinom(n, 1, 0.05)
inst_0_10 <- rbinom(n, 1, 0.10)
inst_0_15 <- rbinom(n, 1, 0.15)
inst_0_20 <- rbinom(n, 1, 0.20)
inst_0_25 <- rbinom(n, 1, 0.25)
inst_0_30 <- rbinom(n, 1, 0.30)
inst_0_35 <- rbinom(n, 1, 0.35)
inst_0_40 <- rbinom(n, 1, 0.40)
inst_0_45 <- rbinom(n, 1, 0.45)
inst_0_50 <- rbinom(n, 1, 0.50)
inst_0_55 <- rbinom(n, 1, 0.55)
inst_0_60 <- rbinom(n, 1, 0.60)
inst_0_65 <- rbinom(n, 1, 0.65)
inst_0_70 <- rbinom(n, 1, 0.70)
inst_0_75 <- rbinom(n, 1, 0.75)
inst_0_80 <- rbinom(n, 1, 0.80)
inst_0_85 <- rbinom(n, 1, 0.85)
inst_0_90 <- rbinom(n, 1, 0.90)
inst_0_95 <- rbinom(n, 1, 0.95)
inst_0_99 <- rbinom(n, 1, 0.99)

inst_list <- list(inst_0_01, inst_0_05, inst_0_10, inst_0_15, inst_0_20, inst_0_25, inst_0_30, inst_0_35, inst_0_40,
                  inst_0_45, inst_0_50, inst_0_55, inst_0_60, inst_0_65, inst_0_70, inst_0_75, inst_0_80, inst_0_85,
                  inst_0_90, inst_0_95, inst_0_99)


