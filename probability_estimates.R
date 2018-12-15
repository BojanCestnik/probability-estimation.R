#####################################################################################
## Project: Probability estimation from small samples
## Script purpose: Experimental framework for analysing errors of probability 
##                 estimation methods on small samples
## Included methods: relative frequency, Laplace's rule, Piegat's formula, m-estimate
## Date: 15. 12. 2018
## Author: Bojan Cestnik
#####################################################################################

# parameters for testing: compare the estimated probability to the known prior
# 1. prior probability of 1's in 0/1 a sample: 01, 05, 10, 15, ..., 90, 95, 99
# 2. number of instances from which the probability is estimated
# 3. method for frobability estimation
# 4. the role of prior probability pa in m-estimate
# 5. the role of parameter m in m-estimate

### Setup of experimental data and functions - BEGIN

# distinct auhentic probabilities for sample generation 
pbase_global <- seq(0.0, 1.0, length = 21)
pbase_global[1] <- 0.01
pbase_global[21] <- 0.99
pbase_global

# generate samples of n (100500) bernoulli binomial instances (outcomes of 0 and 1) 
# for distinc authentic probabilities of obtaining 1
set.seed(123)
n <- 100500
inst_0_01 <- rbinom(n, 1, pbase_global[1])
inst_0_05 <- rbinom(n, 1, pbase_global[2])
inst_0_10 <- rbinom(n, 1, pbase_global[3])
inst_0_15 <- rbinom(n, 1, pbase_global[4])
inst_0_20 <- rbinom(n, 1, pbase_global[5])
inst_0_25 <- rbinom(n, 1, pbase_global[6])
inst_0_30 <- rbinom(n, 1, pbase_global[7])
inst_0_35 <- rbinom(n, 1, pbase_global[8])
inst_0_40 <- rbinom(n, 1, pbase_global[9])
inst_0_45 <- rbinom(n, 1, pbase_global[10])
inst_0_50 <- rbinom(n, 1, pbase_global[11])
inst_0_55 <- rbinom(n, 1, pbase_global[12])
inst_0_60 <- rbinom(n, 1, pbase_global[13])
inst_0_65 <- rbinom(n, 1, pbase_global[14])
inst_0_70 <- rbinom(n, 1, pbase_global[15])
inst_0_75 <- rbinom(n, 1, pbase_global[16])
inst_0_80 <- rbinom(n, 1, pbase_global[17])
inst_0_85 <- rbinom(n, 1, pbase_global[18])
inst_0_90 <- rbinom(n, 1, pbase_global[19])
inst_0_95 <- rbinom(n, 1, pbase_global[20])
inst_0_99 <- rbinom(n, 1, pbase_global[21])

# include the generated samples in a list
inst_list <- list(inst_0_01, inst_0_05, inst_0_10, inst_0_15, inst_0_20, inst_0_25, inst_0_30, inst_0_35, inst_0_40,
                  inst_0_45, inst_0_50, inst_0_55, inst_0_60, inst_0_65, inst_0_70, inst_0_75, inst_0_80, inst_0_85,
                  inst_0_90, inst_0_95, inst_0_99)

# generate a new sample as a union of two sub-samples
generate.inst <- function(inst_1, start_inst_1, len_inst_1, inst_2, start_inst_2, len_inst_2) 
{
  my_inst <- array()
  for (i in 1:len_inst_1) {
    my_inst[i] <- inst_1[start_inst_1-1+i]
  }
  for (i in 1:len_inst_2) {
    my_inst[len_inst_1+i] <- inst_2[start_inst_2-1+i]
  }
  my_inst
}

# display relative frequencies of the generated samples
t <- table(inst_0_01); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_05); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_10); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_15); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_20); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_25); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_30); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_35); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_40); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_45); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_50); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_55); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_60); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_65); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_70); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_75); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_80); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_85); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_90); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_95); t[[2]]/(t[[1]]+t[[2]]);
t <- table(inst_0_99); t[[2]]/(t[[1]]+t[[2]]);

## helper function for m-estimate of probability (parameters p_a and m) from a sample inst_p, 
# starting at start_inst, taking len_inst samples
mprob <- function(inst_p, start_inst, len_inst, pa, m) 
{
  inst <- inst_p[start_inst:(start_inst+len_inst-1)]
  t0 <- length(which(inst==0))
  t1 <- length(which(inst==1))
  px <- (t1 + (pa * m)) / (t0 + t1 + m)
  px
}

# methods for estimating probabilities: relfr, flat, laplace, piegat, cestnik, zadrozny

prob.relfr <- function(inst_p, start_inst, len_inst, pa=0.5, m=0) 
{
  mprob(inst_p, start_inst, len_inst, pa, m)
}

prob.flat <- function(inst_p, start_inst, len_inst, pa=0.5, m=0) 
{
  est <- mprob(inst_p, start_inst, len_inst, pa, m)
  dif <- 0.5 * (est - 0.5)
  est <- est - dif
  est
}

prob.laplace <- function(inst_p, start_inst, len_inst, pa=0.5, m=2) 
{
  mprob(inst_p, start_inst, len_inst, pa, m)
}

prob.piegat <- function(inst_p, start_inst, len_inst, pa=0.5, m=sqrt(2)) 
{
  mprob(inst_p, start_inst, len_inst, pa, m)
}

prob.cestnik <- function(inst_p, start_inst, len_inst, pa, m) 
{
  mprob(inst_p, start_inst, len_inst, pa, m)
}

prob.zadrozny <- function(inst_p, start_inst, len_inst, pa, m) 
{
  new_m <- 10.0 / pa
  mprob(inst_p, start_inst, len_inst, pa, new_m)
}

## helper fuction for randomly adding or substracting a small number 0.0001 to argument x
aprox <- function(x) {
  y <- sample(1:2, 1)
  if (y == 1) {
    dx <- -0.0001
  } else {
    dx <- 0.0001
  }
  res <- x + dx
  res
}

# classifier measures accuracy, sensitivity and specificity
accuracy <- function(cm) {
  # cm = array(TP, FP, FN, TN) 
  acc <- (cm[1]+cm[4])/(cm[1]+cm[2]+cm[3]+cm[4])
  acc
}

sensitivity <- function(cm) {
  # cm = array(TP, FP, FN, TN) 
  acc <- (cm[1])/(cm[1]+cm[3])
  acc
}

specificity <- function(cm) {
  # cm = array(TP, FP, FN, TN) 
  acc <- (cm[4])/(cm[2]+cm[4])
  acc
}

## calculate absolute error AE of the estimation_method in the experimental framework
# parameters: 
#   estimation_method: string name of the estimation method 
#   sample_testing: integer index of the sample for error testing (1 - 21)
#   sample_training: integer index of the sample for performing the estimation
#   sample_training_distort: integer sample ditortion ds for introducing noise in the sample
#   start_inst: integer index of the instance in the sample that is first in the training sub-sample
#   len_inst: integer length of the training sub-sample
#   pa: float prior probability for the estimation
#   m: float m for m-estimate
#   padistort: float distortion of prior probability
prob.AE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                     start_inst, len_inst, pa, m, padistort) {
  # determine actual_sample_training from sample_training and taking sample_training_distort into account
  actual_sample_training <- sample_training
  if (sample_training_distort > 0) { # random distortion
    line_min <- sample_training - sample_training_distort
    line_min <- max(1, line_min)
    line_max <- sample_training + sample_training_distort
    line_max <- min(21, line_max)
    actual_sample_training <- sample(c(line_min:line_max), 1)
  } else if (sample_training_distort < 0) { # extreme distortion
    line_min <- sample_training + sample_training_distort
    line_min <- max(1, line_min)
    line_max <- sample_training - sample_training_distort
    line_max <- min(21, line_max)
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      actual_sample_training <- line_min
    } else {
      actual_sample_training <- line_max
    }
  }
  # determine actual panew from pa and padistort
  if (padistort < -0.0001) { # use extreme distortion padistort
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      panew <- pa - padistort
    } else {
      panew <- pa + padistort
    }
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  } else if (padistort < 0.0001) { # no distortion of pa
    panew <- pa
  } else { # radmoly select panew from the inteval
    minpa <- pa - padistort
    maxpa <- pa + padistort
    if  (minpa < 0.01) {
      minpa <- 0.01
    }
    if (maxpa > 0.99) {
      maxpa <- 0.99
    }
    panew <- runif(1, min = minpa, max = maxpa)
  }

  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.flat") {
    px <- prob.flat(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst_list[[actual_sample_training]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.zadrozny") {
    px <- prob.zadrozny(inst_list[[actual_sample_training]], start_inst, len_inst, panew, m)
  }
  errb <- abs(px - pbase_global[sample_testing]) 
  errb
}

## length of subsample that contains 4s+4f in sample inst starting at fromi 
comp.dist.fromi <- function(inst, fromi) {
  s <- 0
  f <- 0
  found <- FALSE
  for (j in fromi:(fromi+1500)) {
    if (!found) {
      if (inst[j] == 0) {
        f <- f + 1
      } else {
        s <- s + 1 
      }
      if ((f >= 4) & (s >= 4)) {
        found <- TRUE
        dd <- j - fromi + 1
      } else if (j >= 100500) {
        found <- TRUE
        dd <- j - fromi + 1
      }
    }
  }
  dd
}

## calculate mean absolute error MAE of estimation_method in the experimental framework
prob.sum.MAE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                         len_inst, pa, m, padistort = 0.0)
{
  err <- 0 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      if (len_inst == 0) {
        len_inst_new <- comp.dist.fromi(inst_list[[sample_training]], i)
      } else {
        len_inst_new <- len_inst
      }
      err <- err + prob.AE(estimation_method, sample_testing, sample_training, sample_training_distort, 
                            i, len_inst_new, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  err <- err / cnt
  err
}

## calculate SD of estimation_method in the experimental framework
prob.sd.MAE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                         len_inst, pa, m, padistort = 0.0)
{
  err <- array() 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      if (len_inst == 0) {
        len_inst_new <- comp.dist.fromi(inst_list[[sample_training]], i)
      } else {
        len_inst_new <- len_inst
      }
      err[i] <- prob.AE(estimation_method, sample_testing, sample_training, sample_training_distort, 
                            i, len_inst_new, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  errx <- sd(err)
  errx
}

### computation for MSE and RMSE
## calculate mean squared error AE of the estimation_method in the experimental framework
# parameters: see prob.AE
prob.MSE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                     start_inst, len_inst, pa, m, padistort) {
  # determine actual_sample_training 
  actual_sample_training <- sample_training
  if (sample_training_distort > 0) {
    line_min <- sample_training - sample_training_distort
    line_min <- max(1, line_min)
    line_max <- sample_training + sample_training_distort
    line_max <- min(21, line_max)
    actual_sample_training <- sample(c(line_min:line_max), 1)
  } else if (sample_training_distort < 0) {
    line_min <- sample_training + sample_training_distort
    line_min <- max(1, line_min)
    line_max <- sample_training - sample_training_distort
    line_max <- min(21, line_max)
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      actual_sample_training <- line_min
    } else {
      actual_sample_training <- line_max
    }
  }
  # determine actual panew
  if (padistort < -0.0001) { # use extreme distortion padistort
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      panew <- pa - padistort
    } else {
      panew <- pa + padistort
    }
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  } else if (padistort < 0.0001) { # no distortion of pa
    panew <- pa
  } else { # radmoly select panew from the inteval
    minpa <- pa - padistort
    maxpa <- pa + padistort
    if  (minpa < 0.01) {
      minpa <- 0.01
    }
    if (maxpa > 0.99) {
      maxpa <- 0.99
    }
    panew <- runif(1, min = minpa, max = maxpa)
  }
  
  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.flat") {
    px <- prob.flat(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst_list[[actual_sample_training]], start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst_list[[actual_sample_training]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.zadrozny") {
    px <- prob.zadrozny(inst_list[[actual_sample_training]], start_inst, len_inst, panew, m)
  }
  errb <- (px - pbase_global[sample_testing])^2 
  errb
}

prob.sum.RMSE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                          len_inst, pa, m, padistort = 0.0)
{
  err <- 0 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      err <- err + prob.MSE(estimation_method, sample_testing, sample_training, sample_training_distort, 
                            i, len_inst, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  err <- sqrt(err / cnt)
  err
}

prob.sd.RMSE <- function(estimation_method, sample_testing, sample_training, sample_training_distort, 
                          len_inst, pa, m, padistort = 0.0)
{
  err <- array()
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      err[i] <- prob.MSE(estimation_method, sample_testing, sample_training, sample_training_distort, 
                            i, len_inst, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  errx <- sqrt(sd(err))
  errx
}

## ADIFF: average difference between two MAE arrays 
prob.diff.MAE <- function(MAE1, MAE2)
{
  myDiff <- array()
  for (i in 1:21) {
    myDiff[i] <- abs(MAE1[i]-MAE2[i])
  }
  myMinDiff <- min(myDiff)
  myMeanDiff <- mean(myDiff)
  mySdDiff <- sd(myDiff)
  myMaxDiff <- max(myDiff)
  res <- c(myMinDiff, myMeanDiff, mySdDiff, myMaxDiff)
  res
}

### Setup of experimental data and functions - END

### Figure: Comparison of the probability estimations with relative frequency and Lapalace's rule
n <- 100500
set.seed(123)
inst_70 <- rbinom(n, 1, 0.70)
inst_70[1:7]

n <- 7
rf <- array()
lap <- array()
target <- array()
for (l in 1:n) {
  rf[l] <- prob.relfr(inst_70, 1, l)
  lap[l] <- prob.laplace(inst_70, 1, l)
  target[l] <- 0.7
}

p <- seq(1, n)
par(pch=21, col="black", bg="white") 
plot(p, rf, type="n", ylim=c(0, 1.0), xlab="sample size", ylab="probability", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, target, pch=23, col="black", type="l", lty=3, lwd=1) 
lines(p, rf, pch=21, col="black", type="l", lty=1, lwd=2) 
lines(p, lap, pch=21, col="black", type="l", lty=2, lwd=2) 

legend("topright", legend=c("relfr", "laplace"),
       col=c("black", "black"), lty=c(1,2), lwd=2,
       pch=c(NA, NA), cex=1.2)


### Figure The estimation with m-estimate
n <- 100
nn <- 10
ss <- 3
p <- 0.8
rf <- array()
pa <- array()
mm <- array()
mmm <- array()
for (m in 0:n) {
  cnt <- m+1
  rf[cnt] <- ss/nn
  pa[cnt] <- p
  mm[cnt] <- ss/nn * nn/(nn+m) + p * m/(nn+m)
  mmm[cnt] <- NA
  if (m == 10) {
    mmm[cnt] <- mm[cnt]
  }
}

p <- seq(0, n)
par(pch=21, col="black", bg="white") 
plot(p, rf, type="n", ylim=c(0, 1.0), xlab="m", xaxt="n", ylab="probability", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 
axis(side = 1, at = c(0, 10, 20, 40, 60, 80, 100), cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0)

lines(p, rf, pch=23, col="black", type="l", lty=3, lwd=2) 
lines(p, pa, pch=21, col="black", type="l", lty=2, lwd=2) 
lines(p, mm, pch=21, col="black", type="l", lty=1, lwd=2)
lines(p, mmm, pch=15, col="black", type="p", lty=1, lwd=1)
abline(v = 10)

legend("bottomright", legend=c("relative frequency", "prior probability", "m-estimate"),
       col=c("black", "black", "black"), lty=c(3,2,1), lwd=2,
       pch=c(NA, NA, NA), cex=1.2)

### Figure Smoothing with Laplace's rule and simples flat
method_set <- c("prob.relfr", "prob.flat", "prob.laplace")
distort_sample <- 0
distort_pa <- 0.0
l_set <- c(1)
m <- 2
err.list <- list()
sd.list <- list()
cnt <- 0
for (method in method_set) {
  for (l in l_set) {
    cnt <- cnt + 1
    err.new <- array()
    sd.new <- array()
    for (line in c(1:21)) {
      err.new[line] <- prob.sum.MAE(method, line, line, distort_sample, l, pbase_global[line], m, distort_pa)
      sd.new[line] <- prob.sd.MAE(method, line, line, distort_sample, l, pbase_global[line], m, distort_pa)
      print(paste(method, " sample:", line, " distort_sample:", distort_sample, " sample size:", l, " m:", m, " pa:", pbase_global[line], " distort_pa:", distort_pa, " > "))
    }
    err.list[[cnt]] <- err.new
    sd.list[[cnt]] <- sd.new
  }
}

mean(err.list[[1]]); sd(err.list[[1]]); max(err.list[[1]]); 
mean(err.list[[2]]); sd(err.list[[2]]); max(err.list[[2]]); 
mean(err.list[[3]]); sd(err.list[[3]]); max(err.list[[3]]); 
err.list[[1]]
err.list[[2]]
err.list[[3]]

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=21, col="black", bg="white") 
plot(p, err.new, type="n", ylim=c(0, 0.50), ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[2]], pch=22, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[3]], pch=23, col="black", type="o", lty=2, lwd=1) 

legend("topright", legend=c("relfr", "flat", "laplace"),
       col=c("black", "black", "black"), lty=c(2,2,2), 
       pch=c(21, 22, 23), cex=1.2)


### Figures 1, 2, 3, 4, 5,  6,  7,  8,   10
myFigure <- "10"
method_set <- c("prob.relfr", "prob.laplace", "prob.piegat", "prob.cestnik")
l_set <- switch(myFigure, "1" = c(1), "2" = c(2), "3" = c(3), "4" = c(10), "5" = c(100), 
                "6" = c(1), "7" = c(1), "8" = c(1), 
                "10" = c(0), c(1))
distort_sample <- 0 
distort_pa <- switch(myFigure, "6" = 0.3, "7" = 0.5, "8" = 1.0, 0.0)
m <- 2
err.list <- list()
sd.list <- list()
cnt <- 0
for (method in method_set) {
  for (l in l_set) {
    cnt <- cnt + 1
    err.new <- array()
    sd.new <- array()
    for (line in c(1:21)) {
      err.new[line] <- prob.sum.MAE(method, line, line, distort_sample, l, pbase_global[line], m, distort_pa)
      sd.new[line] <- prob.sd.MAE(method, line, line, distort_sample, l, pbase_global[line], m, distort_pa)
      print(paste(method, " sample:", line, " distort_sample:", distort_sample, " sample size:", l, " m:", m, " pa:", pbase_global[line], " distort_pa:", distort_pa, " > "))
    }
    err.list[[cnt]] <- err.new
    sd.list[[cnt]] <- sd.new
  }
}

mean(err.list[[1]]); sd(err.list[[1]]); max(err.list[[1]]); 
mean(err.list[[2]]); sd(err.list[[2]]); max(err.list[[2]]); 
mean(err.list[[3]]); sd(err.list[[3]]); max(err.list[[3]]); 
mean(err.list[[4]]); sd(err.list[[4]]); max(err.list[[4]]); 

d1 <- c(mean(err.list[[1]]), 0, 0, 0)
d2 <- c(mean(err.list[[2]]), prob.diff.MAE(err.list[[2]], err.list[[1]])[2],  0, 0)
d3 <- c(mean(err.list[[3]]), prob.diff.MAE(err.list[[3]], err.list[[1]])[2],  prob.diff.MAE(err.list[[3]], err.list[[2]])[2], 0)
d4 <- c(mean(err.list[[4]]), prob.diff.MAE(err.list[[4]], err.list[[1]])[2],  prob.diff.MAE(err.list[[4]], err.list[[2]])[2], prob.diff.MAE(err.list[[4]], err.list[[3]])[2])
# AMAE and ADIFF
round(d1, 4)
round(d2, 4)
round(d3, 4)
round(d4, 4)

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=21, col="black", bg="white") 
mylim <- switch(myFigure, "7" = c(0, 0.15), "8" = c(0, 0.05), 
                "33" = c(0, 0.18), "34" = c(0, 0.17), "35" = c(0, 0.15), "36" = c(0, 0.15), "37" = c(0, 0.15), 
                c(0, 0.50))  
plot(p, err.new, type="n", ylim=mylim, ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[2]], pch=22, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[3]], pch=23, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[4]], pch=24, col="black", type="o", lty=2, lwd=1) 

if (myFigure == "1") {
  abline(v=c(0.14,0.29), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  abline(v=c(1-0.14,1-0.29), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  text(x=0.07, y=0.01, labels="extra-low", cex=1.1)
  text(x=0.22, y=0.01, labels="low", cex=1.1)
  text(x=0.50, y=0.01, labels="medium", cex=1.1)
  text(x=0.78, y=0.01, labels="high", cex=1.1)
  text(x=0.93, y=0.01, labels="extra-high", cex=1.1)
}

if (myFigure == "2") {
  abline(v=c(0.12,0.22), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  abline(v=c(1-0.12,1-0.22), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  text(x=0.06, y=0.01, labels="extra-low", cex=1.1)
  text(x=0.17, y=0.01, labels="low", cex=1.1)
  text(x=0.50, y=0.01, labels="medium", cex=1.1)
  text(x=0.83, y=0.01, labels="high", cex=1.1)
  text(x=0.94, y=0.01, labels="extra-high", cex=1.1)
}

if (myFigure == "3") {
  abline(v=c(0.11,0.18), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  abline(v=c(1-0.11,1-0.18), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  text(x=0.04, y=0.01, labels="extra-low", cex=1.1)
  text(x=0.14, y=0.01, labels="low", cex=1.1)
  text(x=0.50, y=0.01, labels="medium", cex=1.1)
  text(x=0.86, y=0.01, labels="high", cex=1.1)
  text(x=0.96, y=0.01, labels="extra-high", cex=1.1)
}

if (myFigure == "4") {
  abline(v=c(0.06,0.08), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  abline(v=c(1-0.06,1-0.08), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
  text(x=0.02, y=0.01, labels="extra-", cex=1.1)
  text(x=0.02, y=0.00, labels="low", cex=1.1)
  text(x=0.07, y=0.02, labels="l", cex=1.1)
  text(x=0.07, y=0.01, labels="o", cex=1.1)
  text(x=0.07, y=0.00, labels="w", cex=1.1)
  text(x=0.50, y=0.01, labels="medium", cex=1.1)
  text(x=0.93, y=0.03, labels="h", cex=1.1)
  text(x=0.93, y=0.02, labels="i", cex=1.1)
  text(x=0.93, y=0.01, labels="g", cex=1.1)
  text(x=0.93, y=0.00, labels="h", cex=1.1)
  text(x=0.98, y=0.01, labels="extra-", cex=1.1)
  text(x=0.98, y=0.00, labels="high", cex=1.1)
}

if (myFigure %in% c("1", "2", "3", "4", "5", "6", "7", "8")) {
  legend("topleft", legend=c("relfr", "laplace", "piegat", "cestnik"),
         col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
         pch=c(21, 22, 23, 24), cex=1.1)
}

if (myFigure %in% c("10")) {
  legend("topright", legend=c("relfr", "laplace", "piegat", "cestnik"),
       col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
       pch=c(21, 22, 23, 24), cex=1.1)
}


library(ggplot2)
library(forcats)
library(reshape2)
library(scmamp)

### Figure 9
# comparison of performance wit the approach by Demšar & Garcia et al.: 
# https://cran.r-project.org/web/packages/scmamp/vignettes/Statistical_assessment_of_the_differences.html
# Statistical Assessment of the Differences
CDy <- read.csv("W:/CDy99.csv", stringsAsFactors = FALSE, sep=";", dec=",")
CDy
CDy$ind <- NULL
CDx <- -CDy
names(CDx)
mean(CDx$relfr)
mean(CDx$laplace)
mean(CDx$piegat)
mean(CDx$cestnik)
mean(CDx$cestnik_pd0.3)
mean(CDx$cestnik_pd0.5)
mean(CDx$cestnik_pd1.0)

plotDensities(data=CDx, size=1.1)

qqplot <- qqplotGaussian (CDx[,"cestnik"], size=5 , col="orchid")
qqplot + theme_classic()

friedmanTest(CDx)

imanDavenportTest(CDx)

friedmanAlignedRanksTest(CDx)

quadeTest(CDx)

#Pairwise differences
test <- nemenyiTest (CDx, alpha=0.01)
test

test$diff.matrix

abs(test$diff.matrix) > test$statistic

plotCD (CDx, alpha=0.01, cex=0.9) # Figure 9

