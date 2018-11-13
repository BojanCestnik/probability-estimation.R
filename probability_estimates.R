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

# helper function for m-estimate of probability (parameters p_a and m) from a sample inst_p, 
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
# cestnik_optimal, proxim, proxim_1..9

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

prob.cestnik_optimal <- function(inst_p, start_inst, len_inst, pa, pa_proxim) 
{
  # pa_proxim = 0 (m), 1 0.3, 2 0.5, 3 1.0
  # set m depending on pa, len_inst and pa_proxim
  if (pa_proxim == 0) {
    m_new <- 2.0
  } else {
    pa_border <- pa
    if (pa_border > 0.5) {
      pa_border <- 1.0 - pa_border
    }
    if (pa_border < 0.175) {
      if (pa_border < 0.025) {
        if (pa_proxim == 1) {
          if (len_inst <= 5) {
            m_new <- 0.1
          } else {
            m_new <- 0.5
          }
        } else {
          m_new <- 0.1
        } 
      } else if (pa_border < 0.075) {
        if (pa_proxim == 1) {
          if (len_inst <= 2) {
            m_new <- 0.5
          } else if (len_inst <= 5) {
            m_new <- 1
          } else {
            m_new <- 2
          }
        } else if (pa_proxim == 2) {
          if (len_inst <= 1) {
            m_new <- 0.1
          } else if (len_inst <= 5) {
            m_new <- 0.5
          } else {
            m_new <- 1
          }
        } else if (pa_proxim == 3) {
          if (len_inst <= 3) {
            m_new <- 0.1
          } else {
            m_new <- 0.5
          }
        }
      } else if (pa_border < 0.135) {
        if (pa_proxim == 1) {
          if (len_inst <= 2) {
            m_new <- 1.5
          } else if (len_inst <= 5) {
            m_new <- 2.5
          } else {
            m_new <- 3.5
          }
        } else if (pa_proxim == 2) {
          if (len_inst <= 2) {
            m_new <- 0.5
          } else if (len_inst <= 5) {
            m_new <- 1.5
          } else {
            m_new <- 2
          }
        } else if (pa_proxim == 3) {
          if (len_inst <= 1) {
            m_new <- 0.1
          } else if (len_inst <= 5) {
            m_new <- 0.5
          } else {
            m_new <- 1.0
          }
        }
      } else {
        if (pa_proxim == 1) {
          if (len_inst <= 1) {
            m_new <- 3
          } else if (len_inst <= 5) {
            m_new <- 4.5
          } else {
            m_new <- 5
          }
        } else if (pa_proxim == 2) {
          if (len_inst <= 2) {
            m_new <- 1
          } else if (len_inst <= 6) {
            m_new <- 2
          } else {
            m_new <- 2.5
          }
        } else if (pa_proxim == 3) {
          if (len_inst <= 2) {
            m_new <- 0.5
          } else {
            m_new <- 1.0
          }
        }
      }
      m_new <- 2
    } else {
      if (pa_proxim == 1) {
        m_new <- 6
      } else if (pa_proxim == 2) {
        m_new <- 3
      } else if (pa_proxim == 3) {
        m_new <- 2
      }
    }
  }
  mprob(inst_p, start_inst, len_inst, pa, m_new)
}

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

prob.proxim <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < 0.41) {
    new_pa <- 0.3
  } else  if (pa > 0.59) {
    new_pa <- 0.7
  } else {
    new_pa <- 0.5
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_1 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  new_pa <- 0.5
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_2 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < aprox(0.5)) {
    new_pa <- 0.333333
  } else {
    new_pa <- 0.666667
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_3 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < aprox(0.375)) {
    new_pa <- 0.25
  } else  if (pa < aprox(0.625)) {
    new_pa <- 0.5
  } else {
    new_pa <- 0.75
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_4 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < aprox(0.3)) {
    new_pa <- 0.2
  } else  if (pa < aprox(0.5)) {
    new_pa <- 0.4
  } else  if (pa < aprox(0.7)) {
    new_pa <- 0.6
  } else {
    new_pa <- 0.8
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_5 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < aprox(0.25)) {
    new_pa <- 0.166667
  } else  if (pa < aprox(0.416667)) {
    new_pa <- 0.333333
  } else  if (pa < aprox(0.583333)) {
    new_pa <- 0.5
  } else  if (pa < aprox(0.75)) {
    new_pa <- 0.666667
  } else {
    new_pa <- 0.833333
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

prob.proxim_9 <- function(inst_p, start_inst, len_inst, pa, m) 
{
  if (pa < aprox(0.15)) {
    new_pa <- 0.1
  } else  if (pa < aprox(0.25)) {
    new_pa <- 0.2
  } else  if (pa < aprox(0.35)) {
    new_pa <- 0.3
  } else  if (pa < aprox(0.45)) {
    new_pa <- 0.4
  } else  if (pa < aprox(0.55)) {
    new_pa <- 0.5
  } else  if (pa < aprox(0.65)) {
    new_pa <- 0.6
  } else  if (pa < aprox(0.75)) {
    new_pa <- 0.7
  } else  if (pa < aprox(0.85)) {
    new_pa <- 0.8
  } else {
    new_pa <- 0.9
  }
  mprob(inst_p, start_inst, len_inst, new_pa, m)
}

# classifier measures accuracy, sensitivity and specificity
accuracy <- function(cm) {
  #cm = array(TP, FP, FN, TN) 
  acc <- (cm[1]+cm[4])/(cm[1]+cm[2]+cm[3]+cm[4])
  acc
}

sensitivity <- function(cm) {
  #cm = array(TP, FP, FN, TN) 
  acc <- (cm[1])/(cm[1]+cm[3])
  acc
}

specificity <- function(cm) {
  #cm = array(TP, FP, FN, TN) 
  acc <- (cm[4])/(cm[2]+cm[4])
  acc
}

### computation of absolute errors AE and mean absolute errors MAE of various methods

prob.AE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                     start_inst, len_inst, pa, m, padistort) {
  #determine line_learning (taking line_learning_distort into account)
  actual_line_learning <- line_learning
  if (line_learning_distort > 0) {
    line_min <- line_learning - line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning + line_learning_distort
    line_max <- min(21, line_max)
    actual_line_learning <- sample(c(line_min:line_max), 1)
  } else if (line_learning_distort < 0) {
    line_min <- line_learning + line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning - line_learning_distort
    line_max <- min(21, line_max)
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      actual_line_learning <- line_min
    } else {
      actual_line_learning <- line_max
    }
  }
  #determine actual pa (taking padistort into account)
  if (padistort < -0.0001) {
    #panew <- 0.5 + runif(1, min = padistort, max = -padistort)
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
  } else if (padistort < 0.0001) {
    panew <- pa
  } else {
    minpa <- pa - padistort
    maxpa <- pa + padistort
    if  (minpa < 0.01) {
      minpa <- 0.01
    }
    if (maxpa > 0.99) {
      maxpa <- 0.99
    }
    panew <- runif(1, min = minpa, max = maxpa)
    #if (pa >= 0.5) {panew <- pa - padistort} else {panew <- pa + padistort}
    #panew <- pa + padistort
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  }

  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.flat") {
    px <- prob.flat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.cestnik_optimal") {
    px <- prob.cestnik_optimal(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.zadrozny") {
    px <- prob.zadrozny(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim") {
    px <- prob.proxim(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_1") {
    px <- prob.proxim_1(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_2") {
    px <- prob.proxim_2(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_3") {
    px <- prob.proxim_3(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_4") {
    px <- prob.proxim_4(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_5") {
    px <- prob.proxim_5(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_9") {
    px <- prob.proxim_9(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  }
  errb <- abs(px - pbase_global[line_testing]) #added as a common sense error estimation
  errb
}

prob.sum.MAE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                         len_inst, pa, m, padistort = 0.0)
{
  err <- 0 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      if (len_inst == 0) {
        len_inst_new <- comp.dist.fromi(inst_list[[line_learning]], i)
      } else {
        len_inst_new <- len_inst
      }
      err <- err + prob.AE(estimation_method, line_testing, line_learning, line_learning_distort, 
                            i, len_inst_new, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  err <- err / cnt
  err
}

prob.sd.MAE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                         len_inst, pa, m, padistort = 0.0)
{
  err <- array() 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      if (len_inst == 0) {
        len_inst_new <- comp.dist.fromi(inst_list[[line_learning]], i)
      } else {
        len_inst_new <- len_inst
      }
      err[i] <- prob.AE(estimation_method, line_testing, line_learning, line_learning_distort, 
                            i, len_inst_new, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  errx <- sd(err)
  errx
}

probi.AE <- function(estimation_method, inst, start_inst, len_inst, pa, m) {
  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst, start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst, start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst, start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst, start_inst, len_inst, pa, m)
  } 
  errb <- abs(px - pa) #added as a common sense error estimation
  errb
}

probi.sum.MAE <- function(estimation_method, inst, 
                         len_inst, pa, m)
{
  err <- 0 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      if (len_inst == 0) {
        len_inst_new <- comp.dist.fromi(inst_list[[line_learning]], i)
      } else {
        len_inst_new <- len_inst
      }
      err <- err + probi.AE(estimation_method, inst, 
                           i, len_inst_new, pa, m)
      cnt <- cnt + 1
    }
  }
  err <- err / cnt
  err
}

### computation for alfa

prob.alfa <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                    start_inst, len_inst, pa, m, padistort, allowed_error) {
  #determine line_learning (taking line_learning_distort into account)
  actual_line_learning <- line_learning
  if (line_learning_distort > 0) {
    line_min <- line_learning - line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning + line_learning_distort
    line_max <- min(21, line_max)
    actual_line_learning <- sample(c(line_min:line_max), 1)
  } else if (line_learning_distort < 0) {
    line_min <- line_learning + line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning - line_learning_distort
    line_max <- min(21, line_max)
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      actual_line_learning <- line_min
    } else {
      actual_line_learning <- line_max
    }
  }
  #determine actual pa (taking padistort into account)
  if (padistort < -0.0001) {
    #panew <- 0.5 + runif(1, min = padistort, max = -padistort)
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
  } else if (padistort < 0.0001) {
    panew <- pa
  } else {
    minpa <- pa - padistort
    maxpa <- pa + padistort
    if  (minpa < 0.01) {
      minpa <- 0.01
    }
    if (maxpa > 0.99) {
      maxpa <- 0.99
    }
    panew <- runif(1, min = minpa, max = maxpa)
    #if (pa >= 0.5) {panew <- pa - padistort} else {panew <- pa + padistort}
    #panew <- pa + padistort
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  }
  
  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.flat") {
    px <- prob.flat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.cestnik_optimal") {
    px <- prob.cestnik_optimal(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.zadrozny") {
    px <- prob.zadrozny(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim") {
    px <- prob.proxim(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_1") {
    px <- prob.proxim_1(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_2") {
    px <- prob.proxim_2(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_3") {
    px <- prob.proxim_3(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_4") {
    px <- prob.proxim_4(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_5") {
    px <- prob.proxim_5(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_9") {
    px <- prob.proxim_9(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  }
  #erra <- c(0, 0, 0, 0) # TP FP FN TN
  #for (i in 1:100) {
  #  if (test_list[[line_testing]][i]==1) {
  #    if ((px > 0.4999) & (px < 0.5001)) {
  #      erra[1] <- erra[1] + 0.5
  #      erra[3] <- erra[3] + 0.5
  #    } else if (px > 0.5) {
  #      erra[1] <- erra[1] + 1
  #    } else {
  #      erra[3] <- erra[3] + 1
  #    }
  #  } else {
  #    if ((px > 0.4999) & (px < 0.5001)) {
  #      erra[2] <- erra[2] + 0.5
  #      erra[4] <- erra[4] + 0.5
  #    } else if (px > 0.5) {
  #      erra[2] <- erra[2] + 1
  #    } else {
  #      erra[4] <- erra[4] + 1
  #    }
  #  }
  #}
  
  old_err <- abs(px - pbase_global[line_testing]) #added as a common sense error estimation
  
  #pbase_global[line_testing] - authentic p of the sample
  #px - estimated p
  #panew - prior probability
  rf <- prob.relfr(inst_list[[actual_line_learning]], start_inst, len_inst) # relative frequency
  #len_inst - number of instances (n)
  
  if (abs(panew - rf) < 0.00001) {
    #alfa1 <- runif(1, min = 0.0, max = 1.0)
    alfa1 <- 0.5
    #alfa1 <- 1.0
    alfa2 <- 2.0 # - allowed_error
    alfa3 <- 2.0 # + allowed_error
  } else {
    alfa1 <- (pbase_global[line_testing]-rf)/(panew-rf)
    alfa2 <- (pbase_global[line_testing]-allowed_error-rf)/(panew-rf)
    alfa3 <- (pbase_global[line_testing]+allowed_error-rf)/(panew-rf)
    #alfa <- 0.0
  }
  alfa <- min(c(alfa1, alfa2, alfa3))
  if (alfa > 1.0) {
    alfa <- 1.0
    #alfa <- 0.0
  } else if (alfa < 0.0) {
    alfa <- 0.0
  }
  #alfa
  
  if (alfa > 0.99) {
    alfa <- 0.99
  }
    
  mmm <- (alfa * len_inst)/(1.0 - alfa)
  if (mmm > 10.0) {
    mmm <- 10.0
  }
  px <- prob.cestnik(inst_list[[actual_line_learning]], start_inst, len_inst, panew, mmm)
  new_err <- abs(px - pbase_global[line_testing]) #added as a common sense error estimation
  c(mmm, old_err, new_err)
  #alfa
  
  #intrinsic error - can not be avoided
  #if (rf < panew) {
  #  if (pbase_global[line_testing] < rf) {
  #    mmm <- rf - pbase_global[line_testing]
  #  } else if (pbase_global[line_testing] > panew) {
  #    mmm <- pbase_global[line_testing] - panew
  #  } else {
  #    mmm <- 0.0
  #  }
  #} else {
  #  if (pbase_global[line_testing] < panew) {
  #    mmm <- panew - pbase_global[line_testing]
  #  } else if (pbase_global[line_testing] > rf) {
  #    mmm <- pbase_global[line_testing] - rf
  #  } else {
  #    mmm <- 0.0
  #  }
  #}
  #mmm  
}

prob.sum.alfa <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                         len_inst, pa, m, padistort = 0.0, allowed_error = 0.0)
{
  mmm <- 0
  old_err <- 0
  new_err <- 0
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      tmp <- prob.alfa(estimation_method, line_testing, line_learning, line_learning_distort, 
                             i, len_inst, pa, m, padistort, allowed_error)
      mmm <- mmm + tmp[1]
      old_err <- old_err + tmp[2]
      new_err <- new_err + tmp[3]
      cnt <- cnt + 1
    }
  }
  mmm <- mmm / cnt
  old_err <- old_err / cnt
  new_err <- new_err / cnt
  c(mmm, old_err, new_err)
}

prob.sd.alfa <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                        len_inst, pa, m, padistort = 0.0, allowed_error)
{
  mmm <- array()
  old_err <- array()
  new_err <- array()
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      tmp <- prob.alfa(estimation_method, line_testing, line_learning, line_learning_distort, 
                       i, len_inst, pa, m, padistort, allowed_error)
      mmm[i] <- tmp[1]
      old_err[i] <- tmp[2]
      new_err[i] <- tmp[3]
      cnt <- cnt + 1
    }
  }
  c(sd(mmm), sd(old_err), sd(new_err))
}

### computation for MSE and RMSE

prob.MSE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                     start_inst, len_inst, pa, m, padistort) {
  #determine line_learning (taking line_lerning_distort into account)
  actual_line_learning <- line_learning
  if (line_learning_distort > 0) {
    line_min <- line_learning - line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning + line_learning_distort
    line_max <- min(21, line_max)
    actual_line_learning <- sample(c(line_min:line_max), 1)
  } else if (line_learning_distort < 0) {
    line_min <- line_learning + line_learning_distort
    line_min <- max(1, line_min)
    line_max <- line_learning - line_learning_distort
    line_max <- min(21, line_max)
    tmp <- sample(c(1:2), 1)
    if (tmp == 1) {
      actual_line_learning <- line_min
    } else {
      actual_line_learning <- line_max
    }
  }
  #determine actual pa (taking padistort into account)
  if (padistort < -0.0001) {
    panew <- 0.5 + runif(1, min = padistort, max = -padistort)
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  } else if (padistort < 0.0001) {
    panew <- pa
  } else {
    minpa <- pa - padistort
    maxpa <- pa + padistort
    if  (minpa < 0.01) {
      minpa <- 0.01
    }
    if (maxpa > 0.99) {
      maxpa <- 0.99
    }
    panew <- runif(1, min = minpa, max = maxpa)
    #if (pa >= 0.5) {panew <- pa - padistort} else {panew <- pa + padistort}
    #panew <- pa + padistort
    if (panew < 0.01) {
      panew <- 0.01
    } else if (panew > 0.99) {
      panew <- 0.99
    }
  }
  
  if (estimation_method == "prob.relfr") {
    px <- prob.relfr(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.flat") {
    px <- prob.flat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.laplace") {
    px <- prob.laplace(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.piegat") {
    px <- prob.piegat(inst_list[[actual_line_learning]], start_inst, len_inst)
  } else if (estimation_method == "prob.cestnik") {
    px <- prob.cestnik(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.cestnik_optimal") {
    px <- prob.cestnik_optimal(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.zadrozny") {
    px <- prob.zadrozny(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim") {
    px <- prob.proxim(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_1") {
    px <- prob.proxim_1(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_2") {
    px <- prob.proxim_2(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_3") {
    px <- prob.proxim_3(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_4") {
    px <- prob.proxim_4(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_5") {
    px <- prob.proxim_5(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  } else if (estimation_method == "prob.proxim_9") {
    px <- prob.proxim_9(inst_list[[actual_line_learning]], start_inst, len_inst, panew, m)
  }
  #erra <- c(0, 0, 0, 0) # TP FP FN TN
  #for (i in 1:100) {
  #  if (test_list[[line_testing]][i]==1) {
  #    if ((px > 0.4999) & (px < 0.5001)) {
  #      erra[1] <- erra[1] + 0.5
  #      erra[3] <- erra[3] + 0.5
  #    } else if (px > 0.5) {
  #      erra[1] <- erra[1] + 1
  #    } else {
  #      erra[3] <- erra[3] + 1
  #    }
  #  } else {
  #    if ((px > 0.4999) & (px < 0.5001)) {
  #      erra[2] <- erra[2] + 0.5
  #      erra[4] <- erra[4] + 0.5
  #    } else if (px > 0.5) {
  #      erra[2] <- erra[2] + 1
  #    } else {
  #      erra[4] <- erra[4] + 1
  #    }
  #  }
  #}
  
  errb <- (px - pbase_global[line_testing])^2 #added as a common sense error estimation
  errb
}

prob.sum.RMSE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                          len_inst, pa, m, padistort = 0.0)
{
  err <- 0 
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      err <- err + prob.MSE(estimation_method, line_testing, line_learning, line_learning_distort, 
                            i, len_inst, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  err <- sqrt(err / cnt)
  err
}

prob.sd.RMSE <- function(estimation_method, line_testing, line_learning, line_learning_distort, 
                          len_inst, pa, m, padistort = 0.0)
{
  err <- array()
  cnt <- 0
  for (i in 1:100000) {
    if (TRUE) {
      err[i] <- prob.MSE(estimation_method, line_testing, line_learning, line_learning_distort, 
                            i, len_inst, pa, m, padistort)
      cnt <- cnt + 1
    }
  }
  errx <- sqrt(sd(err))
  errx
}

prob.min.MAE <- function(a1, a2, a3, a4)
{
  myMin <- array()
  for (i in 1:21) {
    myMin[i] <- min(c(a1[i], a2[i], a3[i], a4[i]))
  }
  myMin
}

prob.max.MAE <- function(a1, a2, a3, a4)
{
  myMax <- array()
  for (i in 1:21) {
    myMax[i] <- max(c(a1[i], a2[i], a3[i], a4[i]))
  }
  myMax
}

prob.diff.MAE <- function(a1, a2)
{
  myDiff <- array()
  for (i in 1:21) {
    myDiff[i] <- abs(a1[i]-a2[i])
  }
  myMinDiff <- min(myDiff)
  myMeanDiff <- mean(myDiff)
  mySdDiff <- sd(myDiff)
  myMaxDiff <- max(myDiff)
  res <- c(myMinDiff, myMeanDiff, mySdDiff, myMaxDiff)
  res
}

# length of subsample that contains 4s+4f in sample inst starting at fromi 
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
  # dd <- dd - 1
  dd
}

### Setup of experimental data and functions - END

### figure 1
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


### figure 2
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

### figure 3
method_set <- c("prob.relfr", "prob.flat", "prob.laplace")
distort_line <- 0
distort <- 0.0
l_set <- c(1)
m <- 2
err.list <- list()
sd.list <- list()
cnt <- 0
for (method in method_set) {
  if (method == "prob.cestnik_optimal") {
    m <- 3
  } else {
    m <- 2
  }
  for (l in l_set) {
    cnt <- cnt + 1
    err.new <- array()
    sd.new <- array()
    for (line in c(1:21)) {
      err.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      sd.new[line] <- prob.sd.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      print(paste(method, " line:", line, " size:", l, " m:", m, " pa:", pbase_global[line], " distort:", distort, " > "))
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


### figure 4, 5, 6, 7, 8,   25, 26, 27,   33, 34, 35, 36, 37
method_set <- c("prob.relfr", "prob.laplace", "prob.piegat", "prob.cestnik")
distort_line <- 0
distort <- 0.0
l_set <- c(1) #figure 4
l_set <- c(2) #figure 5
l_set <- c(3) #figure 6
l_set <- c(4)
l_set <- c(5)
l_set <- c(6)
l_set <- c(7)
l_set <- c(10) #figure 7
l_set <- c(100) #figure 8
l_set <- c(0) #figure 33, 34, 35, 36, 37
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
      err.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      sd.new[line] <- prob.sd.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      print(paste(method, " line:", line, " size:", l, " m:", m, " pa:", pbase_global[line], " distort:", distort, " > "))
    }
    err.list[[cnt]] <- err.new
    sd.list[[cnt]] <- sd.new
  }
}

mean(err.list[[1]]); sd(err.list[[1]]); max(err.list[[1]]); 
mean(err.list[[2]]); sd(err.list[[2]]); max(err.list[[2]]); 
mean(err.list[[3]]); sd(err.list[[3]]); max(err.list[[3]]); 
mean(err.list[[4]]); sd(err.list[[4]]); max(err.list[[4]]); 
err.list[[4]]
err.list[[1]][11]
err.list[[2]][11]
err.list[[3]][11]
err.list[[4]][11]
err.list[[1]][11]-err.list[[2]][11]
sd.list[[1]][11]
sd.list[[2]][11]
sd.list[[3]][11]
sd.list[[4]][11]

myMin <- prob.min.MAE(err.list[[1]], err.list[[2]], err.list[[3]], err.list[[4]])
myMax <- prob.max.MAE(err.list[[1]], err.list[[2]], err.list[[3]], err.list[[4]])
round(prob.diff.MAE(myMin, myMax), 4)
d1 <- c(mean(err.list[[1]]), 0, 0, 0)
d2 <- c(mean(err.list[[2]]), prob.diff.MAE(err.list[[2]], err.list[[1]])[2],  0, 0)
d3 <- c(mean(err.list[[3]]), prob.diff.MAE(err.list[[3]], err.list[[1]])[2],  prob.diff.MAE(err.list[[3]], err.list[[2]])[2], 0)
d4 <- c(mean(err.list[[4]]), prob.diff.MAE(err.list[[4]], err.list[[1]])[2],  prob.diff.MAE(err.list[[4]], err.list[[2]])[2], prob.diff.MAE(err.list[[4]], err.list[[3]])[2])
round(d1, 4)
round(d2, 4)
round(d3, 4)
round(d4, 4)

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=21, col="black", bg="white") 
plot(p, err.new, type="n", ylim=c(0, 0.50), ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[2]], pch=22, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[3]], pch=23, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[4]], pch=24, col="black", type="o", lty=2, lwd=1) 

#1
abline(v=c(0.14,0.29), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(1-0.14,1-0.29), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
text(x=0.07, y=0.01, labels="extra-low", cex=1.1)
text(x=0.22, y=0.01, labels="low", cex=1.1)
text(x=0.50, y=0.01, labels="medium", cex=1.1)
text(x=0.78, y=0.01, labels="high", cex=1.1)
text(x=0.93, y=0.01, labels="extra-high", cex=1.1)

#2
abline(v=c(0.12,0.22), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(1-0.12,1-0.22), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
text(x=0.06, y=0.01, labels="extra-low", cex=1.1)
text(x=0.17, y=0.01, labels="low", cex=1.1)
text(x=0.50, y=0.01, labels="medium", cex=1.1)
text(x=0.83, y=0.01, labels="high", cex=1.1)
text(x=0.94, y=0.01, labels="extra-high", cex=1.1)

#3
abline(v=c(0.11,0.18), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(1-0.11,1-0.18), col=c("gray", "gray"), lty=c(2,2), lwd=c(3, 3))
text(x=0.04, y=0.01, labels="extra-low", cex=1.1)
text(x=0.14, y=0.01, labels="low", cex=1.1)
text(x=0.50, y=0.01, labels="medium", cex=1.1)
text(x=0.86, y=0.01, labels="high", cex=1.1)
text(x=0.96, y=0.01, labels="extra-high", cex=1.1)

#10
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

legend("topleft", legend=c("relfr", "laplace", "piegat", "cestnik"),
       col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
       pch=c(21, 22, 23, 24), cex=1.1)

legend("topright", legend=c("relfr", "laplace", "piegat", "cestnik"),
       col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
       pch=c(21, 22, 23, 24), cex=1.1)

legend("bottomright", legend=c("relfr", "laplace", "piegat", "cestnik"),
       col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
       pch=c(21, 22, 23, 24), cex=1.1)

legend("bottomleft", legend=c("relfr", "laplace", "piegat", "cestnik"),
       col=c("black", "black", "black", "black"), lty=c(2,2,2,2), 
       pch=c(21, 22, 23, 24), cex=1.1)


### figures 9, 10, 11, 12, 13, 14, 15,   16, 17 (zadrozny),   18, 19, 20,   21, 22, 23, 24
method_set <- c("prob.relfr")
method_set <- c("prob.laplace")
method_set <- c("prob.piegat")
method_set <- c("prob.cestnik")
method_set <- c("prob.zadrozny")
method_set <- c("prob.proxim_2")
method_set <- c("prob.proxim_3")
distort_line <- 10
distort <- 1.0
l_set <- c(1,2,3,4,5,6,7)
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
      err.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      sd.new[line] <- prob.sd.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      print(paste(method, " line:", line, " size:", l, " m:", m, " pa:", pbase_global[line], " distort:", distort, " > "))
    }
    err.list[[cnt]] <- err.new
    sd.list[[cnt]] <- sd.new
  }
}

mean(err.list[[1]]); sd(err.list[[1]]); max(err.list[[1]]); mean(err.list[[1]]) + sd(err.list[[1]])
mean(err.list[[2]]); sd(err.list[[2]]); max(err.list[[2]]); mean(err.list[[2]]) + sd(err.list[[2]])
mean(err.list[[3]]); sd(err.list[[3]]); max(err.list[[3]]); mean(err.list[[3]]) + sd(err.list[[3]])
mean(err.list[[4]]); sd(err.list[[4]]); max(err.list[[4]]); mean(err.list[[4]]) + sd(err.list[[4]])
mean(err.list[[5]]); sd(err.list[[5]]); max(err.list[[5]]); mean(err.list[[5]]) + sd(err.list[[5]])
mean(err.list[[6]]); sd(err.list[[6]]); max(err.list[[6]]); mean(err.list[[6]]) + sd(err.list[[6]])
mean(err.list[[7]]); sd(err.list[[7]]); max(err.list[[7]]); mean(err.list[[7]]) + sd(err.list[[7]])

err.list[[2]][11]
sd.list[[2]][11]
err.list[[3]][11]
sd.list[[3]][11]

err.list[[4]][11]
sd.list[[4]][11]
err.list[[5]][11]
sd.list[[5]][11]

err.list[[6]][11]
sd.list[[6]][11]
err.list[[7]][11]
sd.list[[7]][11]

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=21, col="black", bg="white") 
plot(p, err.new, type="n", ylim=c(0, 0.5), ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[2]], pch=22, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[3]], pch=23, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[4]], pch=24, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[5]], pch=25, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[6]], pch=15, col="black", type="o", lty=2, lwd=1) 
lines(p, err.list[[7]], pch=16, col="black", type="o", lty=2, lwd=1) 

legend("topleft", legend=c("1", "2", "3", "4", "5", "6", "7"),
       col=c("black", "black", "black", "black", "black", "black", "black"), lty=c(2,2,2,2,2,2,2), 
       pch=c(21, 22, 23, 24, 25, 15, 16), cex=1.2)

legend("bottomleft", legend=c("1", "2", "3", "4", "5", "6", "7"),
       col=c("black", "black", "black", "black", "black", "black", "black"), lty=c(2,2,2,2,2,2,2), 
       pch=c(21, 22, 23, 24, 25, 15, 16), cex=1.2)


### figure 29 - single size sample
relfr.MAE <- c(0.3184851, 0.2191685, 0.1775729, 0.1525493, 0.1359476, 0.1245166, 0.1151591)
laplace.MAE <- c(0.2037209, 0.1708885, 0.1499013, 0.1358202, 0.1245937, 0.1154182, 0.1083315)
piegat.MAE <- c(0.1993492, 0.1670021, 0.1468770, 0.1325229, 0.1216835, 0.1127385, 0.1058672)
cestnik.MAE <- c(0.10616171, 0.10958427, 0.10654376, 0.10169955, 0.09710539, 0.09338742, 0.08956818)
cestnik_d3.MAE <- c(0.13515895, 0.12633960, 0.11783315, 0.11031257, 0.10391853, 0.09849350, 0.09374724)
cestnik_d5.MAE <- c(0.1773416, 0.1523275, 0.1361324, 0.1240818, 0.1148086, 0.1073346, 0.1011626)
l_set <- c(1, 2, 3, 4, 5, 6, 7) 
par(pch=21, col="black", bg="white") 
plot(l_set, relfr.MAE, type="n", ylim=c(0, 0.35), xlab="sample size", ylab="average MAE", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(l_set, relfr.MAE, pch=21, col="black", type="o", lty=2, lwd=1) 
lines(l_set, laplace.MAE, pch=22, col="black", type="o", lty=2, lwd=1) 
lines(l_set, piegat.MAE, pch=23, col="black", type="o", lty=2, lwd=1) 
lines(l_set, cestnik.MAE, pch=24, col="black", type="o", lty=2, lwd=1) 
lines(l_set, cestnik_d3.MAE, pch=15, col="black", type="o", lty=2, lwd=1) 
lines(l_set, cestnik_d5.MAE, pch=16, col="black", type="o", lty=2, lwd=1) 

legend("topright", legend=c("relfr", "laplace", "piegat", "cestnik", "cestnik_d.3", "cestnik_d.5"),
       col=c("black", "black", "black", "black", "black", "black"), lty=c(2,2,2,2,2,2), 
       pch=c(21, 22, 23, 24, 15, 16), cex=1.2)


### figure 30, 31, 32
#                                              *         *    11
mmm03 <- c(0.2, 1.1, 2.5, 4.4, 5.6, 6.0, 6.1, 6.5, 6.9, 6.9,  6.9,  6.9, 6.9, 6.5, 6.1, 6.0, 5.6, 4.4, 2.5,  1.1, 0.2)
mmm05 <- c(0.1, 0.7, 1.3, 1.6, 2.2, 2.7, 3.0, 3.1, 3.2, 3.0,  2.9,  3.0, 3.2, 3.1, 3.0, 2.7, 2.2, 1.6, 1.3,  0.7, 0.1)
mmm06 <- c(0.1, 0.5, 1.0, 1.2, 1.6, 2.1, 2.1, 2.3, 2.5, 2.7,  2.9,  2.7, 2.5, 2.3, 2.1, 2.1, 1.6, 1.2, 1.0,  0.5, 0.1)
mmm10 <- c(0.1, 0.3, 0.6, 0.7, 0.9, 1.3, 1.7, 2.3, 2.6, 2.8,  2.9,  2.8, 2.6, 2.3, 1.7, 1.3, 0.9, 0.7, 0.6,  0.3, 0.1)

method_set <- c("prob.cestnik")
distort_line <- 0
distort <- 0.3
l_set <- c(1)
m <- 2
err.list <- list()
err2.list <- list()
err3.list <- list()
cnt <- 0
for (method in method_set) {
  for (l in l_set) {
    cnt <- cnt + 1
    err.new <- array()
    err2.new <- array()
    err3.new <- array()
    for (line in c(1:21)) {
      #compute absolute errors
      err.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      err2.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], mmm03[line], distort)
      err3.new[line] <- prob.sum.MAE("prob.laplace", line, line, distort_line, l, pbase_global[line], m, distort)
      
      print(paste(method, " line:", line, " size:", l, " m:", m, " pa:", pbase_global[line], " distort:", distort, " > "))
    }
    err.list[[cnt]] <- err.new
    err2.list[[cnt]] <- err2.new
    err3.list[[cnt]] <- err3.new
  }
}

mean(err.list[[1]])
mean(err2.list[[1]])
mean(err3.list[[1]])

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=22, col="black") 
plot(p, err.new, type="n", ylim=c(0, 0.4), ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", bg="grey", type="o", lty=2, lwd=1) 
lines(p, err2.list[[1]], pch=16, col="black", bg="grey", type="o", lty=2, lwd=1) 

lines(p, err3.list[[1]], pch=17, col="black", bg="grey", type="o", lty=2, lwd=1) 

legend("topright", legend=c("m=2", "m=optimal" ),
       col=c("black", "black"), lty=c(2,2), 
       pch=c(21, 16), cex=1.2)



### figure 30, 31, 32 *** - with a single m 6,3,1.7 for each distortion 0.3, 0.5, 1.0
#                                              *         *    11
mmm03 <- c(0.2, 1.1, 2.5, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0,  6.0,  6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 2.5, 1.1, 0.2)
mmm05 <- c(0.1, 0.7, 1.3, 1.6, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,  3.0,  3.0, 3.0, 3.0, 3.0, 3.0, 2.2, 1.6, 1.3, 0.7, 0.1)
mmm10 <- c(0.1, 0.3, 0.6, 0.7, 0.9, 1.3, 1.7, 2.0, 2.0, 2.0,  2.0,  2.0, 2.0, 2.0, 1.7, 1.3, 0.9, 0.7, 0.6, 0.3, 0.1)

method_set <- c("prob.cestnik")
distort_line <- 0
distort <- 0.3
l_set <- c(1)
m <- 2
err.list <- list()
err2.list <- list()
err3.list <- list()
cnt <- 0
for (method in method_set) {
  for (l in l_set) {
    cnt <- cnt + 1
    err.new <- array()
    err2.new <- array()
    err3.new <- array()
    for (line in c(1:21)) {
      #compute absolute errors
      err.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], m, distort)
      err2.new[line] <- prob.sum.MAE(method, line, line, distort_line, l, pbase_global[line], mmm03[line], distort)
      err3.new[line] <- prob.sum.MAE("prob.laplace", line, line, distort_line, l, pbase_global[line], m, distort)
      
      print(paste(method, " line:", line, " size:", l, " m:", m, " pa:", pbase_global[line], " distort:", distort, " > "))
    }
    err.list[[cnt]] <- err.new
    err2.list[[cnt]] <- err2.new
    err3.list[[cnt]] <- err3.new
  }
}

mean(err.list[[1]])
mean(err2.list[[1]])
mean(err3.list[[1]])

p <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99) 
par(pch=22, col="black") 
plot(p, err.new, type="n", ylim=c(0, 0.5), ylab="mean absolute error", cex.lab=1.5, cex.axis=1.2, cex.main=1.0, cex.sub=1.0) 

lines(p, err.list[[1]], pch=21, col="black", bg="grey", type="o", lty=2, lwd=1) 
lines(p, err2.list[[1]], pch=16, col="black", bg="grey", type="o", lty=2, lwd=1) 

lines(p, err3.list[[1]], pch=17, col="black", bg="grey", type="o", lty=2, lwd=1) 

legend("topright", legend=c("m=2", "m=optimal", "laplace" ),
       col=c("black", "black", "black"), lty=c(2,2,2), 
       pch=c(21, 16, 17), cex=1.2)




