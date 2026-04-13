library(Rcpp)
library(pg)
library(pisaRT)

sourceCpp('RtIrtGibbs.cpp')
source('RtIrtGibbs_helpers.R')

# 設定 cpp 路徑（fork 子程序需要）
.rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")



fit <- gibbs_rtirt_null(
  Y = as.matrix(pisaW[2:13]),
  logT = as.matrix(pisaW[26:37]),
  
  n_chains = 4,
  n_iter   = 8000,
  n_burnin = 2000,
  n_cores  = 4
  
  
)


posterior_means(fit)

get_dic(fit = fit, loglik_at_mean = )

compute_dic()


#' testing 0413
#' 
library(Rcpp)
.rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")
sourceCpp("RtIrtGibbs.cpp")
source("RtIrtGibbs_helpers.R")

sim <- sim_rtirt(n_subj = 200, n_item = 12, n_feat = 2, seed = 11)

out <- compare_theta_mean_vs_quantile(
  sim,
  n_iter = 2000,
  n_burnin = 1000,
  n_chains = 2,
  q_rt = 0.75,
  n_cores = 2
)

out$theta_compare






simgen <- function(i) sim_rtirt(n_subj=300, n_item=15, n_feat=2, seed=100+i)

out <- run_q_grid_theta_compare(
  sim_generator = simgen,
  q_grid = c(0.25, 0.5, 0.75, 0.85),
  n_rep = 20,
  n_iter = 2000,
  n_burnin = 1000,
  n_chains = 2,
  n_cores = 2
)

out

plot_q_grid_theta(out)


#'' from this
library(Rcpp)
.rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")
sourceCpp("RtIrtGibbs.cpp")
source("RtIrtGibbs_helpers.R")
simgen <- function(i) sim_rtirt(
  n_subj = 300, n_item = 15, n_feat = 2,
  seed = 100 + i,
  rt_error = "t",
  df_t = 4
)

out <- run_q_grid_ability_compare(
  sim_generator = simgen,
  q_grid = c(0.25, 0.5, 0.75, 0.85),
  n_rep = 20,
  n_iter = 2000,
  n_burnin = 1000,
  n_chains = 2,
  n_cores = 2,
  parallel_backend = "psock",   # 關鍵
  quantile_estimate_rt_slope = TRUE,
  quantile_hierarchical_rt_slope = TRUE,
  quantile_adaptive_slope_mh = TRUE
)

out$summary
plot_q_grid_theta(out)
