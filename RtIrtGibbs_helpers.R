# =============================================================================
# RtIrtGibbs_helpers.R
# R-side helpers for the RT-IRT Gibbs sampler (Rcpp/RcppArmadillo backend)
#
# 功能：
#   1. 資料模擬        sim_rtirt_null(), sim_rtirt()
#   2. 多鏈並行執行    sample_rtirt()
#   3. coda 轉換       to_coda()
#   4. 快速診斷        diagnose()
#   5. 後處理          posterior_means(), evaluate_recovery(), get_dic()
#
# 使用流程：
#   Rcpp::sourceCpp("RtIrtGibbs.cpp")
#   source("RtIrtGibbs_helpers.R")
#   .rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")  # 並行時需要
#
#   sim  <- sim_rtirt_null(300, 15)
#   fit  <- sample_rtirt("null", Y=sim$Y, logT=sim$logT, n_chains=4)
#   mc   <- to_coda(fit)
#   summary(mc)
#   coda::gelman.diag(mc)
#   coda::traceplot(mc[, c("a[1]","b[1]"), drop=FALSE])
# =============================================================================

# --- 載入 C++ 端 ---
#' 編譯並載入 Gibbs sampler
#' @param cpp_file RtIrtGibbs.cpp 的路徑
load_rtirt <- function(cpp_file = "RtIrtGibbs.cpp") {
  Rcpp::sourceCpp(cpp_file)
  invisible(NULL)
}

# =============================================================================
# SECTION 1: 資料模擬
# =============================================================================

#' 模擬 Null RT-IRT 資料（無共變數）
#' @param n_subj  受試者數
#' @param n_item  題目數
#' @param true_corr  theta 與 tau 的真實相關
#' @param seed    亂數種子
#' @return list: Y, logT, true_theta, true_tau, true_a, true_b, true_xi, true_sigma_t, Sigma_p
sim_rtirt_null <- function(n_subj = 500, n_item = 20,
                            true_corr = 0.5, seed = 42,
                            rt_error = c("normal", "t", "contaminated"),
                            df_t = 4,
                            contam_prob = 0.10,
                            contam_scale = 3) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("請安裝 MASS 套件: install.packages('MASS')")
  set.seed(seed)
  rt_error <- match.arg(rt_error)
  if (rt_error == "t" && df_t <= 2) stop("df_t 必須 > 2 才有有限變異數")
  if (rt_error == "contaminated" && (contam_prob < 0 || contam_prob > 1))
    stop("contam_prob 必須在 [0,1]")

  true_a  <- pmax(rnorm(n_item, 1.0, 0.2), 0.3)
  true_b  <- rnorm(n_item, 0, 0.5)
  true_xi <- pmax(rnorm(n_item, 4.0, 0.2), 0.5)
  true_st <- sqrt(pmax(rnorm(n_item, 0.3, 0.1), 0.05))

  Sigma_p  <- matrix(c(1, true_corr, true_corr, 0.4), 2, 2)
  persons  <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma_p)
  true_theta <- persons[, 1]
  true_tau   <- persons[, 2]

  eta <- outer(true_theta, true_a) - matrix(true_a * true_b, n_subj, n_item, byrow = TRUE)
  Y   <- matrix(rbinom(n_subj * n_item, 1, plogis(eta)), n_subj, n_item)

  mu_t <- matrix(true_xi, n_subj, n_item, byrow = TRUE) -
    matrix(true_tau, n_subj, n_item)
  sd_mat <- matrix(true_st, n_subj, n_item, byrow = TRUE)
  eps <- switch(
    rt_error,
    "normal" = matrix(rnorm(n_subj * n_item), n_subj, n_item),
    "t" = {
      z <- matrix(rt(n_subj * n_item, df = df_t), n_subj, n_item)
      z / sqrt(df_t / (df_t - 2))  # unit variance for df>2
    },
    "contaminated" = {
      is_contam <- matrix(rbinom(n_subj * n_item, 1, contam_prob), n_subj, n_item)
      z1 <- matrix(rnorm(n_subj * n_item), n_subj, n_item)
      z2 <- matrix(rnorm(n_subj * n_item, sd = contam_scale), n_subj, n_item)
      (1 - is_contam) * z1 + is_contam * z2
    }
  )
  logT <- mu_t + sd_mat * eps

  list(Y = Y, logT = logT,
       true_theta = true_theta, true_tau = true_tau,
       true_a = true_a, true_b = true_b,
       true_xi = true_xi, true_sigma_t = true_st,
       Sigma_p = Sigma_p,
       sim_setting = list(rt_error = rt_error, df_t = df_t,
                          contam_prob = contam_prob, contam_scale = contam_scale))
}

#' 模擬 Structural RT-IRT 資料（含共變數迴歸）
#' @param n_subj   受試者數
#' @param n_item   題目數
#' @param n_feat   共變數個數（不含截距）
#' @param true_beta  真實迴歸係數 (n_feat x 2)，NULL 時隨機生成
#' @param true_corr  殘差相關
#' @param seed     亂數種子
sim_rtirt <- function(n_subj = 500, n_item = 20, n_feat = 2,
                       true_beta = NULL, true_corr = 0.5, seed = 42,
                       rt_error = c("normal", "t", "contaminated"),
                       df_t = 4,
                       contam_prob = 0.10,
                       contam_scale = 3) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("請安裝 MASS 套件: install.packages('MASS')")
  set.seed(seed)
  rt_error <- match.arg(rt_error)
  if (rt_error == "t" && df_t <= 2) stop("df_t 必須 > 2 才有有限變異數")
  if (rt_error == "contaminated" && (contam_prob < 0 || contam_prob > 1))
    stop("contam_prob 必須在 [0,1]")

  if (is.null(true_beta))
    true_beta <- matrix(rnorm(n_feat * 2, 0, 0.5), n_feat, 2)

  X      <- matrix(rnorm(n_subj * n_feat), n_subj, n_feat)
  X_full <- cbind(1, X)

  true_a  <- pmax(rnorm(n_item, 1.0, 0.2), 0.3)
  true_b  <- rnorm(n_item, 0, 0.5)
  true_xi <- pmax(rnorm(n_item, 4.0, 0.2), 0.5)
  true_st <- sqrt(pmax(rnorm(n_item, 0.3, 0.1), 0.05))

  beta_full  <- rbind(c(0, 0), true_beta)
  mu_persons <- X_full %*% beta_full

  Sigma_resid <- matrix(c(1, true_corr, true_corr, 0.4), 2, 2)
  L     <- t(chol(Sigma_resid))
  noise <- t(L %*% matrix(rnorm(2 * n_subj), 2, n_subj))
  persons <- mu_persons + noise

  true_theta <- persons[, 1]
  true_tau   <- persons[, 2]

  eta <- outer(true_theta, true_a) - matrix(true_a * true_b, n_subj, n_item, byrow = TRUE)
  Y   <- matrix(rbinom(n_subj * n_item, 1, plogis(eta)), n_subj, n_item)

  mu_t <- matrix(true_xi, n_subj, n_item, byrow = TRUE) -
    matrix(true_tau, n_subj, n_item)
  sd_mat <- matrix(true_st, n_subj, n_item, byrow = TRUE)
  eps <- switch(
    rt_error,
    "normal" = matrix(rnorm(n_subj * n_item), n_subj, n_item),
    "t" = {
      z <- matrix(rt(n_subj * n_item, df = df_t), n_subj, n_item)
      z / sqrt(df_t / (df_t - 2))
    },
    "contaminated" = {
      is_contam <- matrix(rbinom(n_subj * n_item, 1, contam_prob), n_subj, n_item)
      z1 <- matrix(rnorm(n_subj * n_item), n_subj, n_item)
      z2 <- matrix(rnorm(n_subj * n_item, sd = contam_scale), n_subj, n_item)
      (1 - is_contam) * z1 + is_contam * z2
    }
  )
  logT <- mu_t + sd_mat * eps

  list(Y = Y, logT = logT, X = X,
       true_theta = true_theta, true_tau = true_tau,
       true_a = true_a, true_b = true_b,
       true_xi = true_xi, true_sigma_t = true_st,
       true_beta = true_beta, Sigma_resid = Sigma_resid,
       sim_setting = list(rt_error = rt_error, df_t = df_t,
                          contam_prob = contam_prob, contam_scale = contam_scale))
}

# =============================================================================
# SECTION 2: 多鏈並行執行
# =============================================================================

#' 執行多鏈 RT-IRT Gibbs sampler（支援自動並行）
#'
#' 每條鏈使用獨立亂數種子，透過 parallel::mclapply 在 Linux/macOS 上並行。
#' Windows 上 mclapply 不支援 fork，自動退回序列執行。
#'
#' @section 並行注意事項:
#' 因為 fork 子程序不繼承已載入的 shared library，每個子程序需要重新
#' sourceCpp()。請在呼叫前設定：
#'   .rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")
#' 若沒設定，會自動在工作目錄找 "RtIrtGibbs.cpp"。
#'
#' @param model   模型："null" | "structural" | "quantile" | "ml_irt"
#' @param Y       n_subj x n_item 作答矩陣 (0/1)
#' @param logT    n_subj x n_item log 反應時間（"null"/"structural"/"quantile" 需要）
#' @param X_cov   n_subj x n_feat 共變數（"structural"/"quantile"/"ml_irt" 需要）
#' @param n_chains  鏈數，預設 4
#' @param n_iter    每鏈總迭代數，預設 5000
#' @param n_burnin  每鏈 burnin，預設 n_iter/2
#' @param n_cores   CPU 核心數，NULL = 自動偵測（最多 n_chains）
#' @param base_seed 基礎種子，第 i 鏈種子 = base_seed + i
#' @param verbose   TRUE = 顯示進度條（序列模式）；並行模式下自動關閉
#' @param ...       傳遞給底層 gibbs_*() 的其他參數
#'                  例如: one_pl=TRUE, fix_intercept=FALSE, q_ra=0.25, mu_xi=4.0
#'
#' @return class "rtirt_chains" 的 list：
#'   $chains   各鏈原始輸出的 list
#'   $model    模型名稱
#'   $n_chains, $n_iter, $n_burnin, $n_save
sample_rtirt <- function(model = c("null", "structural", "quantile", "cross_quantile", "ml_irt"),
                          Y,
                          logT    = NULL,
                          X_cov   = NULL,
                          n_chains  = 4L,
                          n_iter    = 5000L,
                          n_burnin  = NULL,
                          n_cores   = NULL,
                          parallel_backend = c("auto", "multicore", "psock", "sequential"),
                          base_seed = 1234L,
                          laplace_init = FALSE,
                          verbose   = TRUE,
                          ...) {
  model <- match.arg(model)
  parallel_backend <- match.arg(parallel_backend)
  extra_args <- list(...)

  # quantile path now supports hierarchical_rt_slope as well

  # --- 驗證 ---
  if (model != "ml_irt" && is.null(logT))
    stop("model='", model, "' 需要提供 logT 參數")
  if (model %in% c("structural", "quantile", "cross_quantile", "ml_irt") && is.null(X_cov))
    stop("model='", model, "' 需要提供 X_cov 參數")
  if (is.null(n_burnin))
    n_burnin <- as.integer(n_iter / 2)

  n_chains <- as.integer(n_chains)
  n_iter   <- as.integer(n_iter)
  n_burnin <- as.integer(n_burnin)

  # --- 決定並行策略 ---
  on_windows <- .Platform$OS.type == "windows"
  in_rstudio <- nzchar(Sys.getenv("RSTUDIO"))
  if (is.null(n_cores))
    n_cores <- min(n_chains,
                   max(1L, parallel::detectCores(logical = FALSE) - 1L))

  backend <- parallel_backend
  if (backend == "auto") {
    if (on_windows) {
      backend <- if (n_cores > 1L && n_chains > 1L) "psock" else "sequential"
    } else if (in_rstudio) {
      backend <- if (n_cores > 1L && n_chains > 1L) "psock" else "sequential"
    } else {
      backend <- if (n_cores > 1L && n_chains > 1L) "multicore" else "sequential"
    }
  }
  use_parallel <- backend %in% c("multicore", "psock") && n_cores > 1L && n_chains > 1L

  if (use_parallel) {
    message(sprintf(
      "[sample_rtirt] model=%s | %d chains x %d iter (burnin=%d) | %d cores (%s)",
      model, n_chains, n_iter, n_burnin, n_cores, backend))
  } else {
    if (on_windows && n_chains > 1 && parallel_backend %in% c("auto", "multicore"))
      message("[sample_rtirt] Windows 不支援 fork；改用 psock/sequential")
    if (in_rstudio && n_chains > 1 && parallel_backend == "auto")
      message("[sample_rtirt] 偵測到 RStudio；auto 模式避免 fork，改用 psock/sequential")
    message(sprintf(
      "[sample_rtirt] model=%s | %d chains x %d iter (burnin=%d) | sequential",
      model, n_chains, n_iter, n_burnin))
  }

  laplace_ab_init <- function(Y, theta0 = NULL, mu_a = 1, sigma_a = 1,
                              mu_b = 0, sigma_b = 1) {
    n <- nrow(Y); p <- ncol(Y)
    if (is.null(theta0)) theta0 <- qnorm((rowMeans(Y) + 0.5) / 2)
    theta0 <- as.numeric(scale(theta0, scale = FALSE))
    a <- rep(1, p); b <- rep(0, p)
    prop_sd_loga <- rep(0.08, p); prop_sd_b <- rep(0.12, p)

    for (j in seq_len(p)) {
      yj <- Y[, j]
      nll <- function(par) {
        loga <- par[1]
        bj <- par[2]
        aj <- exp(loga)
        eta <- aj * (theta0 - bj)
        ll <- sum(yj * eta - log1p(exp(eta)))
        lp <- dnorm(aj, mu_a, sigma_a, log = TRUE) + loga +
          dnorm(bj, mu_b, sigma_b, log = TRUE)
        -(ll + lp)
      }
      fit <- optim(c(0, 0), nll, method = "BFGS", hessian = TRUE,
                   control = list(maxit = 200, reltol = 1e-8))
      loga_hat <- fit$par[1]
      b_hat <- fit$par[2]
      a[j] <- exp(loga_hat)
      b[j] <- b_hat
      if (!is.null(fit$hessian) && all(is.finite(fit$hessian))) {
        v <- tryCatch(solve(fit$hessian), error = function(e) diag(c(0.0064, 0.0144)))
        prop_sd_loga[j] <- sqrt(max(v[1, 1], 1e-4))
        prop_sd_b[j] <- sqrt(max(v[2, 2], 1e-4))
      }
    }
    list(a = a, b = b,
         prop_sd_loga = median(prop_sd_loga),
         prop_sd_b = median(prop_sd_b))
  }

  # --- 單鏈執行函數 ---
  run_one_chain <- function(chain_id) {
    set.seed(base_seed + chain_id)

    # fork 子程序：必須先載入 pg，再 sourceCpp
    # 順序很重要：pg 要先 library()，sourceCpp 才找得到 <pg.h> 的 symbols
    if (use_parallel) {
      cpp_path <- tryCatch(
        get(".rtirt_cpp_path", envir = .GlobalEnv),
        error = function(e) "RtIrtGibbs.cpp"
      )
      if (!file.exists(cpp_path)) {
        cpp_fallback <- normalizePath("RtIrtGibbs.cpp", mustWork = FALSE)
        if (file.exists(cpp_fallback)) {
          cpp_path <- cpp_fallback
        } else {
          stop(
            "找不到 RtIrtGibbs.cpp。請在目前工作目錄放置該檔案，",
            "或設定 .rtirt_cpp_path <- normalizePath('RtIrtGibbs.cpp')"
          )
        }
      }
      # 1. 確保 pg 在子程序裡是載入的
      if (!requireNamespace("pg", quietly = TRUE))
        stop("pg 套件未安裝，請先 R CMD INSTALL /path/to/pg")
      library(pg, quietly = TRUE)

      # 2. 重新載入 C++：
      #    先用預設 cache（快）；若遇到 zero-byte 類競態錯誤，再 fallback 到隔離 cache（穩）
      fast_res <- try({
        Rcpp::sourceCpp(cpp_path, verbose = FALSE, rebuild = FALSE)
      }, silent = TRUE)
      if (inherits(fast_res, "try-error")) {
        cache_dir <- file.path(tempdir(), sprintf("rcpp-cache-%d-%d", Sys.getpid(), chain_id))
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        safe_res <- try({
          Rcpp::sourceCpp(cpp_path, verbose = FALSE, rebuild = FALSE, cacheDir = cache_dir)
        }, silent = TRUE)
        if (inherits(safe_res, "try-error")) {
          stop("sourceCpp 載入失敗（子程序）: ",
               "fast=", as.character(fast_res), " | safe=", as.character(safe_res))
        }
      }
    }

    # 函數名稱在 sourceCpp 之後才解析（用字串 + get() 取函數物件）
    fn_name <- switch(model,
      "null"       = "gibbs_rtirt_null",
      "structural" = "gibbs_rtirt",
      "quantile"   = "gibbs_rtirt_quantile",
      "cross_quantile" = "gibbs_rtirt_cross_quantile",
      "ml_irt"     = "gibbs_ml_irt"
    )
    # 在 global env 找函數（sourceCpp 把函數放在 global env）
    fn <- tryCatch(
      get(fn_name, envir = .GlobalEnv),
      error = function(e)
        stop(sprintf("找不到 %s，請確認 sourceCpp() 已執行", fn_name))
    )

    data_args <- switch(model,
      "null"       = list(Y = Y, logT = logT),
      "structural" = list(Y = Y, logT = logT, X_cov = X_cov),
      "quantile"   = list(Y = Y, logT = logT, X_cov = X_cov),
      "cross_quantile" = list(Y = Y, logT = logT, X_cov = X_cov),
      "ml_irt"     = list(Y = Y, X_cov = X_cov)
    )

    args_use <- extra_args
    if (isTRUE(laplace_init) && model %in% c("null", "structural", "ml_irt")) {
      la <- laplace_ab_init(Y)
      args_use$use_ab_mh <- TRUE
      args_use$adaptive_mh <- TRUE
      args_use$prop_sd_loga <- la$prop_sd_loga
      args_use$prop_sd_b <- la$prop_sd_b
    }

    do.call(fn, c(
      data_args,
      list(n_iter   = n_iter,
           n_burnin = n_burnin,
           verbose  = verbose && !use_parallel,
           chain_id = chain_id),
      args_use
    ))
  }

  # --- 執行 ---
  t_start <- proc.time()

  if (use_parallel && backend == "multicore") {
    chains <- parallel::mclapply(
      seq_len(n_chains),
      run_one_chain,
      mc.cores    = n_cores,
      mc.set.seed = FALSE,   # 各鏈在 run_one_chain 裡自行 set.seed
      mc.preschedule = FALSE # 減少單一子程序錯誤影響全部工作
    )
    err_idx <- which(vapply(chains, inherits, logical(1L), "try-error"))
    if (length(err_idx) > 0) {
      # try-error 物件的錯誤訊息存在 attr(x, "condition")$message 或直接轉 character
      msgs <- vapply(chains[err_idx], function(e) {
        cond <- attr(e, "condition")
        if (!is.null(cond)) conditionMessage(cond)
        else as.character(e)
      }, character(1L))
      stop("以下鏈發生錯誤:\n",
           paste0("  Chain ", err_idx, ": ", msgs, collapse = "\n"))
    }
    bad_idx <- which(vapply(chains, function(ch) {
      is.null(ch) || !is.list(ch) || is.null(names(ch))
    }, logical(1L)))
    if (length(bad_idx) > 0) {
      stop("以下鏈沒有回傳有效結果（可能是子程序中斷）:\n",
           paste0("  Chain ", bad_idx, collapse = "\n"),
           "\n建議降低 n_cores、縮短 n_iter 測試，或先以 n_chains=1 驗證設定。")
    }
  } else if (use_parallel && backend == "psock") {
    cl <- parallel::makeCluster(min(n_cores, n_chains))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    cpp_path <- tryCatch(
      get(".rtirt_cpp_path", envir = .GlobalEnv),
      error = function(e) "RtIrtGibbs.cpp"
    )
    if (!file.exists(cpp_path)) {
      cpp_fallback <- normalizePath("RtIrtGibbs.cpp", mustWork = FALSE)
      if (file.exists(cpp_fallback)) cpp_path <- cpp_fallback
    }

    parallel::clusterExport(cl, varlist = c(
      "model", "Y", "logT", "X_cov", "n_iter", "n_burnin",
      "base_seed", "verbose", "laplace_init", "extra_args",
      "laplace_ab_init", "cpp_path"
    ), envir = environment())

    parallel::clusterEvalQ(cl, {
      if (!requireNamespace("pg", quietly = TRUE)) {
        stop("pg 套件未安裝，請先 R CMD INSTALL /path/to/pg")
      }
      library(pg, quietly = TRUE)
      suppressMessages(Rcpp::sourceCpp(cpp_path, verbose = FALSE, rebuild = FALSE))
      NULL
    })

    worker_run_chain <- function(chain_id) {
      set.seed(base_seed + chain_id)
      fn_name <- switch(model,
                        "null" = "gibbs_rtirt_null",
                        "structural" = "gibbs_rtirt",
                        "quantile" = "gibbs_rtirt_quantile",
                        "cross_quantile" = "gibbs_rtirt_cross_quantile",
                        "ml_irt" = "gibbs_ml_irt")
      fn <- get(fn_name, envir = .GlobalEnv)
      data_args <- switch(model,
                          "null" = list(Y = Y, logT = logT),
                          "structural" = list(Y = Y, logT = logT, X_cov = X_cov),
                          "quantile" = list(Y = Y, logT = logT, X_cov = X_cov),
                          "cross_quantile" = list(Y = Y, logT = logT, X_cov = X_cov),
                          "ml_irt" = list(Y = Y, X_cov = X_cov))
      args_use <- extra_args
      if (isTRUE(laplace_init) && model %in% c("null", "structural", "ml_irt")) {
        la <- laplace_ab_init(Y)
        args_use$use_ab_mh <- TRUE
        args_use$adaptive_mh <- TRUE
        args_use$prop_sd_loga <- la$prop_sd_loga
        args_use$prop_sd_b <- la$prop_sd_b
      }
      do.call(fn, c(
        data_args,
        list(n_iter = n_iter,
             n_burnin = n_burnin,
             verbose = FALSE,
             chain_id = chain_id),
        args_use
      ))
    }

    parallel::clusterExport(cl, varlist = c("worker_run_chain"), envir = environment())
    chains <- parallel::parLapply(cl, seq_len(n_chains), worker_run_chain)
    bad_idx <- which(vapply(chains, function(ch) {
      is.null(ch) || !is.list(ch) || is.null(names(ch))
    }, logical(1L)))
    if (length(bad_idx) > 0) {
      stop("以下鏈沒有回傳有效結果（psock worker 可能失敗）:\n",
           paste0("  Chain ", bad_idx, collapse = "\n"))
    }
  } else {
    chains <- vector("list", n_chains)
    for (cid in seq_len(n_chains)) {
      if (verbose) cat(sprintf("\n--- Chain %d / %d ---\n", cid, n_chains))
      chains[[cid]] <- run_one_chain(cid)
    }
  }

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  message(sprintf("[sample_rtirt] 完成，耗時 %.1f 秒", elapsed))

  structure(
    list(chains   = chains,
         model    = model,
         n_iter   = n_iter,
         n_burnin = n_burnin,
         n_chains = n_chains,
         n_save   = n_iter - n_burnin),
    class = "rtirt_chains"
  )
}

# --- print method ---
print.rtirt_chains <- function(x, ...) {
  cat(sprintf(
    "rtirt_chains  [model=%s | chains=%d | n_save=%d per chain]\n",
    x$model, x$n_chains, x$n_save
  ))
  params <- setdiff(names(x$chains[[1]]), "loglik")
  cat("Stored parameters:\n")
  for (p in params) {
    m <- x$chains[[1]][[p]]
    dim_str <- if (is.matrix(m)) sprintf("%d columns", ncol(m)) else "1 value"
    cat(sprintf("  $%-10s  %s\n", p, dim_str))
  }
  invisible(x)
}

# =============================================================================
# SECTION 3: coda 轉換
# =============================================================================

#' 將 rtirt_chains 或單鏈 list 轉換為 coda::mcmc.list
#'
#' 把各鏈的所有參數矩陣橫向拼接，形成帶有語意欄位名稱的寬矩陣，
#' 再包裝為 coda::mcmc.list。轉換後可直接使用所有 coda 診斷工具。
#'
#' @section 欄位命名規則:
#' \itemize{
#'   \item 一般參數：\code{theta[1]}, \code{theta[2]}, ..., \code{a[1]}, ...
#'   \item 2x2 person covariance：\code{Sigma_p[1,1]}, \code{Sigma_p[2,1]},
#'         \code{Sigma_p[1,2]}, \code{Sigma_p[2,2]}（column-major 順序）
#'   \item 迴歸係數：\code{beta[1,1]}, \code{beta[2,1]}, ... （行=covariate，列=outcome）
#' }
#'
#' @param fit    sample_rtirt() 的輸出，或單條鏈的 list
#' @param params 要轉換的欄位名稱，NULL = 所有非 loglik 欄位
#' @param thin   thinning 間隔，預設 1
#'
#' @return coda::mcmc.list，可用：
#'   \code{summary()}, \code{coda::gelman.diag()}, \code{coda::effectiveSize()},
#'   \code{coda::traceplot()}, \code{coda::densplot()}, \code{coda::autocorr.plot()}
to_coda <- function(fit, params = NULL, thin = 1L) {
  if (!requireNamespace("coda", quietly = TRUE))
    stop("請安裝 coda 套件: install.packages('coda')")

  # 支援單鏈 list 直接傳入
  if (!inherits(fit, "rtirt_chains")) {
    if (is.list(fit) && !is.null(names(fit))) {
      fit <- structure(
        list(chains  = list(fit),
             model   = "unknown",
             n_chains = 1L,
             n_save  = if (is.matrix(fit[[1]])) nrow(fit[[1]]) else length(fit[[1]])),
        class = "rtirt_chains"
      )
    } else {
      stop("fit 必須是 sample_rtirt() 的輸出，或單條鏈的 named list")
    }
  }

  all_params <- names(fit$chains[[1]])
  if (is.null(params)) {
    preferred <- c(
      "person_ability", "person_speed",
      "item_discrimination", "item_difficulty",
      "item_time_intensity", "rt_resid_sd", "rt_speed_loading",
      "person_regression", "person_cov",
      "accept_ab", "accept_slope", "slope_hyper", "log_likelihood"
    )
    preferred <- intersect(preferred, all_params)
    if (length(preferred) > 0) {
      params <- setdiff(preferred, c("log_likelihood"))
    } else {
      params <- setdiff(all_params, "loglik")
    }
  } else {
    missing_p <- setdiff(params, all_params)
    if (length(missing_p) > 0)
      warning("找不到欄位，已忽略: ", paste(missing_p, collapse = ", "))
    params <- intersect(params, all_params)
  }
  if (length(params) == 0) stop("沒有可轉換的參數")

  # --- 欄位命名 ---
  make_colnames <- function(pname, ncols) {
    if (ncols == 1) return(pname)

    if (pname == "Sigma_p") {
      # 2x2，column-major: [1,1],[2,1],[1,2],[2,2]
      return(c("Sigma_p[1,1]", "Sigma_p[2,1]", "Sigma_p[1,2]", "Sigma_p[2,2]"))
    }
    if (pname == "person_cov") {
      return(c("person_cov[1,1]", "person_cov[2,1]", "person_cov[1,2]", "person_cov[2,2]"))
    }
    if (pname == "beta") {
      # stored as vec(beta) where beta is q x 2, column-major
      # ncols = q * 2; first q cols = column 1 (ra), next q = column 2 (rt)
      q <- ncols / 2
      c(sprintf("beta[%d,ra]", seq_len(q)),
        sprintf("beta[%d,rt]", seq_len(q)))
    } else if (pname == "person_regression") {
      q <- ncols / 2
      c(sprintf("person_regression[%d,ability]", seq_len(q)),
        sprintf("person_regression[%d,speed]", seq_len(q)))
    } else if (pname == "slope_hyper" && ncols == 2) {
      c("slope_hyper[mu_log_s]", "slope_hyper[sigma2_log_s]")
    } else {
      sprintf("%s[%d]", pname, seq_len(ncols))
    }
  }

  # 取第一條鏈計算每個參數的列數，建立全域欄位名稱
  col_names <- unlist(lapply(params, function(p) {
    m <- fit$chains[[1]][[p]]
    nc <- if (is.matrix(m)) ncol(m) else 1L
    make_colnames(p, nc)
  }))

  # --- 每條鏈轉成矩陣 ---
  chain_to_matrix <- function(ch) {
    parts <- lapply(params, function(p) {
      x <- ch[[p]]
      if (!is.matrix(x)) matrix(x, ncol = 1L) else x
    })
    m <- do.call(cbind, parts)
    colnames(m) <- col_names
    m
  }

  mcmc_list <- lapply(seq_len(fit$n_chains), function(i) {
    mat <- chain_to_matrix(fit$chains[[i]])
    if (thin > 1L) {
      keep <- seq(1L, nrow(mat), by = thin)
      mat  <- mat[keep, , drop = FALSE]
    }
    coda::mcmc(mat, thin = thin)
  })

  coda::mcmc.list(mcmc_list)
}

#' 快速印出 Gelman-Rubin 診斷與有效樣本數
#'
#' 自動篩掉高維度的 theta/tau（避免輸出爆炸），
#' 只顯示 item 參數和結構參數。
#'
#' @param fit     sample_rtirt() 的輸出
#' @param params  NULL = 自動選（排除 theta/tau/loglik）
#' @param n_show  每組參數只顯示前幾個，預設 5
diagnose <- function(fit, params = NULL, n_show = 5L) {
  if (!requireNamespace("coda", quietly = TRUE))
    stop("請安裝 coda 套件: install.packages('coda')")

  if (is.null(params))
    params <- setdiff(names(fit$chains[[1]]),
                      c("loglik", "log_likelihood", "theta", "tau",
                        "person_ability", "person_speed"))  # 排除大維度

  mc  <- to_coda(fit, params = params)
  nv  <- coda::nvar(mc)
  idx <- seq_len(min(n_show * length(params), nv))
  mc_sub <- mc[, idx, drop = FALSE]

  cat("=== Gelman-Rubin R-hat（目標 < 1.1） ===\n")
  if (fit$n_chains >= 2L) {
    gr <- coda::gelman.diag(mc_sub, multivariate = FALSE)
    print(round(gr$psrf, 4))
  } else {
    cat("（需要 >= 2 條鏈）\n")
  }

  cat("\n=== 有效樣本數（ESS） ===\n")
  print(round(coda::effectiveSize(mc_sub)))

  mh_fields <- c("accept_ab", "accept_slope")
  has_mh <- any(mh_fields %in% names(fit$chains[[1]]))
  if (has_mh) {
    cat("\n=== MH 診斷（鏈別） ===\n")
    print(diagnose_ab_mh(fit))
  }

  invisible(NULL)
}

#' Diagnose AB Metropolis-Hastings behavior
#' @param fit sample_rtirt() output
#' @return data.frame with chain-level MH diagnostics
diagnose_ab_mh <- function(fit) {
  if (!inherits(fit, "rtirt_chains")) stop("fit 必須是 sample_rtirt() 的輸出")
  out <- lapply(seq_along(fit$chains), function(i) {
    ch <- fit$chains[[i]]
    acc <- ch$accept_ab
    if (is.null(acc)) acc <- NA_real_
    data.frame(
      chain = i,
      mean_accept_ab = mean(acc, na.rm = TRUE),
      p05_accept_ab = stats::quantile(acc, 0.05, na.rm = TRUE, names = FALSE),
      p95_accept_ab = stats::quantile(acc, 0.95, na.rm = TRUE, names = FALSE),
      mean_accept_slope = if (!is.null(ch$accept_slope)) mean(ch$accept_slope, na.rm = TRUE) else NA_real_,
      final_prop_sd_loga = if (!is.null(ch$final_prop_sd_loga)) ch$final_prop_sd_loga else NA_real_,
      final_prop_sd_b = if (!is.null(ch$final_prop_sd_b)) ch$final_prop_sd_b else NA_real_,
      final_prop_sd_log_slope = if (!is.null(ch$final_prop_sd_log_slope)) ch$final_prop_sd_log_slope else NA_real_
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Quantile-regression applicability note
#' @return character vector describing when quantile RT-IRT is suitable
quantile_applicability_note <- function() {
  c(
    "適用：RT 殘差明顯偏態/厚尾，或你關心慢速尾端（如 q=0.75, 0.85）的個體差異。",
    "適用：研究問題是條件分位數效應，而非平均反應時間效應。",
    "不建議：樣本/題數太小時（q 高或低）會讓尾端訊號不足，混合變慢。",
    "目前程式：quantile 路徑尚未接上 RT slope 隨機斜率，建議先在 mean RT 模型驗證 slope。"
  )
}

#' Compare theta recovery: mean RT vs quantile RT
#' @param sim simulation list from sim_rtirt() or sim_rtirt_null()
#' @param n_iter total iterations per chain
#' @param n_burnin burnin iterations
#' @param n_chains number of chains
#' @param q_rt quantile level for RT model
#' @param n_cores cores for parallel chains
#' @return list with fit_mean, fit_quantile, and theta comparison table
compare_theta_mean_vs_quantile <- function(sim,
                                           n_iter = 2000L,
                                           n_burnin = 1000L,
                                           n_chains = 2L,
                                           q_rt = 0.75,
                                           n_cores = NULL,
                                           parallel_backend = "auto",
                                           quantile_estimate_rt_slope = TRUE,
                                           quantile_hierarchical_rt_slope = FALSE,
                                           quantile_adaptive_slope_mh = FALSE) {
  if (is.null(sim$true_theta)) stop("sim 必須包含 true_theta")
  is_struct <- !is.null(sim$X)
  model_mean <- if (is_struct) "structural" else "null"

  common_args <- list(
    Y = sim$Y,
    logT = sim$logT,
    n_iter = as.integer(n_iter),
    n_burnin = as.integer(n_burnin),
    n_chains = as.integer(n_chains),
    n_cores = n_cores,
    parallel_backend = parallel_backend,
    verbose = FALSE
  )
  if (is_struct) common_args$X_cov <- sim$X

  fit_mean <- do.call(sample_rtirt, c(list(model = model_mean), common_args))
  fit_q <- do.call(sample_rtirt, c(
    list(
      model = "quantile",
      q_rt = q_rt,
      estimate_rt_slope = quantile_estimate_rt_slope,
      hierarchical_rt_slope = quantile_hierarchical_rt_slope,
      adaptive_slope_mh = quantile_adaptive_slope_mh
    ),
    common_args
  ))

  pm_mean <- posterior_means(fit_mean)
  pm_q <- posterior_means(fit_q)
  truth <- sim$true_theta

  res <- rbind(
    data.frame(
      model = "mean_rt",
      rmse_theta = rmse(pm_mean$theta, truth),
      bias_theta = bias(pm_mean$theta, truth),
      cor_theta = suppressWarnings(stats::cor(pm_mean$theta, truth))
    ),
    data.frame(
      model = sprintf("quantile_rt_q%.2f", q_rt),
      rmse_theta = rmse(pm_q$theta, truth),
      bias_theta = bias(pm_q$theta, truth),
      cor_theta = suppressWarnings(stats::cor(pm_q$theta, truth))
    )
  )
  rownames(res) <- NULL

  list(
    fit_mean = fit_mean,
    fit_quantile = fit_q,
    theta_compare = res
  )
}

#' Run theta comparison over q-grid (optionally repeated)
#' @param sim_generator function(rep_id) returning simulation list
#' @param q_grid numeric vector of quantile levels
#' @param n_rep number of replications
#' @param n_iter total iterations
#' @param n_burnin burnin iterations
#' @param n_chains chains
#' @param n_cores cores
#' @return list(raw, summary)
run_q_grid_theta_compare <- function(sim_generator,
                                     q_grid = c(0.25, 0.5, 0.75, 0.85),
                                     n_rep = 10L,
                                     n_iter = 2000L,
                                     n_burnin = 1000L,
                                     n_chains = 2L,
                                     n_cores = NULL,
                                     parallel_backend = "auto",
                                     quantile_estimate_rt_slope = TRUE,
                                     quantile_hierarchical_rt_slope = FALSE,
                                     quantile_adaptive_slope_mh = FALSE) {
  if (!is.function(sim_generator)) stop("sim_generator 必須是函數")
  q_grid <- as.numeric(q_grid)
  n_rep <- as.integer(n_rep)

  rows <- vector("list", length = n_rep * (length(q_grid) + 1L))
  idx <- 1L

  for (r in seq_len(n_rep)) {
    sim <- sim_generator(r)
    if (is.null(sim$true_theta)) stop("sim_generator 回傳的資料必須包含 true_theta")
    is_struct <- !is.null(sim$X)
    model_mean <- if (is_struct) "structural" else "null"

    common_args <- list(
      Y = sim$Y,
      logT = sim$logT,
      n_iter = as.integer(n_iter),
      n_burnin = as.integer(n_burnin),
      n_chains = as.integer(n_chains),
      n_cores = n_cores,
      parallel_backend = parallel_backend,
      verbose = FALSE
    )
    if (is_struct) common_args$X_cov <- sim$X

    fit_mean <- do.call(sample_rtirt, c(list(model = model_mean), common_args))
    pm_mean <- posterior_means(fit_mean)
    rows[[idx]] <- data.frame(
      rep = r,
      model = "mean_rt",
      q_rt = NA_real_,
      rmse_theta = rmse(pm_mean$theta, sim$true_theta),
      bias_theta = bias(pm_mean$theta, sim$true_theta),
      cor_theta = suppressWarnings(stats::cor(pm_mean$theta, sim$true_theta))
    )
    idx <- idx + 1L

    for (qv in q_grid) {
      fit_q <- do.call(sample_rtirt, c(
        list(
          model = "quantile",
          q_rt = qv,
          estimate_rt_slope = quantile_estimate_rt_slope,
          hierarchical_rt_slope = quantile_hierarchical_rt_slope,
          adaptive_slope_mh = quantile_adaptive_slope_mh
        ),
        common_args
      ))
      pm_q <- posterior_means(fit_q)
      rows[[idx]] <- data.frame(
        rep = r,
        model = sprintf("quantile_rt_q%.2f", qv),
        q_rt = qv,
        rmse_theta = rmse(pm_q$theta, sim$true_theta),
        bias_theta = bias(pm_q$theta, sim$true_theta),
        cor_theta = suppressWarnings(stats::cor(pm_q$theta, sim$true_theta))
      )
      idx <- idx + 1L
    }
  }

  raw <- do.call(rbind, rows)
  rownames(raw) <- NULL

  by_keys <- list(model = raw$model, q_rt = raw$q_rt)
  mean_tab <- aggregate(
    raw[, c("rmse_theta", "bias_theta", "cor_theta")],
    by = by_keys,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  sd_tab <- aggregate(
    raw[, c("rmse_theta", "bias_theta", "cor_theta")],
    by = by_keys,
    FUN = function(x) stats::sd(x, na.rm = TRUE)
  )
  names(mean_tab)[3:5] <- paste0(names(mean_tab)[3:5], "_mean")
  names(sd_tab)[3:5] <- paste0(names(sd_tab)[3:5], "_sd")
  summary_out <- merge(mean_tab, sd_tab, by = c("model", "q_rt"), sort = FALSE)
  summary_out <- summary_out[order(summary_out$q_rt, na.last = TRUE), ]
  rownames(summary_out) <- NULL

  list(raw = raw, summary = summary_out)
}

#' Plot theta performance across q-grid
#' @param out output from run_q_grid_theta_compare()
#' @return invisible(NULL)
plot_q_grid_theta <- function(out) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Case A: q-grid output
  if (!is.null(out$summary)) {
    s <- out$summary
    mean_row <- s[s$model == "mean_rt", , drop = FALSE]
    q_rows <- s[grepl("^quantile_rt_q", s$model), , drop = FALSE]
    if (nrow(q_rows) == 0) stop("summary 中找不到 quantile 結果")
    q_rows <- q_rows[order(q_rows$q_rt), , drop = FALSE]

    y_rmse <- q_rows$rmse_theta_mean
    x_q <- q_rows$q_rt
    graphics::plot(x_q, y_rmse, type = "b", pch = 19, lwd = 2,
                   xlab = "q_rt", ylab = "Theta RMSE",
                   main = "Theta RMSE vs q_rt")
    if (nrow(mean_row) > 0 && is.finite(mean_row$rmse_theta_mean[1])) {
      graphics::abline(h = mean_row$rmse_theta_mean[1], lty = 2, lwd = 2, col = "gray40")
      graphics::legend("topright",
                       legend = c("quantile", "mean baseline"),
                       lty = c(1, 2), pch = c(19, NA), bty = "n")
    }
    if (!all(is.na(q_rows$rmse_theta_sd))) {
      graphics::arrows(x_q, y_rmse - q_rows$rmse_theta_sd,
                       x_q, y_rmse + q_rows$rmse_theta_sd,
                       angle = 90, code = 3, length = 0.05, col = "gray50")
    }

    y_cor <- q_rows$cor_theta_mean
    graphics::plot(x_q, y_cor, type = "b", pch = 19, lwd = 2,
                   xlab = "q_rt", ylab = "cor(theta_hat, theta_true)",
                   main = "Theta Correlation vs q_rt", ylim = c(min(0.0, min(y_cor, na.rm = TRUE)), 1))
    if (nrow(mean_row) > 0 && is.finite(mean_row$cor_theta_mean[1])) {
      graphics::abline(h = mean_row$cor_theta_mean[1], lty = 2, lwd = 2, col = "gray40")
      graphics::legend("bottomright",
                       legend = c("quantile", "mean baseline"),
                       lty = c(1, 2), pch = c(19, NA), bty = "n")
    }
    if (!all(is.na(q_rows$cor_theta_sd))) {
      graphics::arrows(x_q, y_cor - q_rows$cor_theta_sd,
                       x_q, y_cor + q_rows$cor_theta_sd,
                       angle = 90, code = 3, length = 0.05, col = "gray50")
    }
    return(invisible(NULL))
  }

  # Case B: single compare output
  if (!is.null(out$theta_compare)) {
    d <- out$theta_compare
    x <- seq_len(nrow(d))
    labs <- d$model

    graphics::plot(x, d$rmse_theta, pch = 19, xaxt = "n",
                   xlab = "model", ylab = "Theta RMSE",
                   main = "Theta RMSE (mean vs quantile)")
    graphics::axis(1, at = x, labels = labs, las = 2, cex.axis = 0.8)
    graphics::lines(x, d$rmse_theta, lwd = 2)

    graphics::plot(x, d$cor_theta, pch = 19, xaxt = "n",
                   xlab = "model", ylab = "cor(theta_hat, theta_true)",
                   main = "Theta Correlation (mean vs quantile)", ylim = c(0, 1))
    graphics::axis(1, at = x, labels = labs, las = 2, cex.axis = 0.8)
    graphics::lines(x, d$cor_theta, lwd = 2)
    return(invisible(NULL))
  }

  stop("out 需來自 compare_theta_mean_vs_quantile() 或 run_q_grid_theta_compare()")
}

# =============================================================================
# SECTION 7: Descriptive-name helper aliases (backward compatible)
# =============================================================================

#' Simulate joint RT-IRT null model data
sim_joint_rt_irt_null <- function(...) sim_rtirt_null(...)

#' Simulate joint RT-IRT structural model data
sim_joint_rt_irt <- function(...) sim_rtirt(...)

#' Run multi-chain joint RT-IRT sampler
sample_joint_rt_irt <- function(...) sample_rtirt(...)

#' Convert joint RT-IRT chains to coda format
to_coda_joint_rt_irt <- function(...) to_coda(...)

#' Quick MCMC diagnostics for joint RT-IRT
diagnose_joint_rt_irt <- function(...) diagnose(...)

#' Posterior means for joint RT-IRT outputs
posterior_means_joint_rt_irt <- function(...) posterior_means(...)

#' Compare person ability recovery: mean vs quantile RT
compare_ability_mean_vs_quantile <- function(...) compare_theta_mean_vs_quantile(...)

#' Run q-grid comparison for person ability recovery
run_q_grid_ability_compare <- function(...) run_q_grid_theta_compare(...)

#' Plot q-grid comparison for person ability recovery
plot_q_grid_ability <- function(...) plot_q_grid_theta(...)

# =============================================================================
# SECTION 4: 後處理工具
# =============================================================================

#' 合併所有鏈後計算後驗均值
#'
#' @param fit  sample_rtirt() 的輸出，或單條鏈 list
#' @return 各參數後驗均值的 named list
posterior_means <- function(fit) {
  chains <- if (inherits(fit, "rtirt_chains")) fit$chains else list(fit)
  if (length(chains) == 0) stop("沒有可用鏈")
  bad_idx <- which(vapply(chains, function(ch) {
    is.null(ch) || !is.list(ch) || is.null(names(ch))
  }, logical(1L)))
  if (length(bad_idx) > 0) {
    stop("以下鏈結果無效: ", paste(bad_idx, collapse = ", "))
  }
  nm_list <- lapply(chains, names)
  params <- Reduce(intersect, nm_list)
  if (length(params) == 0) stop("各鏈沒有共同參數欄位；請檢查模型輸出一致性")
  lapply(stats::setNames(params, params), function(p) {
    mats <- lapply(chains, function(ch) {
      x <- ch[[p]]
      if (is.null(x))
        stop("參數欄位缺失: ", p, "（某些鏈不存在）")
      if (is.matrix(x)) x else matrix(x, ncol = 1L)
    })
    cm <- colMeans(do.call(rbind, mats))
    if (length(cm) == 1L) as.numeric(cm) else cm
  })
}

#' RMSE
rmse <- function(est, truth) sqrt(mean((est - truth)^2))

#' Bias
bias <- function(est, truth) mean(est - truth)

#' 計算 DIC
#'
#' @param fit              sample_rtirt() 的輸出
#' @param loglik_at_mean   在後驗均值處評估的 log-likelihood（需自行提供）。
#'                         若 NULL，用所有鏈 loglik 最大值代替（低估 pD，僅供參考）。
get_dic <- function(fit, loglik_at_mean = NULL) {
  loglik_all <- if (inherits(fit, "rtirt_chains")) {
    unlist(lapply(fit$chains, `[[`, "loglik"))
  } else {
    fit$loglik
  }

  if (is.null(loglik_at_mean)) {
    warning("未提供 loglik_at_mean，以 max(loglik) 代替，pD 可能被低估")
    loglik_at_mean <- max(loglik_all, na.rm = TRUE)
  }

  compute_dic(loglik_at_mean, loglik_all)
}

#' 模擬回收評估：RMSE 與 Bias
#'
#' @param fit         sample_rtirt() 的輸出
#' @param sim         sim_rtirt_null() 或 sim_rtirt() 的輸出
#' @param model_type  "null" 或 "structural"
#' @return data.frame: param, rmse, bias
evaluate_recovery <- function(fit, sim, model_type = "structural") {
  pm <- posterior_means(fit)
  rows <- list(
    data.frame(param="a",       rmse=rmse(pm$a,       sim$true_a),       bias=bias(pm$a,       sim$true_a)),
    data.frame(param="b",       rmse=rmse(pm$b,       sim$true_b),       bias=bias(pm$b,       sim$true_b)),
    data.frame(param="theta",   rmse=rmse(pm$theta,   sim$true_theta),   bias=bias(pm$theta,   sim$true_theta)),
    data.frame(param="xi",      rmse=rmse(pm$xi,      sim$true_xi),      bias=bias(pm$xi,      sim$true_xi)),
    data.frame(param="tau",     rmse=rmse(pm$tau,     sim$true_tau),     bias=bias(pm$tau,     sim$true_tau)),
    data.frame(param="sigma_t", rmse=rmse(pm$sigma_t, sim$true_sigma_t), bias=bias(pm$sigma_t, sim$true_sigma_t))
  )
  if (model_type == "structural" && !is.null(sim$true_beta)) {
    true_bv <- as.vector(rbind(c(0, 0), sim$true_beta))
    rows <- c(rows, list(data.frame(param="beta", rmse=rmse(pm$beta, true_bv), bias=bias(pm$beta, true_bv))))
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

# =============================================================================
# SECTION 5: 完整使用範例（if(FALSE) 包裹，不會自動執行）
# =============================================================================
if (FALSE) {

# 0. 安裝套件 ------------------------------------------------------------------
# install.packages(c("Rcpp", "RcppArmadillo", "MASS", "coda"))

# 1. 載入 ---------------------------------------------------------------------
Rcpp::sourceCpp("RtIrtGibbs.cpp")
source("RtIrtGibbs_helpers.R")

# ★ 並行時必須設定這個，讓 fork 子程序知道 cpp 在哪裡
.rtirt_cpp_path <- normalizePath("RtIrtGibbs.cpp")

# 2. 模擬資料 ------------------------------------------------------------------
sim <- sim_rtirt_null(n_subj = 300, n_item = 15, true_corr = 0.5, seed = 1)

# 3. 執行多鏈（4 鏈，2000 iter，burnin 1000）-----------------------------------
fit <- sample_rtirt(
  model    = "null",
  Y        = sim$Y,
  logT     = sim$logT,
  n_chains = 4,
  n_iter   = 2000,
  n_burnin = 1000,
  n_cores  = 4        # NULL = 自動偵測
)

# 4. 印出結構 ------------------------------------------------------------------
print(fit)
# rtirt_chains  [model=null | chains=4 | n_save=1000 per chain]
# Stored parameters:
#   $theta       500 columns
#   $tau         500 columns
#   $a           15 columns
#   ...

# 5. 轉換為 coda::mcmc.list ---------------------------------------------------
mc <- to_coda(fit)
#  mc 現在是標準 coda::mcmc.list，每條鏈是一個 coda::mcmc 物件

# 所有 coda 工具都可用：
summary(mc)                                      # 後驗摘要（均值、SD、HPD）

coda::gelman.diag(mc)                            # Gelman-Rubin R-hat，全部參數
coda::gelman.diag(mc[, c("a[1]","b[1]")])        # 只看特定參數

coda::effectiveSize(mc)                          # 有效樣本數
coda::traceplot(mc[, c("a[1]","b[1]","xi[1]"), drop=FALSE])  # trace plot
coda::densplot(mc[, "theta[1]"])                 # 後驗密度
coda::autocorr.plot(mc[, "a[1]"])                # 自相關

# Sigma_p 的 trace plot（person covariance）
mc_sigma <- mc[, c("Sigma_p[1,1]","Sigma_p[2,1]","Sigma_p[2,2]"), drop=FALSE]
coda::traceplot(mc_sigma)

# 6. 快速診斷（只看 item 參數）-------------------------------------------------
diagnose(fit, params = c("a", "b", "xi", "sigma_t", "Sigma_p"), n_show = 5)
# === Gelman-Rubin R-hat（目標 < 1.1） ===
#            Point est. Upper C.I.
# a[1]           1.002      1.007
# ...
# === 有效樣本數（ESS） ===
# a[1]  a[2]  ...
#  923   887  ...

# 7. 模擬回收 -----------------------------------------------------------------
evaluate_recovery(fit, sim, model_type = "null")
#     param      rmse         bias
# 1       a 0.0812...  0.0023...
# 2       b 0.1234...  0.0011...
# ...

# 8. Structural model ---------------------------------------------------------
sim2 <- sim_rtirt(n_subj = 300, n_item = 15, n_feat = 2, seed = 2)
fit2 <- sample_rtirt(
  model    = "structural",
  Y        = sim2$Y,
  logT     = sim2$logT,
  X_cov    = sim2$X,
  n_chains = 4,
  n_iter   = 2000,
  n_burnin = 1000
)
mc2 <- to_coda(fit2)

# beta 的欄位名稱：beta[1,ra], beta[2,ra], beta[1,rt], beta[2,rt]
# （row = covariate index，不含截距；ra/rt = 能力/速度）
coda::traceplot(mc2[, grep("^beta", coda::varnames(mc2)), drop=FALSE])
evaluate_recovery(fit2, sim2, model_type = "structural")

# 9. Quantile model（下四分位數迴歸）------------------------------------------
fit3 <- sample_rtirt(
  model    = "quantile",
  Y        = sim2$Y,
  logT     = sim2$logT,
  X_cov    = sim2$X,
  n_chains = 4,
  n_iter   = 2000,
  n_burnin = 1000,
  q_ra     = 0.25,
  q_rt     = 0.25
)

# 10. DIC 模型比較 -------------------------------------------------------------
dic2 <- get_dic(fit2)   # structural
dic3 <- get_dic(fit3)   # quantile (q=0.25)
cat(sprintf("Structural DIC=%.1f  |  Quantile(q=0.25) DIC=%.1f\n",
            dic2$DIC, dic3$DIC))

# 11. 只轉換部分參數（減少 coda 物件大小）--------------------------------------
mc_item <- to_coda(fit, params = c("a", "b", "xi", "sigma_t"))
# 只有 item 參數，適合快速診斷和繪圖

}  # end if(FALSE)
