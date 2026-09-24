# px_sampler.R — 群作用尺度 sampler（PX-DA / generalized Gibbs, Liu & Sabatti 2000）
#
# 在未識別（擴張）的模型上做一維的尺度移動：
#   scaleUp   節點 × alpha
#   scaleDown 節點 / alpha
#   covNode   的第 covIndex 列與行 × alpha（對角元素 × alpha^2）
# alpha = exp(u)，u 以對稱 random walk 提議。對乘法群而言 Haar 測度在 u 上即 du，
# 因此接受率為
#   log r = log pi(T_alpha x) - log pi(x) + jacExp * u,
#   jacExp = #scaleUp - #scaleDown + (d + 1)      (d = covNode 維度)
#
# 若似然在此變換下不變（例如 a_j * (theta_i - b_j) 在 theta, b × alpha、a / alpha 下不變），
# 可用 control$skipInvariant = TRUE 只計算先驗；預設計算所有相依節點（較慢但不依賴此假設）。
#
# control:
#   scaleUp, scaleDown : 字元向量（節點或元素，會展開成純量元素）
#   covNode            : 共變矩陣節點名稱（例如 "Sigma[1:2, 1:2]"），可省略
#   covIndex           : 要縮放的列／行（預設 1）
#   skipInvariant      : TRUE 時略過資料節點（需確定似然不變）
#   scale              : log alpha 的初始提議標準差（預設 0.1）
#   adaptive           : 是否在 burn-in 期間調整 scale（預設 TRUE）

sampler_pxScale <- nimbleFunction(
  name = "sampler_pxScale",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    upNodes   <- model$expandNodeNames(control$scaleUp,   returnScalarComponents = TRUE)
    downNodes <- if (is.null(control$scaleDown)) character(0) else
      model$expandNodeNames(control$scaleDown, returnScalarComponents = TRUE)
    covNode   <- if (is.null(control$covNode)) character(0) else control$covNode
    covIndex  <- if (is.null(control$covIndex)) 1 else control$covIndex
    skipInv   <- isTRUE(control$skipInvariant)
    scale     <- if (is.null(control$scale)) 0.1 else control$scale
    adaptive  <- if (is.null(control$adaptive)) TRUE else control$adaptive

    hasDown <- length(downNodes) > 0
    hasCov  <- length(covNode) > 0
    if (!hasDown) downNodes <- upNodes[1]     # 佔位，不會被使用
    covExpo <- 0
    covDim  <- 0
    if (hasCov) {
      cv <- model$expandNodeNames(covNode, returnScalarComponents = TRUE)
      covDim <- round(sqrt(length(cv)))
      rc <- expand.grid(r = 1:covDim, c = 1:covDim)   # column-major，與 values() 一致
      covExpo <- as.numeric((rc$r == covIndex) + (rc$c == covIndex))
    } else {
      cv <- upNodes[1]
      covExpo <- c(0, 0)
    }
    jacExp <- length(upNodes) - (if (hasDown) length(downNodes) else 0) +
      (if (hasCov) covDim + 1 else 0)

    if (!is.null(control$jacExpOverride)) jacExp <- control$jacExpOverride  # 僅供測試
    moved <- c(upNodes, if (hasDown) downNodes, if (hasCov) cv)
    calcNodes <- model$getDependencies(unique(moved))
    if (skipInv) calcNodes <- calcNodes[!model$isData(calcNodes)]

    timesRan <- 0; timesAccepted <- 0; timesAdapted <- 0
    adaptInterval <- 100
  },
  run = function() {
    u <- rnorm(1, 0, scale)
    alpha <- exp(u)
    lp0 <- model$getLogProb(calcNodes)
    values(model, upNodes) <<- values(model, upNodes) * alpha
    if (hasDown) values(model, downNodes) <<- values(model, downNodes) / alpha
    if (hasCov) {
      vc <- values(model, cv)
      for (k in 1:length(vc)) vc[k] <- vc[k] * alpha^covExpo[k]
      values(model, cv) <<- vc
    }
    lp1 <- model$calculate(calcNodes)
    logMHR <- lp1 - lp0 + jacExp * u
    jump <- decide(logMHR)
    if (jump) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
    if (adaptive) {
      timesRan <<- timesRan + 1
      if (jump) timesAccepted <<- timesAccepted + 1
      if (timesRan %% adaptInterval == 0) {
        acc <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        g <- 10 / ((timesAdapted + 3)^0.8)
        scale <<- scale * exp(g * (acc - 0.44))
        timesRan <<- 0; timesAccepted <<- 0
      }
    }
  },
  methods = list(
    reset = function() {
      timesRan <<- 0; timesAccepted <<- 0; timesAdapted <<- 0
    }
  )
)
