# px_sampler.R — 群作用尺度 sampler（PX-DA / generalized Gibbs, Liu & Sabatti 2000）
#
# 在未識別（擴張）的模型上做一維的尺度移動：
#   scaleUp   節點 × alpha
#   scaleDown 節點 / alpha
#   covNode   的第 covIndex 列與行 × alpha（對角元素 × alpha^2）
# alpha = exp(u)，u ~ N(0, scale^2)。正確性來自「對合式 MH」（Tierney 1998；
# Green 1995 的確定性提議）：(x, u) -> (T_alpha x, -u) 是對合，且 q(u) 對稱，故
#   log r = log pi(T_alpha x) - log pi(x) + jacExp * u,
#   jacExp = #scaleUp - #scaleDown + (d + 1)      (d = covNode 維度)
# 其中 (d + 1) 是對稱矩陣第 k 列／行縮放的 Jacobian：Sigma_kk × alpha^2，
# 另外 d - 1 個非對角元素 × alpha。這要求 covNode 的先驗定義在對稱共變矩陣上
# （例如 dinvwish）；若節點是精度矩陣或 Cholesky 因子，數值變換與 Jacobian 都不同，
# 本 sampler 不支援。
#
# 若似然在此變換下不變，可用 control$skipInvariant = TRUE 只計算先驗。
# 這需要模型的參數化在變換下真的不變，例如：
#   a_j * (theta_i - b_j)  → theta, b 放 scaleUp，a 放 scaleDown
#   a_j * theta_i - d_j    → theta 放 scaleUp，a 放 scaleDown，d 不動
#   lambda_j * theta_i     → lambda 也要放 scaleDown
# 預設計算所有相依節點（較慢，但不依賴此假設）。
#
# 注意：scaleUp 節點的先驗變化會與 Jacobian 抵消（例如 theta ~ N(0, Sigma_11) 時，
# 先驗變化 = -n u，Jacobian 貢獻 +n u），所以在 skipInvariant 下接受率只由 Sigma、
# a、b 等的先驗決定，與資料無關。這是 marginal augmentation 的性質：這一步只在
# 未識別的方向上移動，可識別的量（如 a sqrt(Sigma_11)）在軌道上不變，改善只能
# 透過後續的條件更新間接發生（Liu & Wu 1999）。
#
# control:
#   scaleUp, scaleDown : 字元向量（節點或元素，會展開成純量元素；須為非資料的隨機節點）
#   covNode            : 共變矩陣節點名稱（例如 "Sigma[1:2, 1:2]"），可省略
#   covIndex           : 要縮放的列／行（預設 1）
#   skipInvariant      : TRUE 時略過資料節點（需確定似然不變）
#   scale              : log alpha 的初始提議標準差（預設 0.1）
#   adaptive           : 是否調整 scale（預設 TRUE；整條鏈都會調整，步長遞減，
#                        與 nimble 內建 RW sampler 相同）
#   jacExpOverride     : 僅供測試，強制指定 Jacobian 指數

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
    for (nd in c(upNodes, downNodes)) {
      if (!model$isStoch(nd)) stop("pxScale: ", nd, " is not a stochastic node")
      if (model$isData(nd))   stop("pxScale: ", nd, " is a data node")
    }
    if (!hasDown) downNodes <- upNodes[1]     # 佔位，不會被使用
    covExpo <- 0
    covDim  <- 0
    if (hasCov) {
      cv <- model$expandNodeNames(covNode, returnScalarComponents = TRUE)
      covDim <- round(sqrt(length(cv)))
      if (covDim^2 != length(cv)) stop("pxScale: covNode must be a square matrix node")
      if (covIndex < 1 || covIndex > covDim) stop("pxScale: covIndex out of range")
      # 指數 [r == k] + [c == k] 對 (r, c) 對稱，所以不依賴 expandNodeNames 的排列順序
      rc <- expand.grid(r = 1:covDim, c = 1:covDim)
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

    scaleOriginal <- scale
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
      scale <<- scaleOriginal
      timesRan <<- 0; timesAccepted <<- 0; timesAdapted <<- 0
    }
  )
)
