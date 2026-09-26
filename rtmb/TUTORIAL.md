# 逐行教學：用 RTMB 寫半模組推論（SMI）RT-IRT

本教學逐行解說 `rtmb/smi_rtmb.R`（56 行）。讀完後你應該能：

1. 寫出這個模型的目標函數，並知道每一項從哪裡來；
2. 解釋「把作答時間的邊際似然乘上 η」為什麼可行，以及為什麼不能改成對 τ 的先驗調溫；
3. 看懂 RTMB 的四個關鍵動作：`getAll`、`random =`、`map =`、`ADREPORT`／`sdreport`；
4. 自己改寫模型，例如換成 logit、加共變數，或加入其他推論對象。

文中用到的數學恆等式都已在 `rtmb/check_identities.R` 用數值驗證過，每一節最後的「動手驗證」可以直接執行。

---

## 0. 先把模型寫清楚

**作答模組（受信任）**，probit 2PL：

$$Y_{ij} = 1\{a_j\theta_i - d_j + e_{ij} > 0\},\quad e_{ij}\sim N(0,1)
\;\Rightarrow\; P(Y_{ij}=1\mid\theta_i) = \Phi(a_j\theta_i - d_j).$$

**作答時間模組（可能設定錯誤）**，lognormal：

$$\log T_{ij} = \xi_j - \tau_i + \gamma_j\theta_i + \varepsilon_{ij},\quad \varepsilon_{ij}\sim N(0,\sigma_j^2),\quad \tau_i \sim N(m_i, v).$$

這裡和 van der Linden (2007) 的寫法有一個差別：原模型讓 (θ, τ) 二元常態相關；這裡讓 τ 與 θ 獨立，能力與時間的全部關聯都由 γ_j 承擔。兩種寫法在可辨識的部分等價，但後者把「作答時間透過什麼管道影響 θ」集中到 γ 一個地方，後面推導洩漏時會很方便。

**conditioning 模型**：要報告的群體變數 X（例如性別）同時進入兩個平均數：

$$\theta_i \sim N(X_i\beta_\theta, 1),\qquad m_i = X_i\beta_\tau.$$

**把 τ 積分掉。** 固定一個人 i，把 p 題的 log 時間排成向量 T_i。給定 τ_i，各題獨立，T_i ~ N(ξ − τ_i 1 + γθ_i, D)，D = diag(σ²)。τ_i 是常態，而且對每題的作用都一樣（係數都是 −1），積分後得到：

$$T_i \mid \theta_i \sim N\big(\xi - m_i\mathbf 1 + \gamma\theta_i,\; \Omega\big),\qquad \Omega = D + v\,\mathbf 1\mathbf 1^\top.$$

直覺是：τ 是一個「整個人整份測驗都快或都慢」的共同成分，積分後就變成各題時間之間的共同共變。

**SMI 的第一階段目標**（取 log 後加負號，就是程式裡的 `nll`）：

$$\text{nll} = -\sum_{ij}\log\Phi\big(S_{ij}(a_j\theta_i - d_j)\big) \;-\; \eta\sum_i \log N(T_i;\,\mu_i,\Omega) \;-\; \sum_i \log N(\theta_i; X_i\beta_\theta, 1) \;-\; \log(\text{先驗}).$$

η 只乘在作答時間的似然上：η = 0 時作答時間完全不影響 θ（cut），η = 1 時就是完整的聯合模型。

> **動手驗證**：`check_identities.R` 第 3 項用數值積分驗證「τ 積分後就是 N(μ, D + v11ᵀ)」，兩邊都是 −4.89219。

---

## 1. 檔頭與函式介面（第 1–14 行）

```r
10 suppressMessages(library(RTMB))
11 smi_rtmb <- function(Y, logT, eta, X = NULL, W, silent = TRUE) {
```

- `Y`：n × p 的 0/1 作答矩陣；`logT`：n × p 的 log 作答時間（秒）。
- `eta`：影響參數，介於 0 到 1。
- `X`：conditioning 共變數（n × q），要先置中，而且不含截距。θ 的平均數由先驗的 0 固定，τ 的平均數由 ξ_j 吸收，放截距會造成不可辨識。
- `W`：推論對象的權重矩陣（n × k）。每一欄定義一個推論對象 Wᵀθ，例如群體差距或迴歸斜率，見第 54–56 行。

```r
12   n <- nrow(Y); p <- ncol(Y)
13   X <- if (is.null(X)) matrix(0, n, 1) else as.matrix(X)
14   q <- ncol(X)
```

沒有共變數時放一欄 0，讓後面的 `X %*% b_th` 不必分兩種情況處理。這欄對應的 β 只受先驗約束，不影響任何結果。

---

## 2. 資料前處理與參數初始值（第 15–18 行）

```r
15   S <- 2 * Y - 1
```

**probit 的對稱技巧。** P(Y = 1) = Φ(η)，P(Y = 0) = 1 − Φ(η) = Φ(−η)。令 S = 2Y − 1 ∈ {−1, +1}，兩種情況就合成一個式子 P(Y) = Φ(S·η)，後面只需要寫一次。

```r
16   par <- list(theta = rep(0, n), a = rep(1, p), d = rep(0, p), b_th = rep(0, q),
17               xi = colMeans(logT), gam = rep(0, p), ls2 = rep(log(.25), p), lv = log(.3),
18               b_tau = rep(0, q))
```

RTMB 的參數必須放在一個具名 list 裡，名字就是之後 `getAll` 解開時的變數名。

| 參數 | 意義 | 尺度 | 初始值的理由 |
|---|---|---|---|
| `theta` | 能力（隨機效果） | 原尺度 | 先驗平均 0 |
| `a`, `d` | 鑑別度、難度 | 原尺度 | 中性值 |
| `b_th` | θ 對 X 的迴歸係數 | 原尺度 | 0 |
| `xi` | 題目時間強度 | log 秒 | 各題 log 時間平均，最接近解 |
| `gam` | 能力對時間的效果 | 原尺度 | 0，從「無關聯」出發 |
| `ls2` | log σ_j² | **log 尺度** | 殘差 SD 約 0.5 |
| `lv` | log v（速度變異） | **log 尺度** | 與模擬設定同量級 |
| `b_tau` | τ 對 X 的迴歸係數 | 原尺度 | 0 |

**為什麼變異數要放在 log 尺度？** 最佳化器（`nlminb`）是無約束的，直接估 σ² 可能跑到負數；估 log σ² 就沒有這個問題，Laplace 近似在 log 尺度上通常也更接近常態。

---

## 3. 目標函數：作答模組（第 19–25 行）

```r
19   f <- function(par) {
20     getAll(par)
```

`getAll(par)` 把 list 裡的每個元素解開成區域變數（`theta`、`a`、`d`……），讓後面的程式可以照數學式的樣子寫。RTMB 會在這個函數上錄製自動微分（AD）的運算帶（tape），所以函數內只能用 RTMB 支援 AD 的運算：矩陣運算、`pnorm`、`dnorm`、`log`、`exp` 都可以；`if` 只能依賴資料或常數（例如 `eta`），不能依賴參數值。

```r
21     lin <- matrix(theta, n, p) * matrix(a, n, p, byrow = TRUE) - matrix(d, n, p, byrow = TRUE)
```

建立 n × p 的線性預測子 a_jθ_i − d_j。
- `matrix(theta, n, p)`：把 θ 複製成 p 欄，每一**列**是同一個人。
- `matrix(a, n, p, byrow = TRUE)`：把 a 複製成 n 列，每一**欄**是同一題。
- 兩者逐元素相乘，就得到 a_jθ_i。

這種「展開成矩陣再逐元素運算」的寫法 AD 支援得很好。如果用 `outer()`，在某些 RTMB 版本不支援 AD 型別，所以避開。

```r
22     nll <- -sum(log(pnorm(S * lin)))
```

作答模組的負對數似然：−Σ log Φ(S_ij · lin_ij)。

注意：當 S·lin 很負時（例如 −10），`pnorm` 會下溢為 0，`log` 變成 −∞。在正常的 2PL 範圍（|lin| < 6）沒有問題；如果資料有極端題目，可以改成 `pnorm(S * lin, log.p = TRUE)`（需確認你的 RTMB 版本支援），或對 lin 截斷。

```r
23     nll <- nll - sum(dnorm(theta, drop(X %*% b_th), 1, log = TRUE))
```

θ 的先驗：N(X_iβ_θ, 1)，也就是 **conditioning 模型**。它有兩個作用：
1. θ 的尺度由變異數固定為 1 來識別；
2. 群體差距不會被收縮。沒有這一項時，θ 的後驗會把兩組都往 0 拉，真實的群體差距被低估，大約是乘上信度（模擬中 cut 的偏誤是 −0.055，約等於 1 − 0.82）。

`drop()` 把 n × 1 矩陣變成向量，避免 `dnorm` 的長度循環出錯。

```r
24     nll <- nll - sum(dnorm(a, 1, 1, log = TRUE)) - sum(dnorm(d, 0, 2, log = TRUE)) -
25       sum(dnorm(b_th, 0, sqrt(10), log = TRUE))
```

作答模組參數的先驗：a ~ N(1, 1)、d ~ N(0, 2²)、β ~ N(0, 10)，與 Gibbs 版本完全相同。這裡先驗當作懲罰項，所以 RTMB 求的是**後驗眾數**（MAP），不是最大概似估計。這是刻意的：要和 Gibbs 的後驗比較，先驗必須一致。

---

## 4. 目標函數：作答時間模組（第 26–38 行）

```r
26     if (eta > 0) {
```

η = 0 時整段略過，作答時間對任何參數都沒有貢獻，就是 cut。這個 `if` 依賴的是常數 `eta`，不是參數，所以不影響 AD。

```r
27       w <- exp(-ls2); v <- exp(lv); Sw <- sum(w)
```

- w_j = 1/σ_j²：每題的**精確度**。
- v：速度變異。
- Sw = Σ_j w_j。

```r
28       m <- drop(X %*% b_tau)
```

τ 的平均數 m_i = X_iβ_τ（conditioning 模型的作答時間那一半）。例如 G = 1 那組整體比較慢，就會反映在 β_τ < 0。

```r
29       R <- logT - matrix(xi, n, p, byrow = TRUE) + matrix(m, n, p) -
30         matrix(theta, n, p) * matrix(gam, n, p, byrow = TRUE)
```

殘差 r_ij = log T_ij − ξ_j + m_i − γ_jθ_i，也就是 log T 減去它在 τ 積分後的平均數 ξ_j − m_i + γ_jθ_i。注意 m 前面是「+」號：模型裡 τ 以「−τ」進入 log T，所以平均數是 −m，減掉 −m 就變成 +m。

```r
31       rw <- drop(R %*% w)
32       quad <- drop((R * R) %*% w) - v * rw^2 / (1 + v * Sw)
```

這兩行計算每個人的二次式 r_iᵀΩ⁻¹r_i，卻**完全不做矩陣求逆**。

**Sherman–Morrison 公式。** Ω = D + v11ᵀ 是「對角矩陣加一個秩一矩陣」，它的逆矩陣有封閉解：

$$\Omega^{-1} = D^{-1} - \frac{v\,D^{-1}\mathbf 1\mathbf 1^\top D^{-1}}{1 + v\,\mathbf 1^\top D^{-1}\mathbf 1}.$$

代入二次式：

$$r^\top\Omega^{-1}r = \sum_j w_j r_j^2 \;-\; \frac{v\,(\sum_j w_j r_j)^2}{1 + v\sum_j w_j}.$$

- 第一項 `(R * R) %*% w` = Σ w_j r_j²：把各題殘差平方加權加總，就是「如果 τ 不存在」的二次式。
- 第二項扣掉 v(Σ w_j r_j)²/(1 + vΣw)：扣掉「全部題目一起慢／快」的那部分。這正是 τ 解釋掉的共同成分。

這樣每個人只要 O(p) 的運算，而不是 O(p³)。n = 5000、p = 30 時差別很大，而且 AD 通過矩陣求逆很貴，封閉式對 AD 也更友善。

```r
33       ll <- -0.5 * (p * log(2 * pi) + sum(ls2) + log(1 + v * Sw) + quad)
```

多元常態的對數密度：−½[p log 2π + log|Ω| + rᵀΩ⁻¹r]。

**矩陣行列式引理**：|D + v11ᵀ| = |D|(1 + v1ᵀD⁻¹1)，所以

$$\log|\Omega| = \sum_j \log\sigma_j^2 + \log\Big(1 + v\sum_j w_j\Big) = \texttt{sum(ls2)} + \texttt{log(1 + v * Sw)}.$$

`ll` 是長度 n 的向量，每個元素是一個人的 log L_i。

> **動手驗證**：`check_identities.R` 第 1、2 項用 `solve()` 和 `determinant()` 直接算，與封閉式完全相同（14.4257、−6.223988）。

```r
34       nll <- nll - eta * sum(ll)
```

**整個方法的核心就是這一行**：作答時間的邊際似然乘上 η。

為什麼要乘在**邊際**似然上，而不是分開乘在 p(T | τ) 和 p(τ) 上？
- 如果把整個模組連同 τ 的先驗一起調溫，[p(T | τ)p(τ | v)]^η，τ 的先驗密度被開 η 次方後，對 τ 積分會多出一個 v^{(1−η)/2} 的因子，隱含的 v 分布不可積（不當分布）。`check_identities.R` 第 5 項顯示 ∫N(τ; 0, v)^{0.5}dτ 隨 v 增長（2.24、3.17、4.48，約正比於 v^{1/4}）。
- 只對 p(T | τ) 調溫、先驗不動（Gibbs 版的 `temper = "likelihood"`）雖然合法，但 η 和「借用多少資訊」的對應很不均勻：`smi/README.md` 的 `exp_temper_compare.R` 顯示，η = 0.25 時就已經傳過約 71% 的資訊增益（邊際調溫只有 19%），η 在 0.25 以上幾乎沒有調節作用。兩種調溫在相同資訊量下的偏誤差不多，所以差別主要在 η 好不好用，不在偏誤與精確度的取捨本身。
- 對**邊際**似然調溫，意思清楚：「作答時間這個資料來源，對 θ 的整體證據只算 η 份」。借用的資訊大致隨 η 等比例上升（η = 0.5 時約 38%，0.75 時約 63%），η 因此成為平滑的旋鈕。

Gibbs 版本要用 τ̃ 擴增加上 |Ω| 校正項才能做到同樣的事；RTMB 可以直接寫出邊際似然，這一步反而更簡單。

```r
35       # priors: xi ~ N(4, 10^2), gamma ~ N(0, 1), s^2 ~ IG(1, 1), v ~ IG(1, 1) (on the log scale)
36       nll <- nll - sum(dnorm(xi, 4, 10, log = TRUE)) - sum(dnorm(gam, 0, 1, log = TRUE)) -
37         sum(-ls2 - exp(-ls2)) - (-lv - exp(-lv)) - sum(dnorm(b_tau, 0, sqrt(10), log = TRUE))
```

作答時間模組的先驗，與 Gibbs 版相同。重點在變異數的先驗：σ² ~ IG(1, 1) 的密度是 (σ²)^{−2}e^{−1/σ²}。因為我們估的是 u = log σ²，要乘上 Jacobian |dσ²/du| = σ²：

$$p(u) = (\sigma^2)^{-2}e^{-1/\sigma^2}\cdot\sigma^2 = e^{-u}\,e^{-e^{-u}} \;\Rightarrow\; \log p(u) = -u - e^{-u}.$$

這就是 `-ls2 - exp(-ls2)`。少了 Jacobian，眾數會落在不同的位置，和 Gibbs 的後驗就對不上。

> **動手驗證**：`check_identities.R` 第 4 項確認 exp(−u − e^{−u}) 在整條實數線上積分為 1。

**先驗在 `if` 裡面的原因**：η = 0 時作答時間參數被 map 掉（見第 43 行），不再是自由參數，它們的先驗就不該留在目標函數裡。

---

## 5. 推論對象與回傳（第 39–42 行）

```r
39     tgt <- drop(t(W) %*% theta)
40     ADREPORT(tgt)
41     nll
42   }
```

- `tgt`：k 個推論對象，每個都是 θ 的線性組合 W_kᵀθ。
- `ADREPORT(tgt)`：告訴 RTMB 這些量要在最後報告估計值與標準誤。因為 tgt 依賴隨機效果 θ，`sdreport` 會用**廣義 delta 法**同時考慮兩種不確定性：固定效果估計的不確定性，以及給定固定效果時 θ 的條件不確定性（由 Laplace 的 Hessian 給出）。這讓它的標準誤可以和 Gibbs 的後驗 SD 直接比較；實測中位數比值為 1.02。
- 函數回傳純量 `nll`，RTMB 對它做 Laplace 積分與最佳化。

---

## 6. `map`：怎麼做出 cut（第 43–45 行）

```r
43   map <- if (eta > 0) list() else
44     list(xi = factor(rep(NA, p)), gam = factor(rep(NA, p)), ls2 = factor(rep(NA, p)),
45          lv = factor(NA), b_tau = factor(rep(NA, q)))
```

`map` 是 TMB 固定參數的機制：把參數對應到 `factor(NA)` 就表示「固定在初始值，不估計」。
- η > 0：空的 list，全部參數都估計。
- η = 0：作答時間模組的參數全部固定。它們在 `nll` 裡也完全不出現（第 26 行的 `if` 已略過），所以固定在什麼值都無所謂。這樣 θ 只由作答資料決定，就是 **cut posterior 的第一階段**。

為什麼不直接讓它們自由？η = 0 時這些參數只受先驗約束，沒有資料資訊，最佳化器會在平坦的方向上遊走，Hessian 也可能接近奇異，導致 `sdreport` 出錯。

---

## 7. 建模、最佳化、報告（第 46–51 行）

```r
46   obj <- MakeADFun(f, par, random = "theta", map = map, silent = silent)
```

RTMB 的核心呼叫：
- 錄製 `f` 的 AD 運算帶；
- `random = "theta"` 指定 θ 為**隨機效果**。TMB 會自動做 Laplace 近似：

$$\int e^{-\text{nll}(\phi,\theta)}\,d\theta \;\approx\; e^{-\text{nll}(\phi,\hat\theta)}\,(2\pi)^{n/2}\,\big|H(\hat\theta)\big|^{-1/2},$$

其中 θ̂ 是給定固定效果時的內層最佳解，H 是對 θ 的 Hessian。θ_i 彼此條件獨立，所以 H 是**對角**的；TMB 會自動偵測這種稀疏結構，n 個人的 Laplace 其實是 n 個一維 Laplace，成本線性於 n。

**Laplace 在這裡有多準？** 之前用 AGHQ 比較過：一維 θ 的 probit 模型，Laplace 對每人的對數邊際似然大約差 0.007。作答時間部分在 θ 上是二次的（常態），Laplace 對這部分是精確的；誤差只來自 probit 項。

```r
47   opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 500, eval.max = 1000))
```

在固定效果上最小化 Laplace 近似後的負對數邊際後驗。`obj$gr` 是 AD 算出的精確梯度，而且已經考慮了 θ̂ 會隨固定效果移動（隱函數定理）。一定要檢查 `opt$convergence == 0`；函式把它回傳在 `conv` 裡。

```r
48   sdr <- sdreport(obj)
49   rep <- summary(sdr, "report")
```

`sdreport` 用固定效果的 Hessian 求共變異數，再以 delta 法傳到每個 `ADREPORT` 的量。`summary(sdr, "report")` 取出推論對象的估計與標準誤。

```r
50   list(est = rep[, 1], se = rep[, 2], theta = summary(sdr, "random")[, 1], opt = opt,
51        conv = opt$convergence, eta = eta)
```

回傳推論對象的估計（`est`）與標準誤（`se`），以及每個人的 θ̂（條件眾數，即經驗貝氏分數）。

---

## 8. 推論對象的權重（第 54–56 行）

```r
55 w_gap <- function(b) ifelse(b == 1, 1 / sum(b == 1), -1 / sum(b == 0))
```

群體差距 = 組 1 平均 − 組 0 平均 = Σ_i w_iθ_i，其中組 1 的人 w_i = 1/n_1，組 0 的人 w_i = −1/n_0。

```r
56 w_slope <- function(z) { zc <- z - mean(z); zc / sum(zc^2) }
```

θ 對連續變數 z 的 OLS 斜率 = Σ(z_i − z̄)θ_i / Σ(z_i − z̄)²，也是 θ 的線性組合。

**關鍵觀念：推論對象要是 θ 的線性組合。** 只要能寫成 Wᵀθ（差距、斜率、組內平均、對比），`ADREPORT` 都能直接處理。非線性的推論對象（例如組間 SD 比）也可以放進 `ADREPORT`，但 delta 法的近似會比較粗。

---

## 9. 把它用起來：η 的選擇

```r
source("rtmb/smi_rtmb.R")
X <- cbind(G = dat$G - mean(dat$G))                  # conditioning 變數，置中
W <- cbind(G = w_gap(dat$G), Z = w_gap(dat$Z))       # 兩個推論對象
etas <- c(0, .1, .25, .5, .75, 1)
fits <- lapply(etas, function(e) smi_rtmb(dat$Y, dat$logT, e, X, W))
m <- sapply(fits, function(f) f$est["Z"]); V <- sapply(fits, function(f) f$se["Z"]^2)
D <- m - m[1]                                        # 與 cut 的差
R <- D^2 - V[1] + 2 * V                              # 風險估計
etas[which.min(R)]                                   # 對 Z 差距選的 η
```

**風險公式從哪來？** 以 cut 的估計 Δ̂_0 當作不偏基準。對某個 η：

$$\text{MSE}(\eta) = \text{Var}(\hat\Delta_\eta) + \text{Bias}_\eta^2,\qquad \text{Bias}_\eta = E[D],\quad D = \hat\Delta_\eta - \hat\Delta_0.$$

E[D²] = Bias² + Var(D)，所以 Bias² 的不偏估計是 D² − Var(D)。在模型正確時 Δ̂_η 比 Δ̂_0 有效率，Hausman 的結果給出 Var(D) ≈ V_0 − V_η。代入：

$$\widehat{\text{MSE}}(\eta) = V_\eta + D^2 - (V_0 - V_\eta) = D^2 - V_0 + 2V_\eta.$$

- D 大（作答時間把差距推走）→ 風險高 → 選小 η。
- D 小 → 風險主要由 V_η 決定 → 選大 η，拿精確度。

RTMB 的 V 是確定性的（沒有 MCMC 雜訊），所以 D 和 V_0 − V_η 都比 Gibbs 版更穩定。

---

## 10. 練習

1. **換成 logit。** 把第 22 行改成 `-sum(dbinom(Y, 1, plogis(lin), log = TRUE))`，比較 θ 的尺度變化（約 1.7 倍）。
2. **確認 cut 真的與作答時間無關。** η = 0 時，把 `logT` 整個打亂重排，`est` 應該完全不變。
3. **看 η 的平滑性。** 在 η ∈ seq(0, 1, .05) 上畫 Z 差距的估計與 θ 的平均標準誤；邊際調溫下兩者都應該隨 η 平滑地變化。
4. **驗證 Laplace 的準確度。** 對單一個人，用 `integrate()` 對 θ 做數值積分，與 Laplace 的結果比較。
5. **加入第二個 conditioning 變數。** `X <- cbind(G, Z)`（兩欄都置中），確認 Z 差距的洩漏消失（對應模擬的 full+GZ）。
6. **洩漏係數。** 用擬合結果的 γ、σ²、v 和 θ 的精確度計算 λ(η) = η·1ᵀΩ⁻¹γ / (I + η·γᵀΩ⁻¹γ)，對照不同 η 下 Z 差距的移動量。

## 11. 和 Gibbs 版的對照

| | Gibbs（`smi/smi_gibbs.R`） | RTMB（`rtmb/smi_rtmb.R`） |
|---|---|---|
| τ | 用 τ̃ 擴增抽樣，需要 \|Ω\| 校正項 | 解析積分，直接寫邊際似然 |
| θ | 共軛常態抽樣 | Laplace 積分 |
| 輸出 | 後驗抽樣（完整分布） | 眾數 + delta 法標準誤 |
| 推論對象的變異 | 有 MC 雜訊 | 確定性 |
| 第二階段（RT 參數） | 巢狀 MCMC，C++ 加速 | 要另外擬合 |
| 速度（n = 500、p = 10） | 約 10 秒 | 約 2 秒 |
| 實測一致性 | — | 估計差 ≤ 0.011，SE 比值 1.02 |
