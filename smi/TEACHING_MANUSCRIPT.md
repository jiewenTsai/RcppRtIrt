# 作答時間該借多少？半模組推論與 RT-IRT 的群體比較

> 這份稿子的目的有兩個：**教會你**這個方法在做什麼，以及**說服你**它值得寫成一篇論文。
> 每個數字都來自本 repo 的模擬（出處見 `smi/README.md`）。文獻標記「✓」者已在本次研究中用 Consensus 查證；標記「†」者是憑記憶引用，投稿前需要核對。

---

## 0. 一句話

聯合 RT-IRT 用作答時間提高能力估計的精確度，代價是：只要群體之間有與能力無關的速度差，這個差異就會「洩漏」進能力，推移群體比較。精確度指標看不出來，階層模型的 conditioning 也擋不住題目層次的速度差。**本研究把「用不用作答時間」變成「對這個推論，借多少作答時間」**，並給出一個可以從資料決定借用量的規則。

---

## 1. 心理計量的問題：附帶資訊是雙面刃

### 1.1 附帶資訊的承諾

電腦化測驗讓我們免費拿到每一題的作答時間。van der Linden (2007)† 的階層模型把作答（IRT）和作答時間（lognormal）放在同一個模型裡，第二層讓能力 θ 和速度 τ 相關，於是作答時間成為估計能力的**附帶資訊**（collateral information; van der Linden et al., 2010 ✓）。這個想法很吸引人：題目數量不變，能力估計卻變準了。後續研究把它用在：

- 適性測驗：用作答時間輔助 EAP 估計，或選擇「每單位時間資訊量最大」的題目（Kern et al., 2021 ✓；He et al., 2023 ✓）；
- 更精緻的聯合模型：允許殘差時間與作答之間的交叉負荷，以榨出更多資訊（Bolsinova & Tijmstra, 2018 ✓）。

這些研究評估成效的方式幾乎都一樣：**看能力估計的 RMSE、信度、或資訊量有沒有改善**。在我們的模擬中，改善確實存在：θ 的 RMSE 從只用作答資料的 0.44 降到聯合模型的 0.37。

### 1.2 附帶資訊的代價：借用管道不分青紅皂白

作答時間之所以能幫忙估能力，是因為模型假設「大家共用同一套速度分布和同一組題目時間參數」。在這個假設下，一個人比預期慢，模型就會把它部分解讀為能力訊號（在本研究的參數化中，透過 γ_j 這條管道）。

問題是，模型不問**為什麼慢**。如果某個群體因為與能力無關的原因作答較慢：

- 多語學習者讀題較慢（Park et al., 2024 ✓ 在 PISA 2018 發現多語學生整體較慢，而且有 7–10% 的題目出現差異作答時間）；
- 性別之間的作答節奏不同（Kapoor et al., 2024 ✓ 把 PISA 的性別差距拆成速度與能力兩部分）；
- 延長時間的考場調整、裝置差異、答題策略……

那麼這些差異會沿著同一條管道進入能力估計。我們把它稱為**速度洩漏（speed leakage）**。

### 1.3 洩漏有多大、為什麼看不見

模擬結果（n = 500，p = 10，性別真實差距 0.3 SD）：

| 情境 | 完整聯合模型的差距偏誤 | 95% 區間覆蓋率 |
|---|---|---|
| 模型外變數 Z 的那組整體慢 0.25 | +0.07 SD | 0.42 |
| 模型外變數 Z 的那組整體慢 0.5 | +0.13 SD | 0.06 |
| Z 那組只在 3 題上慢 1，且這些題 γ 較高 | +0.44 SD | 0 |
| 同上情境，但看**已放進模型的性別**差距 | +0.21 SD | 0.06 |

而在所有這些情境裡，**θ 的 RMSE 都仍然偏好完整模型**。原因很簡單：RMSE 是對所有人平均，洩漏只是讓某個群體整體平移一點點，對個人層次的總誤差影響很小；但群體比較恰恰只看這個平移。

**這就是 P（問題）：精確度的評估標準和群體比較的效度要求是兩回事，而現行文獻只檢查前者。**

---

## 2. 文獻地圖：四條線交會的地方

這篇研究站在四條文獻的交會處，理解每一條，就知道貢獻在哪裡。

### 2.1 作答時間當附帶資訊（本研究的對話對象）

- van der Linden (2007)†；van der Linden, Klein Entink & Fox (2010, APM) ✓：階層框架與附帶資訊的兩個來源。**我們的立場：部分同意**。增益是真的，但前提是速度模型對所有人成立；群體速度差和 DRT 會把增益變成群體偏誤。
- Bolsinova & Tijmstra (2018, BJMSP) ✓：用交叉負荷取得更多附帶資訊。**部分同意／反對**：借得越多，洩漏管道越寬；而以 RMSE 評估正好看不到洩漏。
- Kern et al. (2021, APM) ✓：在 CAT 用 J-EAP「控制差異速度」。**部分同意**：它評估的是個人準確度，一旦拿來做群體比較，就需要本研究的檢查。

### 2.2 群體速度差與差異作答時間（證明問題真實存在）

- Kapoor et al. (2024, JEM) ✓：PISA 各國的性別差距中有一部分來自速度差異。
- Park et al. (2024, EMIP) ✓：多語學習者作答較慢，7–10% 的題目有 DRT，而且 DIF 題與 DRT 題重疊很少。
- Duan et al. (2024, EPM) ✓、Molenaar et al. (2024, Psych Methods) ✓：作答時間 DIF 的偵測與解釋。

這條線的角色是證明 P 的前提（群體速度差）在真實資料中普遍存在，不是模擬者自己編出來的。

### 2.3 Plausible values 的 conditioning 模型（最接近的類比，也是最重要的對照）

大型評量用 conditioning 模型產生 plausible values：把背景變數放進 θ 的先驗，次級分析才不會有偏誤（Wu, 2005 ✓；Mislevy, 1991†）。已知的問題是：**沒放進 conditioning 模型的變數，分析結果會有偏誤**（Monseur & Adams, 2009 ✓；Bailey et al., 2023 ✓ 的 Dire 套件甚至點名 process data）。

本研究與它的關係很微妙，也是說服審稿人的關鍵：

| | 遺漏的 conditioning 變數（PV 文獻） | 速度洩漏（本研究） |
|---|---|---|
| 偏誤方向 | 差距往 0 縮小（衰減） | 加法性，方向由 γ 與速度差決定，可以製造或反轉差距 |
| 隨測驗加長 | 消失（Marsman et al.†） | 也會變小，但與題目時間參數有關 |
| 修正方式 | 重新 conditioning（Dire） | 人層次的速度差可以靠 conditioning τ 修正；**題目層次的 DRT 不行** |

所以我們的一個核心發現是：**階層模型的標準解法（把變數放進 conditioning）在作答時間上不夠**。

### 2.4 模組化貝氏推論：cut 與 SMI（方法的來源）

- **cut posterior**：模型由幾個模組組成，其中一個可能設定錯誤時，切斷它對其他模組的回饋（Liu et al., 2025, JRSS-B ✓ 給出一般定義；Plummer, 2015† 指出 WinBUGS 的 `cut()` 抽的不是正確的 cut 分布）。
- **Levy (2024, JEBS; 2026, BJMSP)** ✓：心理計量中的 measurement-preserving 多階段貝氏。它的出發點是 *interpretational confounding*：加入外部變數後，測量模型被改寫。**我們支持並延伸**：作答時間模組正是一種會改寫 θ 意義的外部資訊。
- **Semi-modular inference**（Carmona & Nicholls, 2020 ✓）：用影響參數 η 在 cut（η = 0）與完整模型（η = 1）之間連續內插。
- **Frazier et al. (2025, JASA)** ✓：證明在後驗風險意義下 cut 並非最佳，並提出利用偏誤與變異取捨的 SMI。**我們支持並延伸**：給出一個針對潛在變數模組、針對特定推論對象的可操作規則。

模組化貝氏在統計學界已有十年的發展，但在心理計量中，除了 Levy 的工作，幾乎沒有被用來處理「輔助資料來源會污染測量」這個問題。作答時間是最自然的第一個應用。

---

## 3. 核心概念教學

### 3.1 兩個模組

$$\underbrace{p(Y \mid \theta, a, d)}_{\text{作答模組（信任）}}\qquad\underbrace{p(T \mid \theta, \xi, \gamma, \sigma, v)}_{\text{作答時間模組（可疑）}}$$

θ 是兩個模組共用的參數，這就是資訊（和污染）流動的地方。

### 3.2 三種推論

1. **完整貝氏（η = 1）**：兩個模組一起決定 θ。精確，但作答時間模組的任何錯誤都會流進 θ。
2. **cut（η = 0）**：θ 只由作答決定；作答時間模組之後再以給定的 θ 估計。θ 不受污染，但放棄了全部的精確度增益。
3. **SMI（0 < η < 1）**：作答時間對 θ 的證據只算 η 份：

$$p_\eta(\theta, \ldots) \propto p(Y\mid\theta,\ldots)\; L_{RT}(\theta,\ldots)^{\eta}\; p(\theta)\,p(\ldots).$$

### 3.3 為什麼要對「邊際」似然調溫

L_RT 是把速度 τ 積分掉後的作答時間似然，一個人的 p 題時間服從 N(ξ − m + γθ, Ω)，Ω = diag(σ²) + v11ᵀ。

- 若對含 τ 先驗的整個模組調溫，τ 的先驗被開 η 次方後會產生不可積的 v 分布（不當後驗）。
- 若只對 p(T | τ) 調溫、先驗不動，η 與借用資訊的對應很不均勻（η = 0.25 就借走約 71% 的資訊）。
- 對邊際似然調溫，借用量大致隨 η 等比例上升（η = 0.5 約 38%，0.75 約 63%），η 成為好用的旋鈕。

### 3.4 洩漏的封閉式

在常態近似下，某個群體若有未被模型吸收的速度差 δ，其能力差距會被推移約 λ(η)·δ：

$$\lambda(\eta) = \frac{\eta\,\mathbf 1^\top\Omega^{-1}\gamma}{I_{\text{IRT}} + \eta\,\gamma^\top\Omega^{-1}\gamma}.$$

分子是「速度差透過 γ 流進 θ 的強度」，分母是「θ 的總精確度」。作答資訊越多（I_IRT 越大），洩漏越小；η 越小，洩漏越小。在模擬中，這個公式預測的偏誤與實際偏誤接近，並略為高估（例如 η = 1：預測 0.14，實際 0.12），是保守的界限。

---

## 4. 本研究的做法（細節）

### 4.1 模型

- 作答：probit 2PL，Φ(a_jθ_i − d_j)。
- 作答時間：log T_ij = ξ_j − τ_i + γ_jθ_i + ε_ij，τ ⟂ θ。能力與時間的關聯全部由 γ_j 承擔，消除了剪切方向的不可辨識性。
- conditioning：要報告的群體變數 X 同時進入 θ 與 τ 的平均數。

### 4.2 估計（兩種實作，結果一致）

- **Gibbs**（`smi/smi_gibbs.R`）：Albert–Chib 擴增 + z-scale PX-DA（鑑別度 ESS 中位數 +56%）。邊際調溫用擴增變數 τ̃ ~ N(m, v/η) 加校正項 |Ω|^{n(1−η)/2}，讓所有主要條件分布仍是共軛。每次擬合約 10 秒。
- **RTMB**（`rtmb/smi_rtmb.R`）：直接寫出邊際似然的 η 次方，θ 用 Laplace 積分。每次約 2 秒。與 Gibbs 的估計差 ≤ 0.011，標準誤比值中位數 1.02。

### 4.3 η 的選擇：針對推論對象的風險規則

對每個推論對象 Δ（群體差距、迴歸斜率），在 η 網格上計算：

$$D(\eta) = \hat\Delta_\eta - \hat\Delta_0,\qquad \hat R(\eta) = D(\eta)^2 - V_0 + 2V_\eta,$$

選 R̂ 最小的 η。直覺：D 大代表作答時間正在把估計推離不偏基準，此時應少借；D 小則借越多越好，因為 V_η 會下降。

**為什麼基準是 conditioning 過的 cut？** 沒有 conditioning 時，cut 本身也會把真實差距往 0 縮小（偏誤 −0.055），基準本身就不對。

### 4.4 一次擬合的篩檢

把每人的速度位置分數 E[(τ_i − m_i)/v | data] 對 conditioning 變數投影後，依候選變數排序累加。模型正確時，這個累加過程近似 Brownian bridge，可以做 double-max 或 LM 檢定，用來決定哪些變數需要跑整組 η。

### 4.5 事後校正（負面結果，但有教育意義）

我們曾嘗試讓產分數的單位釋出速度殘差 u 與洩漏比 c，讓次級分析者用 gap(θ) − c·gap(u) 自行修正。在人層次的速度差下效果很好（偏誤 0.12 → 0.00），但在 DRT 下完全失效（0.45 → 0.46），因為題目層次的速度差不會出現在人層次的 u 裡。**這個失敗正是「為什麼需要以 cut 為基準的規則」的最好證明**：任何只在人層次處理速度的方法，都擋不住 DRT。

---

## 5. 證據

### 5.1 風險規則（每種情境 50–100 份資料）

| 比較對象 | 情境 | 完整模型：偏誤／覆蓋率 | 風險規則：偏誤／覆蓋率 | 選出的 η |
|---|---|---|---|---|
| Z | 無速度差 | 0.000 / 0.95 | 0.000 / 0.91 | 0.79 |
| Z | 整體慢 0.25 | 0.073 / 0.42 | 0.006 / 0.90 | 0.13 |
| Z | 整體慢 0.5 | 0.125 / 0.06 | 0.004 / 0.88 | 0.06 |
| Z（連續） | 線性速度漂移 | 0.063 / 0.08 | 0.001 / 0.92 | 0.05 |
| Z（連續） | 門檻型漂移 | 0.043 / 0.26 | 0.002 / 0.92 | 0.11 |
| Z | 3 題 DRT、高 γ | 0.439 / 0 | 0.002 / 0.96 | 0.00 |
| 性別（已在模型中） | 3 題 DRT、高 γ | 0.205 / 0.06 | 0.013 / 0.96 | 0.06 |

### 5.2 精確度保留

沒有速度差時，風險規則平均選 η ≈ 0.6–0.8，θ 的 RMSE 為 0.38–0.39（cut 0.44、完整模型 0.37）。所以規則並不是「一律退回 cut」，而是只在有證據時才少借。

### 5.3 一個關鍵觀察

群體差距的後驗 SD：完整模型 0.042，cut 0.044。**作答時間幾乎沒有提高群體比較的精確度**，它的價值在個人分數。這一點讓實務建議變得很簡單：做群體比較時，少借的代價很小。

### 5.4 篩檢

對線性與門檻型的速度漂移，檢定力都是 1.00；拿無關變數排序的型一錯誤率為 0.00–0.06。

---

## 6. 為什麼你應該相信：逐一回應最強的反對意見

**反對一：「把速度相關的變數放進階層模型就好，何必 SMI？」**
對已知、且作用在人層次的變數，確實可以，而且我們的模擬也證實 full+G 對 G 沒有偏誤。但有兩個情況不行：(1) 估分時還不存在的變數（次級分析）；(2) 題目層次的速度差。在 DRT 情境下，連已經放進模型的性別差距都偏了 0.21 SD。

**反對二：「那乾脆都用 cut，只用作答資料就好。」**
對群體比較，這其實是個不錯的預設，我們也這樣建議釋出資料的單位。但 cut 放棄了個人分數的精確度（RMSE 0.44 對 0.37），而且 Frazier et al. (2025) 已證明 cut 在後驗風險下並非最佳。風險規則在沒有問題時保住大部分增益，有問題時才退回 cut，兩邊都照顧到。

**反對三：「為什麼不直接把 DRT 建模進去？」**
可以，前提是知道哪個變數、哪些題目。風險規則不需要知道機制：它只比較「借作答時間」與「不借」對這個推論對象的差異。這就是它對 DRT 仍然有效的原因。

**反對四：「選完 η 再報告區間，不是 post-selection inference 嗎？」**
是，這也是為什麼覆蓋率只有 0.88–0.96，略低於名目的 0.95。我們如實報告，並且這個幅度遠小於完整模型的覆蓋率崩壞（0–0.68）。

**反對五：「這不就是 DIF 嗎？」**
不一樣。DIF 是作答模組的參數在群體間不同；這裡作答模組完全正確，問題出在作答時間模組把速度差轉嫁給 θ。即使作答資料完全沒有 DIF，洩漏仍然發生。這也呼應 Park et al. (2024) 的發現：DIF 題與 DRT 題重疊很少。

**反對六：「這個問題只在模擬裡出現。」**
模擬中的速度差大小（0.25–0.5 個 log 單位）與實證研究相符：Park et al. (2024) 與 Kapoor et al. (2024) 都在 PISA 中看到系統性的群體速度差。PISA 2018 的實證示範會直接回答這個問題（流程已完成，待取得資料）。

**反對七：「Hausman 檢定不是更標準嗎？」**
我們測過。在 η 很小時，D 與 V_0 − V_η 都小到和 Monte Carlo 誤差同一量級，Hausman 規則在沒有問題時也常停在 η = 0。風險規則不會發生這種失效，所以採用它。改用 RTMB 的確定性標準誤後，Hausman 規則也許能救回來，這是一個待驗證的延伸。

**反對八：「這需要特殊的 Gibbs 程式，應用者用不了。」**
不需要。第一階段只是一個「某項對數似然乘上 η」的目標函數，RTMB 大約 40 行 R，Stan／PyMC 也寫得出來。自訂 Gibbs 的價值在模擬研究的速度，以及需要第二階段推論時的巢狀 MCMC。

---

## 7. 這個框架還能解決哪些心理計量問題

作答時間只是「輔助資料來源會改寫測量」這個一般問題的一個實例。同樣的結構（受信任的測量模組 + 可能設定錯誤的輔助模組 + 共用潛在變數）出現在：

1. **Process data 與 log 資料**：點擊次數、動作序列等特徵如果和群體有關，同樣會洩漏進能力。
2. **Plausible values 的 conditioning 模型**：背景變數模組設定錯誤（例如只用主成分）時，也可以用 η 控制它對 θ 的影響，並以只用作答的結果為基準。
3. **結構方程模型的 interpretational confounding**（Levy, 2024, 2026）：加入外部變數時測量模型被改寫，SMI 提供 cut 與完整模型之間的連續選項。
4. **評分者資料、多來源評量**：可疑的評分者模組可以部分降權，而不必全丟。
5. **題目試測與線上校準**：Jewsbury et al. (2026) ✓ 用模組化貝氏保護正式量尺不受試測題影響，是同一個思路的另一個應用。

我們不在這篇論文中宣稱這些延伸都有效；它們是框架的自然推廣，也是後續研究的方向。

---

## 8. 限制，以及不能宣稱的事

1. **基準假設**：cut 以 conditioning 模型正確為前提。作答模組本身若有 DIF，基準也會偏，那是另一個問題。
2. **選擇後的覆蓋率**：0.88–0.96，略低於名目。
3. **模擬條件**：lognormal 作答時間、單一速度因素、n = 500、p = 10。題目更多時洩漏會變小；速度結構更複雜時，λ 的封閉式只是近似。
4. **不能宣稱**：
   - 事後校正是通用解法（它在 DRT 下失效）；
   - 風險規則在所有情境都達名目覆蓋率；
   - 作答時間對群體比較完全沒有價值（它提高了個人分數的精確度，而且速度本身可能就是研究對象）。

---

## 9. 這篇論文的定位

- **O2 貢獻**：為「作答時間當附帶資訊」加上一個邊界條件（群體速度差與 DRT 會造成群體偏誤），並提供一套程序（邊際調溫 SMI + 針對推論對象的風險規則 + 篩檢）。
- **目標期刊**：JEM（群體比較、公平性，需要 PISA 實證示範）或 APM（方法比較、模擬）。
- **一句話說服審稿人**：現行文獻用精確度評估作答時間的價值，但群體比較需要的是不偏；本研究證明兩者會衝突，並給出一個不需要知道衝突機制的規則。

---

## 參考文獻

✓ = 本研究中已用 Consensus 查證；† = 憑記憶引用，投稿前需核對。

- Bailey, P., et al. (2023). Expanding NAEP and TIMSS analysis to include additional variables or a new scoring model using the R package Dire. *Psych*. ✓
- Bolsinova, M., & Tijmstra, J. (2018). Improving precision of ability estimation: Getting more from response times. *British Journal of Mathematical and Statistical Psychology*. ✓
- Carmona, C. U., & Nicholls, G. K. (2020). Semi-modular inference: Enhanced learning in multi-modular models by tempering the influence of components. arXiv:2003.06804. ✓
- Duan, Q., et al. (2024). Detecting differential item functioning using response time. *Educational and Psychological Measurement*. ✓
- Frazier, D. T., Nott, D. J., et al. (2025). Posterior risk of modular and semi-modular Bayesian inference. *Journal of the American Statistical Association*. ✓
- He, Y.-H., et al. (2023). Using response time in multidimensional computerized adaptive testing. *Journal of Educational Measurement*. ✓
- Jewsbury, P. A., et al. (2026). Analytically corrected Bayesian modularization for local item calibration. arXiv:2608.11542. ✓
- Kapoor, R., et al. (2024). Differences in time usage as a competing hypothesis for observed group differences in accuracy. *Journal of Educational Measurement*. ✓
- Kern, J. L., et al. (2021). Using a response time–based expected a posteriori estimator to control for differential speededness in computerized adaptive test. *Applied Psychological Measurement*. ✓
- Levy, R., et al. (2024). Measurement and uncertainty preserving parametric modeling for continuous latent variables with discrete indicators and external variables. *Journal of Educational and Behavioral Statistics*. ✓
- Levy, R. (2026). Modular item response and structural equation modelling via measurement and uncertainty preserving parametric modelling. *British Journal of Mathematical and Statistical Psychology*. ✓
- Liu, Y., et al. (2025). A general framework for cutting feedback within modularized Bayesian inference. *Journal of the Royal Statistical Society: Series B*. ✓
- Marsman, M., et al. (2016). What can we learn from plausible values? *Psychometrika*. †
- Mislevy, R. J. (1991). Randomization-based inference about latent variables from complex samples. *Psychometrika*. †
- Molenaar, D., et al. (2024). Relating violations of measurement invariance to group differences in response times. *Psychological Methods*. ✓
- Monseur, C., & Adams, R. (2009). Plausible values: How to deal with their limitations. *Journal of Applied Measurement*. ✓
- Park, J., et al. (2024). Measurement invariance for multilingual learners using item response and response time in PISA 2018. *Educational Measurement: Issues and Practice*. ✓
- Plummer, M. (2015). Cuts in Bayesian graphical models. *Statistics and Computing*. †
- van der Linden, W. J. (2007). A hierarchical framework for modeling speed and accuracy on test items. *Psychometrika*. †
- van der Linden, W. J., Klein Entink, R. H., & Fox, J.-P. (2010). IRT parameter estimation with response times as collateral information. *Applied Psychological Measurement*. ✓
- Wu, M. (2005). The role of plausible values in large-scale surveys. *Studies in Educational Evaluation*. ✓
