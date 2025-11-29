# Sasa (2014 PRL) を読むための準備メモ

## 0. このメモの目的

Sasa (2014, PRL) の「局所ギブス分布からオイラー方程式・NS方程式を導く」議論を理解するために、これまで議論したポイントを整理する。

- ミクロ密度とデルタ関数
- インデックス（添字）とベクトルの扱い
- ミクロ保存則（質量・運動量）
- アンサンブル平均と \(P_t(\Gamma)\)
- 局所ギブス分布と“局所平衡”
- 共役場 \(\lambda_\alpha(\vec r)\) の意味
- カノニカル分布の最大エントロピー導出
- 熱浴の物理的イメージと情報論的解釈

---

## 1. ミクロ密度とデルタ関数

### 1.1 粒子系の記述

- 粒子番号：\(i = 1,\dots,N\)
- 粒子の位置：\(\vec r_i(t)\)
- 粒子の運動量：\(\vec p_i(t)\)
- 観測点（空間座標）：\(\vec r = (x,y,z)\)

### 1.2 デルタ関数を用いた「ミクロ密度」

**質量ミクロ密度**：
\[
\hat\rho(\vec r,t)
= \sum_{i=1}^N m_i\,\delta\bigl(\vec r-\vec r_i(t)\bigr).
\]

イメージ：
- 各粒子は「点」に集中した密度を持つ。
- デルタ関数でその「点状の寄与」を表現し、全粒子について足し合わせている。

**重要な理解ポイント**：
- \(\hat\rho(\vec r,t)\) は「連続体の密度」ではなく、「粒子的世界から定義されたランダムな場」。
- 各実現（あるミクロ状態）では非常にスパイクだらけだが、アンサンブル平均すると滑らかなマクロ密度になる。

### 1.3 デルタ関数の時間微分

粒子の軌道 \(\vec r_i(t)\) に沿うデルタ関数
\[
\delta\bigl(\vec r - \vec r_i(t)\bigr)
\]
を時間微分すると、チェーンルールから

\[
\frac{\partial}{\partial t}\delta(\vec r-\vec r_i(t))
= -\dot{\vec r}_i(t)\cdot\nabla_{\vec r}\,\delta(\vec r-\vec r_i(t)).
\]

- \(\dot{\vec r}_i(t) = \vec v_i(t)\)：粒子の速度ベクトル
- \(\nabla_{\vec r} = (\partial/\partial x,\partial/\partial y,\partial/\partial z)\)

ここで
\[
\vec v_i(t)\cdot\nabla_{\vec r}
= v_i^x\frac{\partial}{\partial x}
+ v_i^y\frac{\partial}{\partial y}
+ v_i^z\frac{\partial}{\partial z}
\]
は、「速度方向に沿った方向微分」を表す演算子。

**スカラーか？ベクトルか？**
- \(\delta(\vec r-\vec r_i(t))\)：スカラー
- \(\nabla_{\vec r}\delta\)：勾配 → ベクトル
- \(\vec v_i\cdot\nabla_{\vec r}\delta\)：スカラー
- 両辺ともスカラーで整合している。

---

## 2. ベクトル・インデックスの記法

### 2.1 ベクトルと成分

- 位置ベクトル：
  \[
  \vec r = (r^1,r^2,r^3) = (x,y,z).
  \]
- 速度：
  \[
  \vec v_i(t) = (v_i^1,v_i^2,v_i^3).
  \]

添字 \(a,b\) は
\[
a,b = 1,2,3 \quad\leftrightarrow\quad x,y,z
\]
の「空間方向のラベル」。

### 2.2 ナブラと添字

\[
\nabla_{\vec r}
= \left(
\frac{\partial}{\partial r^1},
\frac{\partial}{\partial r^2},
\frac{\partial}{\partial r^3}
\right),
\quad
\partial_a := \frac{\partial}{\partial r^a}.
\]

内積
\[
\vec v_i\cdot\nabla_{\vec r}
= v_i^a\,\partial_a
\]
は「速度方向に沿った微分」を意味する。

---

## 3. ミクロ質量保存則（ミクロ連続の式）

### 3.1 ミクロ質量密度と質量流束

ミクロ質量密度：
\[
\hat\rho(\vec r,t)
= \sum_{i} m_i\,\delta(\vec r-\vec r_i(t)).
\]

時間微分：
\[
\frac{\partial}{\partial t}\hat\rho(\vec r,t)
= \sum_i m_i\,\partial_t \delta(\vec r-\vec r_i(t))
= -\sum_i m_i\,\vec v_i(t)\cdot\nabla_{\vec r}\,\delta(\vec r-\vec r_i(t)).
\]

ここで、**ミクロ質量流束ベクトル**を
\[
\hat{\vec J}_\rho(\vec r,t)
:= \sum_i m_i\,\vec v_i(t)\,\delta(\vec r-\vec r_i(t))
\]
と定義すると、成分表示から

\[
\vec v_i\cdot\nabla_{\vec r}\delta
= \partial_a\bigl(v_i^a\delta\bigr)
\]

を使って最終的に

\[
\boxed{
\frac{\partial}{\partial t}\hat\rho(\vec r,t)
+ \nabla_{\vec r}\cdot\hat{\vec J}_\rho(\vec r,t)
= 0
}
\]

という「**ミクロな連続の式（質量保存則）**」が得られる。

- \(\hat\rho\)：スカラー場
- \(\hat{\vec J}_\rho\)：ベクトル場
- \(\nabla\cdot\hat{\vec J}_\rho\)：スカラー
- スカラー = スカラー の形。

---

## 4. ミクロ運動量保存則の構造（Irving–Kirkwood 的スケッチ）

### 4.1 ミクロ運動量密度

運動量密度のミクロ版：
\[
\hat\pi^a(\vec r,t)
= \sum_i p_i^a(t)\,\delta(\vec r-\vec r_i(t)).
\]

同様に時間微分すると、

\[
\partial_t\hat\pi^a
= \sum_i \dot p_i^a\,\delta(\vec r-\vec r_i)
  - \partial_b\sum_i p_i^a v_i^b\,\delta(\vec r-\vec r_i).
\]

ここで
\[
\hat J_{\pi,\text{kin}}^{ab}(\vec r,t)
:= \sum_i p_i^a v_i^b\,\delta(\vec r-\vec r_i)
\]
と置くと、「運動量輸送」によるフラックスの発散として

\[
-\partial_b\hat J_{\pi,\text{kin}}^{ab}
\]
が得られる。

### 4.2 力の項と応力テンソル

力 \(\dot p_i^a\) は粒子間の相互作用から来る：
\[
\dot{\vec p}_i = \sum_{j\neq i}\vec F_{ij}.
\]

Irving–Kirkwood の計算（ペア力と線分積分）により、この「力の項」も

\[
\sum_i \dot p_i^a\,\delta(\vec r-\vec r_i)
= -\partial_b \hat J_{\pi,\text{int}}^{ab}(\vec r,t)
\]

と「発散の形」で書ける。ここで \(\hat J_{\pi,\text{int}}^{ab}\) は相互作用由来のポテンシャル応力のミクロ版。

**合計の運動量フラックス（ミクロ応力テンソル）**：
\[
\hat J_\pi^{ab}
:= \hat J_{\pi,\text{kin}}^{ab} + \hat J_{\pi,\text{int}}^{ab}.
\]

すると、ミクロ運動量保存則は

\[
\boxed{
\partial_t\hat\pi^a(\vec r,t)
+ \partial_b\hat J_\pi^{ab}(\vec r,t)
= 0
}
\]

の形になる。

- \(a\)：運動量成分
- \(b\)：流れ方向（空間方向）
- \(\hat J_\pi^{ab}\)：運動量の \(a\) 成分が、空間方向 \(b\) に流れる量 → ミクロ応力テンソル。

### 4.3 参照すべき原著論文

ミクロ応力・ミクロ保存則のきちんとした導出：

- J. H. Irving and J. G. Kirkwood,  
  *“The statistical mechanical theory of transport processes. IV. The equations of hydrodynamics,”*  
  J. Chem. Phys. **18**, 817–829 (1950).
- R. J. Hardy,  
  *“Formulas for determining local properties in molecular-dynamics simulations: Shock waves,”*  
  J. Chem. Phys. **76**, 622–628 (1982).
- D. J. Evans and G. P. Morriss,  
  *Statistical Mechanics of Nonequilibrium Liquids*, 2nd ed., Cambridge Univ. Press.

---

## 5. アンサンブル平均と \(P_t(\Gamma)\) の意味

### 5.1 ミクロ状態 \(\Gamma\) と確率分布 \(P_t(\Gamma)\)

- \(\Gamma\)：系のミクロ状態（すべての粒子の \(\{\vec r_i,\vec p_i\}\) をまとめたもの）
- 実際の世界では、ある時刻 \(t\) の \(\Gamma_t\) は**一意**に決まっている。
- しかし我々はそれを完全には知らないので、

> 「どのミクロ状態にいるか分からない」という **不確実さ** を  
> 確率分布 \(P_t(\Gamma)\) で表す。

典型的な解釈：
- 同じマクロ条件で多数のコピーを用意して考えたとき、  
  その集団の中で「\(\Gamma\) 近辺にある系の割合」が \(P_t(\Gamma)\)。
- 数学的には「エンサンブル（ensemble）」と呼ぶ。

### 5.2 期待値 \(\langle\cdot\rangle_t\)

任意のミクロ量 \(\hat A(\Gamma)\) の時刻 \(t\) における平均値：
\[
\big\langle \hat A\big\rangle_t
:= \int d\Gamma\,\hat A(\Gamma)\,P_t(\Gamma).
\]

サイコロの例：
- サイコロの出目 \(n\) に対する平均
  \[
  \langle n\rangle = \sum_n n\,P(n)
  \]
- ミクロ状態 \(\Gamma\) に対しては、その連続版
  \[
  \langle \hat A\rangle = \int d\Gamma\,\hat A(\Gamma)\,P(\Gamma)
  \]

と同じ構造。

### 5.3 マクロ保存則への移行

ミクロ保存則
\[
\partial_t\hat C_\alpha(\vec r,\Gamma_t)
+ \partial_a\hat J^{\alpha a}(\vec r,\Gamma_t) = 0
\]
を \(P_t(\Gamma)\) で平均すると、

\[
\partial_t C_t^\alpha(\vec r)
+ \partial_a J_t^{\alpha a}(\vec r) = 0,
\]

\[
C_t^\alpha(\vec r)
:= \big\langle \hat C_\alpha(\vec r,\Gamma_t)\big\rangle_t,\quad
J_t^{\alpha a}(\vec r)
:= \big\langle \hat J^{\alpha a}(\vec r,\Gamma_t)\big\rangle_t.
\]

これがマクロな保存則（後に Euler / NS になる）への出発点。

---

## 6. 局所ギブス分布（local Gibbs）と局所平衡

### 6.1 Sasa の local Gibbs 分布

初期時刻 \(t=0\) の分布として

\[
P_0(\Gamma)
= P_{\mathrm{LG}}(\Gamma;\lambda)
= \exp\Big[
 -\,\lambda^\alpha\!\cdot\!\hat C_\alpha(\Gamma)
 - \Psi(\lambda)
\Big]
\]

を仮定。

ここで
\[
\lambda^\alpha\cdot\hat C_\alpha(\Gamma)
:= \sum_\alpha\int d^3r\,\lambda_\alpha(\vec r)\,\hat C_\alpha(\vec r,\Gamma),
\]
\(\Psi(\lambda)\) は正規化条件を保証する \(\log Z\) 的な役割。

### 6.2 局所平衡の物理的イメージ

局所平衡（local equilibrium）とは：

- 空間を小さなセルに分けると、  
  各セル内ではほとんど**平衡状態**（カノニカル／グランドカノニカル）と見なせる。
- ただし、セルごとに
  - 温度 \(T(\vec r)\)
  - 流速 \(\vec u(\vec r)\)
  - 化学ポテンシャル \(\mu(\vec r)\)
  などの値は少しずつ違う。

条件：
- ミクロスケール（衝突時間・平均自由行程）に比べて、
- マクロ場の変化スケール（長さ・時間）がずっと大きい。

この条件のもとで、
「各セルごとに平衡分布だが、パラメータは位置依存」
という状態を local Gibbs が表している。

### 6.3 最大エントロピーから local Gibbs が出る

空間をセル \(i\) に離散化して考えると：

- セル \(i\) の保存量（エネ・運動量・質量など）：
  \(\hat C_{i\alpha}(\Gamma)\)
- その平均値（マクロな局所量）：
  \(C_{i\alpha} = \langle \hat C_{i\alpha}\rangle\)

**制約**：
\[
\langle \hat C_{i\alpha}\rangle = C_{i\alpha}\quad (\forall i,\alpha)
\]

この制約の下で、エントロピー
\[
S[P] = -\int d\Gamma\,P(\Gamma)\log P(\Gamma)
\]
を最大化すると、ラグランジュ乗数 \(\lambda_{i\alpha}\) を導入して

\[
P(\Gamma)
\propto \exp\Big[-\sum_{i,\alpha}\lambda_{i\alpha}\,\hat C_{i\alpha}(\Gamma)\Big]
\]
という形が出る。

連続極限で

- \(i \to \vec r\)
- \(\hat C_{i\alpha} \to \hat C_\alpha(\vec r)\,\Delta V\)
- \(\lambda_{i\alpha} \to \lambda_\alpha(\vec r)\)

とすれば、
\[
P(\Gamma)
\propto \exp\Big[-\sum_\alpha\int d^3r\,
\lambda_\alpha(\vec r)\hat C_\alpha(\vec r,\Gamma)\Big]
\]
となり、これが local Gibbs 分布。

**物理的意味**：
- 「空間各点での局所保存量 \(C^\alpha(\vec r)\) の平均値だけを固定」し、
- それ以外のミクロ情報は一切入れず「最も無知」な分布を選ぶと、
- 必然的に local Gibbs になる → それを局所平衡と呼ぶ。

---

## 7. 共役場 \(\lambda_\alpha(\vec r)\) の意味

### 7.1 一般論：制約とラグランジュ乗数

最大エントロピー原理：

- 制約：\(\langle A(\Gamma)\rangle = A_0\)
- エントロピー：\(S[P] = -\int P\log P\)

を最大化すると、

\[
P(\Gamma)\propto \exp[-\lambda A(\Gamma)]
\]

が最適解で、\(\lambda\) は制約に対するラグランジュ乗数。

### 7.2 熱力学での「共役変数」

有名な関係：

\[
\left(\frac{\partial S}{\partial U}\right)_{V,N}
= \frac{1}{T},\quad
\left(\frac{\partial S}{\partial V}\right)_{U,N}
= \frac{p}{T},\quad
\left(\frac{\partial S}{\partial N}\right)_{U,V}
= -\frac{\mu}{T}.
\]

- \(U\) に共役：逆温度 \(1/T\)
- \(V\) に共役：\(p/T\)
- \(N\) に共役：\(-\mu/T\)

のように、「保存量」と「共役変数」がセットで現れる。

### 7.3 Sasa の共役場 \(\lambda_\alpha(\vec r)\)

Sasa のエントロピー汎関数：
\[
S(C)
= \inf_\lambda\big[\lambda^\alpha\cdot C_\alpha + \Psi(\lambda)\big].
\]

Legendre 変換より
\[
\boxed{
\lambda_\alpha(\vec r)
= \frac{\delta S(C)}{\delta C^\alpha(\vec r)}
}
\]

が成り立ち、これは有限次元の
\(\partial S/\partial U = 1/T\) などの場の一般化。

よって、物理的には

- \(\lambda_0(\vec r)\) → \(\beta(\vec r)=1/T(\vec r)\)
- \(\lambda_4(\vec r)\) → \(-\beta(\vec r)\mu(\vec r)\)
- \(\lambda_a(\vec r)\) → 流速 \(\vec u(\vec r)\) に関連

と解釈できる。

> つまり：  
> **共役場 \(\lambda_\alpha(\vec r)\) は、局所の保存量密度 \(C^\alpha(\vec r)\) に共役な「温度・流速・化学ポテンシャルなどの場」**。

---

## 8. カノニカル分布の最大エントロピー導出とその物理的意味

### 8.1 状況設定

- 小さい系 \(S\) がある。
- 残りの「巨大な環境」（熱浴）とエネルギーをやり取りできる。
- 我々が知っているのは「系の平均エネルギー \(\langle H\rangle = U\) くらい」。
- それ以外のミクロ情報は何も知らない。

### 8.2 最大エントロピー問題

探すもの：確率分布 \(P(\Gamma)\)。

制約：
1. 正規化：
   \[
   \int d\Gamma\,P(\Gamma) = 1
   \]
2. 平均エネルギー：
   \[
   \int d\Gamma\,H(\Gamma) P(\Gamma) = U
   \]

最大化したい量：
\[
S[P] = -\int d\Gamma\,P(\Gamma)\log P(\Gamma).
\]

### 8.3 ラグランジュ乗数法

汎関数：
\[
\mathcal{L}[P]
= -\int d\Gamma\,P\log P
  -\alpha\left(\int d\Gamma\,P - 1\right)
  -\beta\left(\int d\Gamma\,H P - U\right).
\]

変分条件 \(\delta\mathcal{L} = 0\) より

\[
\log P(\Gamma) = -1 - \alpha - \beta H(\Gamma)
\]

→
\[
P(\Gamma) = \exp[-1-\alpha]\,\exp[-\beta H(\Gamma)]
= \frac{1}{Z(\beta)}\exp[-\beta H(\Gamma)],
\]

\[
Z(\beta) = \int d\Gamma\,\exp[-\beta H(\Gamma)].
\]

**これがカノニカル分布**：

\[
\boxed{
P_{\mathrm{can}}(\Gamma)
= \frac{1}{Z(\beta)}\exp[-\beta H(\Gamma)].
}
\]

### 8.4 β = 1/T になる理由

エントロピー：
\[
S = -\langle\log P\rangle
= \beta\langle H\rangle + \log Z(\beta)
= \beta U + \log Z(\beta).
\]

熱力学では
\[
\left(\frac{\partial S}{\partial U}\right)_{V,N}
= \frac{1}{T}.
\]

一方、上式からの計算で
\[
\frac{\partial S}{\partial U} = \beta
\]
となるので、両者を比較すると

\[
\boxed{\beta = \frac{1}{T}} \quad (\text{一般には } \beta = 1/(k_B T)).
\]

**物理的意味**：
- \(\beta\) はエネルギー制約のラグランジュ乗数。
- その熱力学的意味は「エネルギーに共役な変数」＝逆温度。

> **カノニカル分布は**  
> 「平均エネルギー U だけが分かっているときの、  
> 最も“無知”（エントロピー最大）の分布」  
> であり、その共役変数 \(\beta\) が温度に対応する。

---

## 9. 熱浴のイメージと工学的な読み替え

### 9.1 熱浴とは何か？

- 小系 \(S\) と巨大な系 \(R\)（熱浴）を考え、全体 S+R は孤立系。
- 全エネルギー：\(E_{\mathrm{tot}} = E_S + E_R = \text{const}\)。

統計的に：

- 自由度 \(N\) の系では、エネルギーの平均は \(O(N)\)、揺らぎは \(O(\sqrt{N})\)。
- 相対揺らぎ \(\Delta E/E \sim 1/\sqrt{N}\) は \(N\) が大きいと非常に小さい。

したがって：

- 小さい系 \(S\)：エネルギー揺らぎが比較的大きい。
- 巨大な浴 \(R\)：エネルギー揺らぎの相対量が極めて小さい → 温度はほぼ一定。

> 「熱浴の温度 T が変わらない」のは、  
> 制御装置が頑張っているというより、  
> **ただ大きいので相対揺らぎが無視できる**から。

### 9.2 「熱浴前提だと工学的に汎用性が低そう」への回答

カノニカル分布は

1. 「巨大な周囲（熱浴）と接した小系」という**物理モデル**から導くこともできるし、
2. 「平均エネルギーだけ知っていて、それ以外は最大エントロピー」という  
   **情報論的原理**からも導ける。

工学的には、

- 実際に熱浴がある場合は物理モデルとして理解して良いし、
- 熱浴を明示したくない場合でも  
  「平均エネルギー制約だけを使った最も自然な分布」として使える。

Sasa の local Gibbs も同様に：

- 「ミクロ衝突が速くて各セル内は平衡化している」という**物理像**と、
- 「局所保存量 \(C^\alpha(\vec r)\) の平均値だけを指定し、最大エントロピーにした分布」という**情報論的像**

の両方から解釈できる。

---

## 10. 今後のロードマップ（Sasa の本丸へ）

ここまでで整ったもの：

1. ミクロ保存則（質量・運動量・エネルギー）：
   \[
   \partial_t\hat C_\alpha + \partial_a\hat J^{\alpha a} = 0.
   \]
2. アンサンブル平均によるマクロ保存則：
   \[
   \partial_t C_t^\alpha + \partial_a J_t^{\alpha a} = 0.
   \]
3. 初期分布としての局所ギブス分布：
   \[
   P_0(\Gamma) = P_{\mathrm{LG}}(\Gamma;\lambda).
   \]
4. エントロピー汎関数と共役場：
   \[
   S(C) = \inf_\lambda[\lambda^\alpha\cdot C_\alpha + \Psi(\lambda)],
   \quad
   \lambda_\alpha(\vec r) = \frac{\delta S}{\delta C^\alpha(\vec r)}.
   \]

この上で、Sasa は

- 時間発展した分布 \(P_t\) を  
  「そのときの local Gibbs \(P_{\mathrm{LG}}(\lambda_t)\) × 補正因子 \(\exp\hat\Sigma_t\)」  
  で書き直し、
- スケール分離パラメータ \(\epsilon = (\text{ミクロ長})/(\text{マクロ長}) \ll 1\) で展開して、
  - 0次：Euler 方程式（理想流体）
  - 1次：Navier–Stokes（粘性・熱伝導、Green–Kubo で与えられる）

を導く。

---

以上。
この md をベースに、次のステップ（例えば「\(S(C)\) の 0 次近似から Euler を出す」など）をさらに分解して追っていくと、Sasa 論文全体の骨格がかなり見えやすくなるはず。
