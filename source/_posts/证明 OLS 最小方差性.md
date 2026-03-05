---
title: 证明 OLS 最小方差性
---

<meta name="referrer" content="no-referrer"/>

# 证明 OLS 最小方差性

该性质是高斯-马尔可夫定理（Gauss-Markov Theorem）的核心内容，即“在满足经典线性回归模型假设的条件下，OLS 估计量是最优线性无偏估计（BLUE）”。

## 前提假设与模型设定
### 线性回归模型
设一元线性回归模型为： 
$$
Y_i = \beta_0 + \beta_1 X_i + \mu_i \quad (i=1,2,\dots,n)
$$ 

其中：
- $Y_i$ 为因变量，$X_i$ 为自变量（非随机变量）， 
- $\beta_0$（截距项）和 $\beta_1$（斜率项）为待估参数， 
- $\mu_i$ 为随机误差项，满足 **高斯-马尔可夫假设** 

### 关键假设
- **假设 1（线性性）**：$Y_i$ 是 $\beta_0, \beta_1$ 的线性函数。  
- **假设 2（无偏性）**：$E(\mu_i) = 0$，故 $E(Y_i) = \beta_0 + \beta_1 X_i$。  
- **假设 3（同方差性）**：$\text{Var}(\mu_i) = \sigma^2$（常数，与 $i$ 无关）。  
- **假设 4（无自相关性）**：$\text{Cov}(\mu_i, \mu_j) = 0$（$i \neq j$，误差项互不相关）。  


## OLS 估计量的表达式
OLS 估计量通过 **最小化残差平方和** $\sum_{i=1}^n (Y_i - \hat{\beta}_0 - \hat{\beta}_1 X_i)^2$ 得到，其表达式为：  

### 斜率估计量 $\hat{\beta}_1$ 

$$
\hat{\beta}_1 = \frac{\sum_{i = 1}^n (X_i - \bar{X})(Y_i - \bar{Y})}{\sum_{i = 1}^n (X_i - \bar{X})^2} = \sum_{i = 1}^n k_i Y_i
$$

其中：
- $\bar{X} = \frac{1}{n}\sum X_i$，$\bar{Y} = \frac{1}{n}\sum Y_i$（样本均值），  
- 为便于后续推导，我们设定 $k_i = \dfrac{X_i - \bar{X}}{\sum_{j=1}^n (X_j - \bar{X})^2}$，它是 $X_i$ 的函数。

  可见，$\hat{\beta}_1$ 是 $Y_i$ 的 **线性组合**（线性估计量）。

### 截距估计量 $\hat{\beta}_0$ 

$$
\hat{\beta}_0 = \bar{Y} - \hat{\beta}_1 \bar{X} = \sum_{i=1}^n \left( \frac{1}{n} - \bar{X} k_i \right) Y_i
$$

也是 $Y_i$ 的线性组合（线性估计量）。 

## 最小方差性的证明
需证明：**在所有线性无偏估计量中，OLS 估计量 $\hat{\beta}_1$ 的方差最小**。

### 步骤 1：定义任意线性无偏估计量  
设 $\hat{\beta}_1^*$ 是 $\beta_1$ 的任意线性估计量，可表示为：

$$
\hat{\beta}_1^* = \sum_{i=1}^n (k_i + d_i) Y_i
$$
其中 $d_i$ 是任意常数（非随机，因 $k_i$ 非随机），故 $\hat{\beta}_1^*$ 可分解为： 

$$
\hat{\beta}_1^* = \sum_{i=1}^n k_i Y_i + \sum_{i=1}^n d_i Y_i = \hat{\beta}_1 + \sum_{i=1}^n d_i Y_i
$$

### 步骤 2：利用无偏性约束推导 $d_i$ 的条件 

估计量无偏性要求 $E(\hat{\beta}_1^*) = \beta_1$，代入 $Y_i = \beta_0 + \beta_1 X_i + \mu_i$： 

$$
\begin{align*}
E(\hat{\beta}_1^*) &= E\left [ \sum (k_i + d_i)(\beta_0 + \beta_1 X_i + \mu_i) \right] \\
&= \beta_0 \sum (k_i + d_i) + \beta_1 \sum (k_i + d_i) X_i + \sum (k_i + d_i) E(\mu_i) \\
&= \beta_0 \sum (k_i + d_i) + \beta_1 \sum (k_i + d_i) X_i \quad (\text{因 } E(\mu_i)= 0)
\end{align*}
$$

由于 $\hat{\beta}_1$ 是无偏估计量，已知： $$\sum k_i = 0, \quad \sum k_i X_i = 1$$ 

因此，$\hat{\beta}_1^*$ 无偏性要求： $$\sum d_i = 0, \quad \sum d_i X_i = 0$$

### 步骤 3：计算 $\hat{\beta}_1^*$ 的方差 

利用方差 $\text{Var}(a + bZ) = \text{Var}(bZ) = b^2 \text{Var}(Z)$ 的性质，且 $\mu_i$ 无自相关：

$$
\begin{align*}
\text{Var}(\hat{\beta}_1^*) &= \text{Var}\left( \hat{\beta}_1 + \sum d_i Y_i \right) \\
&= \text{Var}(\hat{\beta}_1) + \text{Var}\left( \sum d_i Y_i \right) + 2\text{Cov}\left( \hat{\beta}_1, \sum d_i Y_i \right)
\end{align*}
$$

### 步骤 4：化简协方差项 

利用协方差的双线性性质和 $\hat{\beta}_1 = \sum k_i Y_i$，结合无自相关、同方差假定，得：

$$
\begin{split}\text{Cov}\left( \hat{\beta}_1, \sum d_i Y_i \right) &= \text{Cov}\left( \sum k_i Y_i, \sum d_j Y_j \right) \\&= \sum_{i = 1}^n \sum_{j = 1}^n k_i d_j \text{Cov}(Y_i, Y_j)
\\\xlongequal[\text{when} \ i = j,\ \text{Cov}(Y_i, Y_i) = \text{Var}(Y_i) = \text{Var}(\mu_i) = \sigma^2]{\text{when} \ i \neq j,\ \text{Cov}(Y_i, Y_j) = \text{Cov}(\mu_i, \mu_j) = 0} & = \sigma^2 \cdot \sum_{i = 1}^n k_i d_i
\end{split}
$$

代入 $k_i = \dfrac{X_i - \bar{X}}{\sum (X_j - \bar{X})^2}$，并利用无偏性得到的  $\sum d_i = 0$ 且 $\sum d_i X_i = 0$，得：

$$
\sum k_i d_i = \frac{\sum (X_i - \bar{X}) d_i}{\sum (X_j - \bar{X})^2} =\frac{\sum d_i X_i- \bar{X}\sum d_i}{\sum (X_j - \bar{X})^2}  = 0
$$

因此，协方差项为 $0$，简化得到： 
$$
\text{Var}(\hat{\beta}_1^*) = \text{Var}(\hat{\beta}_1) + \text{Var}\left( \sum d_i Y_i \right)
$$


### 步骤 5：证明方差最小性  
由于方差非负（$\text{Var}(\cdot) \geq 0$），故： 
$$
\text{Var}(\hat{\beta}_1^*) = \text{Var}(\hat{\beta}_1) + \text{Var}\left( \sum d_i Y_i \right) \geq \text{Var}(\hat{\beta}_1)
$$

当且仅当 $d_i = 0$（对所有 $i$）时，等号成立，此时 $\hat{\beta}_1^* = \hat{\beta}_1$ 同理可证， $\hat{\beta}_0$ 也满足最小方差性。


## **结论**
在满足高斯-马尔可夫假设的条件下，**OLS 估计量 $\hat{\beta}_1$ 的方差小于任何其他线性无偏估计量的方差**，即 OLS 估计量具有 **最小方差性**。  

