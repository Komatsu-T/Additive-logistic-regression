rm(list=ls())
library("tidyverse")
library("mgcv")
library('splines')
library("caret")
library("gridExtra")
set.seed(143)

# データ読み込み
test_data <- read.csv("logistic_univariate.csv") %>% 
  as_tibble()
head(test_data)

# ===================================================
# additive logistic regression: mgcvを用いる
# ===================================================
model_gam <- gam(y ~ s(x, bs = "ps", sp = 5), data = test_data, family = binomial(link = logit))
summary(model_gam)
model_gam$coefficients

# 予測確率plot
pred_data <- tibble(x = seq(min(test_data$x), max(test_data$x), length = 1000))
pred_data <- pred_data %>% 
  mutate(predicted_prob = predict(model_gam, pred_data, type = "response"))
g1 <- ggplot(NULL) +
  geom_point(data = test_data, aes(x = x, y = y, color = factor(y))) +
  geom_line(data = pred_data, aes(x = x, y = predicted_prob), size = 1.5, alpha = 0.7) +
  xlab("x") + scale_colour_discrete("y") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
        legend.title = element_text(size = 25), legend.text = element_text(size = 20))
g1

# partial residuals plot
plot <- plot(model_gam, residual = TRUE)
g2 <- ggplot(NULL) +
  geom_point(data = tibble(x = plot[[1]]$raw, y = plot[[1]]$p.resid), aes(x = x, y = y)) +
  geom_line(data = tibble(x = plot[[1]]$x, y = plot[[1]]$fit), aes(x = x, y = y), color = "orange", size = 1, alpha = 0.7) +
  xlab("x") + ylab("s(x)") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25))
g2

# confusion matrix
y_pred_prob <- predict(model_gam, test_data, type = "response")
y_pred <- c()
for (i in 1:length(y_pred_prob)){
  if (y_pred_prob[i] >= 0.5){
    y_pred[i] <- 1
  }else{
    y_pred[i] <- 0
  }
}
m1 <- confusionMatrix(data = factor(y_pred), reference = factor(test_data$y))
m1

# ===================================================
# additive logistic regression: mgcvを用いずに推定
# ===================================================

# データ
x <- test_data$x
y <- test_data$y
N <- nrow(test_data)

# 平滑化パラメータ
lambda <- 5

# ノットベクトル
knot <- c(-9.3042857,-7.8728571,-6.4414286,-5.0100000,-3.5785714,-2.1471429,-0.7157143,
          0.7157143,2.1471429,3.5785714,5.0100000,6.4414286,7.8728571,9.3042857)

# デザイン行列の作成
X <- splineDesign(knot, x)
QR <- qr(t(matrix(1, nrow = 1, ncol = N) %*% X))
Z <- qr.Q(QR, complete = TRUE)[, 2:ncol(X)]
X <- X %*% Z
X <- cbind(1, X)

# ペナルティ行列（2階差分ペナルティ）
S <- crossprod(diff(diag(length(knot)-4), differences = 2))
S <- (t(Z) %*% S %*% Z)/norm(S)
S <- lambda * S
S <- cbind(0, S)
S <- rbind(0, S)

# パラメータ推定
# 初期値
mu <- y*0.8 + 0.1
eta <- log(mu/(1-mu))

# P-IRLS
for (i in 1:10){
  # ①ベクトルzの計算
  z <- ((y-mu)/(mu*(1-mu))) + eta
  z <- matrix(z)
  
  # ①行列Wの計算
  W <- mu*(1-mu)
  W <- diag(W)
  
  # ②未知パラメータの推定
  M <- sqrt(W)
  
  SVD <- svd(S)
  D <- SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
  
  XD <- rbind(X, D)
  z0 <- rbind(z, matrix(0, ncol = 1, nrow = nrow(D)))
  WI <- rbind(cbind(M, matrix(0, ncol = nrow(D), nrow = nrow(M))),
               cbind(matrix(0, ncol = ncol(M), nrow = nrow(D)), diag(1, nrow = nrow(D), ncol = nrow(D))))
  
  QR <- qr(WI%*%XD)
  beta <- backsolve(qr.R(QR), qr.qty(QR, WI%*%z0))

  # ③eta, muの更新
  eta <- as.numeric(X%*%beta)
  mu <- 1/(1+exp(-1*eta))
}

as.numeric(beta)
beta
model_gam$coefficients
round(as.numeric(beta), 6) == round(model_gam$coefficients, 6)

# 予測確率plot
X_pred <- splineDesign(knot, seq(-5, 5, length = 1000))
X_pred <- X_pred %*% Z
X_pred <- cbind(1, X_pred)

linear_predictor <- as.numeric(X_pred %*% beta)
linear_predictor_partial <- as.numeric(X_pred[,2:10] %*% beta[2:10,,drop=F])
predicted_probability <- 1/(1 + exp(-1*linear_predictor))

# Predicted prob plot
g3 <- ggplot(NULL) +
  geom_point(data = test_data, aes(x = x, y = y, color = factor(y))) +
  geom_line(data = tibble(x = seq(-5, 5, length = 1000), y = predicted_probability), aes(x = x, y = y), size = 1.5, alpha = 0.7) +
  xlab("x") + scale_colour_discrete("y") + ggtitle("Predicted probability ") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
        legend.title = element_text(size = 25), legend.text = element_text(size = 20), title = element_text(size = 20))
g3

# Partial residual plot
predicted_prob <- 1/(1+exp(-1*as.numeric(X%*%beta)))
pearson_residuals <- (y-predicted_prob)/sqrt(predicted_prob*(1-predicted_prob))
partial_residuals <- as.numeric(X[,2:10] %*% beta[2:10,,drop=F]) + pearson_residuals

g4 <- ggplot(NULL) +
  geom_point(data = tibble(x = x, y = partial_residuals), aes(x = x, y = y)) +
  geom_line(data = tibble(x = seq(-5, 5, length = 1000), y = linear_predictor_partial), aes(x = x, y = y), color = "orange", size = 1, alpha = 0.7) +
  xlab("x") + ylab("s(x)") + ggtitle("Partial residual plot") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), title = element_text(size = 20))
g4

g <- gridExtra::grid.arrange(g3, g4, nrow = 1)

# confusion matrix
predicted_prob <- 1/(1+exp(-1*as.numeric(X%*%beta)))
y_pred <- c()
for (i in 1:length(predicted_prob)){
  if (predicted_prob[i] >= 0.5){
    y_pred[i] <- 1
  }else{
    y_pred[i] <- 0
  }
}
m2 <- confusionMatrix(data = factor(y_pred), reference = factor(test_data$y))
m2










