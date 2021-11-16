rm(list=ls())
library("tidyverse")
library("mgcv")
library('splines')
library("caret")
set.seed(143)

# Load Data
test_data <- read.csv("Moon_data.csv") %>% 
  as_tibble() %>% 
  rename(y = target)
head(test_data)

# ===================================================
# mgcv
# ===================================================
model_gam <- gam(
  y ~ s(Var1, bs = "ps", k = 10, sp = 1) + s(Var2, bs = "ps", k = 10, sp = 2), 
  data = test_data, 
  family = binomial(link = logit)
)
summary(model_gam)
model_gam$coefficients

# ===================================================
# manual
# ===================================================

# データ
var1 <- test_data$Var1
var2 <- test_data$Var2
y <- test_data$y
N <- nrow(test_data)

# 平滑化パラメータ
lambda1 <- 1
lambda2 <- 2

# ノットベクトル
knot1 <- c(-2.7156358,-2.2285810,-1.7415261,-1.2544713,-0.7674164,-0.2803616,0.2066933,
           0.6937481,1.1808030,1.6678578,2.1549126,2.6419675,3.1290223,3.6160772)
knot2 <- c(-1.37318817,-1.11812781,-0.86306744,-0.60800707,-0.35294671,-0.09788634,0.15717402,
           0.41223439,0.66729475,0.92235512,1.17741548,1.43247585,1.68753621,1.94259658)

# デザイン行列の作成
X1 <- splineDesign(knot1, var1)
QR1 <- qr(t(matrix(1, nrow = 1, ncol = N) %*% X1))
Z1 <- qr.Q(QR1, complete = TRUE)[, 2:ncol(X1)]
X1 <- X1 %*% Z1

X2 <- splineDesign(knot2, var2)
QR2 <- qr(t(matrix(1, nrow = 1, ncol = N) %*% X2))
Z2 <- qr.Q(QR2, complete = TRUE)[, 2:ncol(X2)]
X2 <- X2 %*% Z2

X <- cbind(1, X1, X2)

# ペナルティ行列
S1 <- crossprod(diff(diag(length(knot1)-4), differences = 2))
S1 <- (t(Z1) %*% S1 %*% Z1)/norm(S1)
S1 <- lambda1 * S1

S2 <- crossprod(diff(diag(length(knot2)-4), differences = 2))
S2 <- (t(Z2) %*% S2 %*% Z2)/norm(S2)
S2 <- lambda2 * S2

S <- rbind(cbind(S1, matrix(0, ncol = ncol(S2), nrow = nrow(S1))),
           cbind(matrix(0, ncol = ncol(S1), nrow = nrow(S2)), S2))
S <- cbind(0, S)
S <- rbind(0, S)

# 初期値
mu <- y*0.8 + 0.1
eta <- log(mu/(1-mu))

# P-IRLS
for (i in 1:10){
  # ①ベクトルzの計算
  z <- (((1/mu)-(1/(mu-1))) * (y - mu)) + eta
  z <- matrix(z)
  
  # ①行列Wの計算
  W <- 1/((((1/mu)-(1/(mu-1)))**2)*mu*(1-mu))
  W <- diag(W)
  
  # ②未知パラメータの推定
  SVD <- svd(W)
  M <- SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
  
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

# determinant surface
pred_data <- expand.grid(seq(min(test_data$Var1), max(test_data$Var1), length = 500), 
                         seq(min(test_data$Var2), max(test_data$Var2), length = 500))

X1_pred <- splineDesign(knot1, pred_data$Var1)
X2_pred <- splineDesign(knot2, pred_data$Var2)

X1_pred <- X1_pred %*% Z1
X2_pred <- X2_pred %*% Z2
X_pred <- cbind(1, X1_pred, X2_pred)

pred_data <- pred_data %>% 
  mutate(linear_predictor = as.numeric(X_pred %*% beta)) %>% 
  mutate(predicted_prob = 1/(1+exp(-1*linear_predictor))) %>% 
  mutate(y = ifelse(predicted_prob>=0.5, 1, 0)) %>% 
  as_tibble()

g <- ggplot(NULL) +
  geom_point(data = test_data, aes(x = Var1, y = Var2, color = factor(y)), size = 4) +
  geom_point(data = pred_data, aes(x = Var1, y = Var2, color = factor(y)), alpha = 0.03, shape = 15) +
  scale_colour_discrete("y") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
        legend.title = element_text(size = 25), legend.text = element_text(size = 20))
g

# partial residuals
predicted_prob <- 1/(1+exp(-1*as.numeric(X%*%beta)))
pearson_residuals <- (y-predicted_prob)/sqrt(predicted_prob*(1-predicted_prob))

partial_residuals1 <- as.numeric(X[,2:10] %*% beta[2:10,,drop=F]) + pearson_residuals
partial_residuals2 <- as.numeric(X[,11:19] %*% beta[11:19,,drop=F]) + pearson_residuals

line1 <- as.numeric(X_pred[,2:10] %*% beta[2:10,,drop=F])
line2 <- as.numeric(X_pred[,11:19] %*% beta[11:19,,drop=F])

g1 <- ggplot(NULL) +
  geom_point(data = tibble(Var1 = test_data$Var1, partial = partial_residuals1), aes(x = Var1, y = partial)) +
  geom_line(data = tibble(Var1 = pred_data$Var1, line = line1), aes(x = Var1, y = line), color = "orange", size = 1, alpha = 0.7) +
  xlab("Var1") + ylab("s(Var1)") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25))
g1

g2 <- ggplot(NULL) +
  geom_point(data = tibble(Var2 = test_data$Var2, partial = partial_residuals2), aes(x = Var2, y = partial)) +
  geom_line(data = tibble(Var2 = pred_data$Var2, line = line2), aes(x = Var2, y = line), color = "orange", size = 1, alpha = 0.7) +
  xlab("Var2") + ylab("s(Var2)") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25))
g2

g <- gridExtra::grid.arrange(g1, g2, nrow = 1,
                             top = textGrob("Partial residual plot",gp=gpar(fontsize=20,font=3)))

















