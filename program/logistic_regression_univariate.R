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

# ロジスティック回帰実行
model_glm <- glm(y ~ x, data = test_data, family = binomial(link = "logit"))
model_glm$coefficients

# 予測確率plot
pred_data <- tibble(x = seq(min(test_data$x), max(test_data$x), length = 1000)) %>% 
  mutate(linear_predictor = model_glm$coefficients[1] + model_glm$coefficients[2]*x) %>% 
  mutate(predicted_prob = exp(linear_predictor)/(1+exp(linear_predictor)))
g1 <- ggplot(NULL) +
  geom_point(data = test_data, aes(x = x, y = y, color = factor(y))) +
  geom_line(data = pred_data, aes(x = x, y = predicted_prob), size = 1.5, alpha = 0.7) +
  xlab("x") + scale_colour_discrete("y") + ggtitle("Predicted probability ") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
        legend.title = element_text(size = 25), legend.text = element_text(size = 20), title = element_text(size = 20))
g1

# partial residuals plot
linear_predictor <- model_glm$coefficients[2]*test_data$x
pearson_residuals <- residuals(model_glm, type = "pearson")
partial_residuals <- linear_predictor + pearson_residuals

g2 <- ggplot(NULL) +
  geom_point(data = tibble(x = test_data$x, y = partial_residuals), aes(x = x, y = y)) +
  geom_line(data = tibble(x = test_data$x, y = linear_predictor), aes(x = x, y = y), color = "orange", size = 1, alpha = 0.7) +
  xlab("x") + ylab(bquote(~beta[1]*'x')) + ggtitle("Partial residual plot") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), title = element_text(size = 20))
g2

g <- gridExtra::grid.arrange(g1, g2, nrow = 1)

# confusion matrix
y_pred_prob <- predict(model_glm, test_data, type = "response")
y_pred <- c()
for (i in 1:length(y_pred_prob)){
  if (y_pred_prob[i] >= 0.5){
    y_pred[i] <- 1
  }else{
    y_pred[i] <- 0
  }
}
confusionMatrix(data = factor(y_pred), reference = factor(test_data$y))







