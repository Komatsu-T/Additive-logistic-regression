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

# ロジスティック回帰実行
model_glm <- glm(y ~ Var1 + Var2, data = test_data, family = binomial(link = "logit"))
model_glm$coefficients

# prediction
pred_data <- expand.grid(seq(min(test_data$Var1), max(test_data$Var1), length = 500), 
                         seq(min(test_data$Var2), max(test_data$Var2), length = 500)) %>% 
  mutate(linear_predictor = model_glm$coefficients[1] + model_glm$coefficients[2]*Var1 + model_glm$coefficients[3]*Var2) %>% 
  mutate(predict_prob = 1/(1+exp(-1*linear_predictor))) %>%
  mutate(y = ifelse(predict_prob>=0.5, 1, 0)) %>% 
  as_tibble()

# determinant surface
g <- ggplot(NULL) +
  geom_point(data = test_data, aes(x = Var1, y = Var2, color = factor(y)), size = 4) +
  geom_point(data = pred_data, aes(x = Var1, y = Var2, color = factor(y)), alpha = 0.03, shape = 15) +
  scale_colour_discrete("y") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
        legend.title = element_text(size = 25), legend.text = element_text(size = 20))
g

















