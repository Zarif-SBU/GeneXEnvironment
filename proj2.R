library(MASS)
library(leaps)
library(knitr)
library(glmnet)

setwd("~/AMS 315 Project 2")
Dat <- read.csv('P2_422817.csv', header=TRUE)

# Look for environmental variables with the most impact
plot(Dat$E1, Dat$Y, main = "Y vs E1")
plot(Dat$E2, Dat$Y, main = "Y vs E2")
plot(Dat$E3, Dat$Y, main = "Y vs E3")
plot(Dat$E4, Dat$Y, main = "Y vs E4")

# Find the max lambda value from M_raw that will give an estimated transformation for y
M_raw <- lm(Y ~ E1 + I(E1^2) + E2 + E3 + E4, data = Dat)
plot(resid(M_raw) ~ fitted(M_raw), main='Residual Plot')
summary(M_raw)

bc <- boxcox(M_raw)
lambda <- bc$x[which.max(bc$y)] # Optimal lambda found to be .383838

# See the new transformed model
M_trans <- lm(I(Y^.38) ~ E1 + I(E1^2) + E2 + E3 + E4, data = Dat)
plot(resid(M_trans) ~ fitted(M_trans), main='Residual Plot')
summary(M_trans)

# Find correlation between Y and each IV
all_vars <- c("E1","E2","E3","E4", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", 
              "G11", "G12", "G13", "G14", "G15", "G16", "G17", "G18", "G19", "G20")
correlations <- cor((Dat$Y)^.38, Dat[,all_vars])
cor_vector <- as.vector(correlations)
names(cor_vector) <- colnames(correlations)
sorted_cors <- cor_vector[order(abs(cor_vector), decreasing = TRUE)]
top_10 <- names(sorted_cors)[1:10]
top_10

# Run model with the top 10 most correlated variables plus E1^2
M_main <- lm(I(Y^.38) ~ E1 + I(E1^2) + E2 + E4 + G7 + G11 + G12 + G13 + G14 + G15 + G19, data=Dat)
temp <- summary(M_main)
kable(temp$coefficients[ abs(temp$coefficients[,4]) <= 0.1, ], caption='Sig Coefficients')

# View bonferroni results for the most significant terms
bonferroni_results <- data.frame(
  Variable = rownames(temp$coefficients),
  Coefficient = temp$coefficients[, 1],
  Original_p = temp$coefficients[, 4],
  Bonferroni_p = temp$coefficients[, 4] * 11,
  Significant = temp$coefficients[, 4] * 11 <= .5
)

# Sort by adjusted p-value
bonferroni_results <- bonferroni_results[order(bonferroni_results$Bonferroni_p), ]
kable(bonferroni_results)

# Perform Lasso technique on the interactions
X <- model.matrix(~ (E1 + I(E1^2) + G11 + G13 + G12 + G14 + G15 + E4 + G7 + E2 + G19)^2, data = Dat)[, -1]
y <- (Dat$Y)^0.38
cv_fit <- cv.glmnet(X, y, alpha = 1)
coefs <- coef(cv_fit, s = "lambda.min")
sig_vars <- rownames(coefs)[which(coefs != 0)]
sig_interactions <- grep(":", sig_vars, value = TRUE)
sig_interactions
coefs_df <- as.matrix(coefs)
interactions <- coefs_df[grep(":", rownames(coefs_df)), , drop = FALSE]
interactions <- interactions[interaction_coefs != 0, , drop = FALSE]
sorted_interactions <- sort(abs(interactions[,1]), decreasing = TRUE)
sorted_interactions

# Model with most significant variables and interactions
M_almost <- lm(I(Y^0.38) ~ I(E1^2)+ G15 + G13 +
                G11:G15 + G11:G14, data = Dat)
summary(M_almost)
semi_final_summary <- summary(M_almost)

M_2stage <- lm( I(Y^.38) ~ (I(E1^2)+G1+G11+G15+G13)^2, data=Dat)
temp <- summary(M_2stage)
kable(temp$coefficients[ abs(temp$coefficients[,3]) >= 4, ])

M_final <- lm(I(Y^0.38) ~ I(E1^2) +
                 G11:G15, data = Dat)
summary(M_final)
