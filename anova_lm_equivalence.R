# Step 1: Define the data
group <- factor(rep(c("A", "B", "C"), each = 3))
y <- c(4, 5, 6,   # Group A
       7, 8, 9,   # Group B
       10, 11, 12)  # Group C

# Step 2: ANOVA by hand ---------------------------------------

# Group means and grand mean
group_means <- tapply(y, group, mean)
grand_mean <- mean(y)

# SSB: Between-group sum of squares
n_per_group <- table(group)
SSB <- sum(n_per_group * (group_means - grand_mean)^2)

# SSW: Within-group sum of squares
SSW <- sum(tapply(y, group, function(g) sum((g - mean(g))^2)))

# SST: Total sum of squares
SST <- sum((y - grand_mean)^2)

# Degrees of freedom
df_between <- length(unique(group)) - 1
df_within <- length(y) - length(unique(group))

# Mean squares and F
MSB <- SSB / df_between
MSW <- SSW / df_within
F_anova <- MSB / MSW

# Step 3: Linear model manually --------------------------------

# Construct design matrix for treatment-coded regression
X <- model.matrix(~ group)  # Intercept + groupB + groupC
# Fit beta: (X^T X)^-1 X^T y
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

# Predicted values and residuals
y_hat <- X %*% beta_hat
residuals <- y - y_hat

# RSS (same as SSW)
RSS <- sum(residuals^2)
# ESS (explained sum of squares)
ESS <- sum((y_hat - grand_mean)^2)

# F-stat from regression
MSR <- ESS / df_between
MSE <- RSS / df_within
F_lm <- MSR / MSE

# Step 4: Compare Results --------------------------------------

cat("ANOVA F-statistic:", F_anova, "\n")
cat("Regression F-statistic:", F_lm, "\n")

# Confirm equivalence
all.equal(F_anova, F_lm)  # Should return TRUE

# Optional: Compare to built-in
summary(aov(y ~ group))
summary(lm(y ~ group))
