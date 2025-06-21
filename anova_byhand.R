# Step 1: Simulate data
set.seed(123)
group1 <- c(8, 9, 6, 7, 10)
group2 <- c(4, 5, 6, 5, 5)
group3 <- c(9, 10, 11, 9, 12)

data <- c(group1, group2, group3)
groups <- factor(rep(c("A", "B", "C"), each = 5))

# Step 2: Compute group statistics
group_levels <- levels(groups)
k <- length(group_levels)

group_means <- sapply(group_levels, function(g) mean(data[groups == g]))
group_ns    <- sapply(group_levels, function(g) sum(groups == g))
grand_mean  <- mean(data)

# Step 3: Compute Sum of Squares
# Between-group SS
SSB <- sum(group_ns * (group_means - grand_mean)^2)

# Within-group SS
SSW <- sum(sapply(group_levels, function(g) {
  sum((data[groups == g] - mean(data[groups == g]))^2)
}))

# Total SS
SST <- sum((data - grand_mean)^2)

# Check decomposition
all.equal(SST, SSB + SSW)  # Should now return TRUE (or very close)
df_between <- k - 1
df_within <- length(data) - k
MSB <- SSB / df_between
MSW <- SSW / df_within
Fstat <- MSB / MSW
p_value <- pf(Fstat, df_between, df_within, lower.tail = FALSE)
anova_table <- data.frame(
  Source = c("Between", "Within", "Total"),
  SS = c(SSB, SSW, SST),
  df = c(df_between, df_within, df_between + df_within),
  MS = c(MSB, MSW, NA),
  F = c(Fstat, NA, NA),
  p = c(p_value, NA, NA)
)
print(anova_table)
summary(aov(data ~ groups))
