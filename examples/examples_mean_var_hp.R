# Example 1

# Theoretical values
mean_var_hp(mu=5.5, sigma=0.1)

# Using simulated values
y <- rHYPERPO(n=1000, mu=5.5, sigma=0.1)
mean(y)
var(y)


# Example 2

# Theoretical values
mean_var_hp2(mu=5.5, sigma=1.9)

# Using simulated values
y <- rHYPERPO2(n=1000, mu=5.5, sigma=1.9)
mean(y)
var(y)
