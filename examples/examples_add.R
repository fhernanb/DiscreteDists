# Poisson expected value
add(f=function(x, lambda) x*dpois(x, lambda), lower=0, upper=Inf,
    lambda=7.5)

# Binomial expected value
add(f=function(x, size, prob) x*dbinom(x, size, prob), lower=0, upper=20,
    size=20, prob=0.5)

# Examples with infinite series
add(f=function(x) 0.5^x, lower=0, upper=100) # Ans=2
add(f=function(x) (1/3)^(x-1), lower=1, upper=Inf) # Ans=1.5
add(f=function(x) 4/(x^2+3*x+2), lower=0, upper=Inf, abs.tol=0.001) # Ans=4.0
add(f=function(x) 1/(x*(log(x)^2)), lower=2, upper=Inf, abs.tol=0.000001) # Ans=2.02
add(f=function(x) 3*0.7^(x-1), lower=1, upper=Inf) # Ans=10
add(f=function(x, a, b) a*b^(x-1), lower=1, upper=Inf, a=3, b=0.7) # Ans=10
add(f=function(x, a=3, b=0.7) a*b^(x-1), lower=1, upper=Inf) # Ans=10

