# Example 1
x_max <- 30
probs1 <- dHYPERPO(x=0:x_max, sigma=0.01, mu=5)
probs2 <- dHYPERPO(x=0:x_max, sigma=0.50, mu=5)
probs3 <- dHYPERPO(x=0:x_max, sigma=1.00, mu=5)
# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probabilities for hyper-Poisson",
     ylim=c(0, 0.20))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="green4")
legend("topright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("sigma=0.01, mu=5",
                "sigma=0.50, mu=5",
                "sigma=1.00, mu=5"))
