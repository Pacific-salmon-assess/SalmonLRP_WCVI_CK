library(invgamma)
x<-seq(0,6,by=0.001)
invgamma0.01<-dinvgamma(x,0.01,0.01)
invgamma0.1<-dinvgamma(x,0.1,0.1)
invgamma1<-dinvgamma(x,1,1)
plot(x,invgamma0.01, typ="l", col="red", ylim=c(0,0.6), xlim=c(0,6),
     ylab="Density",xlab=expression(sigma^2), lwd=2)
lines(x,invgamma0.1, col="grey70",lty=1,lwd=2)
lines(x,invgamma1, col="red", lty=1,lwd=2)
lines(x,invgamma0.01, col="purple", lty=1,lwd=2)
lines(x, dunif(x,0.001, 10), col="black", lty=1, lwd=2)

legend(x=2.2, y=0.59, legend=c("inv.gamma(0.01,0.01)", "inv.gamma(0.1,0.1)", "inv.gamma(1,1)", "Uniform (1,10)"),
  lty=1, col=c("purple", "grey70", "red", "black"), bty="n") 

## Extra code to plot gamma

gamma0.01<-1/dgamma(x,0.01,rate=0.01)
gamma0.1<-1/dgamma(x,0.1,rate=0.1)
gamma1<-1/dgamma(x,1,rate=1)

plot(x,gamma0.01,typ="l", col="red", ylim = c(0,0.6), xlim=c(0,10))
lines(x,gamma0.1, col = "blue")
lines(x,gamma1, col = "steelblue2")