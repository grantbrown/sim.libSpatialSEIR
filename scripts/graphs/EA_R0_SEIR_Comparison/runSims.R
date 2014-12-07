library(lssSimulate)

a = runSimulation.eaR0(underspecified=TRUE)
dat1 = read.csv("./simR0_1.txt")
underSpecItrs = max(dat1$Iteration)
underSpecTime = max(dat1$Time)
a = runSimulation.eaR0(underspecified=FALSE)
dat1 = read.csv("./simR0_1.txt")
correctSpecItrs = max(dat1$Iteration)
correctSpecTime = max(dat1$Time)

out = matrix(c(underSpecItrs, correctSpecItrs, underSpecTime, correctSpecTime), ncol = 2)
colnames(out) = c("Iterations", "Time")
write.csv(out, file="./IterationsTime.csv")


#betaDens = density(exp(dat1$BetaP_SE_0[floor(nrow(dat1)/2):nrow(dat1)]))
#plot(betaDens, 
#     main = "Intensity Intercept\n Posterior Density", xlab = expression(beta))
#maxIdx = which.max(betaDens$y)
#lines(x = rep(betaDens$x[maxIdx],2), c(0, betaDens$y[maxIdx]), lty = 2)
#density(exp(dat1$BetaP_SE_0[floor(nrow(dat1)/2):nrow(dat1)])
#text(x=0.25, y = betaDens$y[maxIdx], labels=paste("HPD Estimate: ", 
#                                                  round(betaDens$x[maxIdx],2),sep=""))


