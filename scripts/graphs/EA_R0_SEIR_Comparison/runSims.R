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


