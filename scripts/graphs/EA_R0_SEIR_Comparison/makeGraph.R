load("./simR0_results_0.tmp")




r01 = outData$R0Estimates[[1]][[1]]
r02 = outData$R0Estimates[[1]][[2]]

min.r0 = min(min(outData$R0Estimates[[1]][[1]]$LB), min(outData$R0Estimates[[1]][[2]]$LB))
max.r0 = max(max(outData$R0Estimates[[1]][[1]]$UB), max(outData$R0Estimates[[1]][[2]]$UB) + 0.1)
r0.ylim = c(min.r0, max.r0)
pdf(file="./EA_RO_Comparison.pdf", width = 12, height = 8)
    #par(mfrow = c(2,1))
    layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
    


    plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0: Underspecified Intensity\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l")
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = 1.5, lty = 2)
    lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2)
    lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2)

    barplot(t(outData$simResults$I_star), main = "New Cases", xlab = "Day", ylab = "Cases")
    axis(side = 1, at = seq(0, (length(r01$mean)), 50))

    plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R0: Underspecified Intensity\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l")
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = 1.5, lty = 2)
    lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2)
    lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty = 2)

dev.off()

#barplot(t(outData$simResults$I_star))

