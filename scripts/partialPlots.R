readData = function(sim = 1)
{
    dat1 <<- read.csv(paste("./sim", sim, "_1.txt", sep = ""))
    dat2 <<- read.csv(paste("./sim", sim, "_2.txt", sep = ""))
    dat3 <<- read.csv(paste("./sim", sim, "_3.txt", sep = ""))
}


plotData = function(dataName, ...)
{
    plot(dat1[[dataName]], type = "l", main = dataName, ...)
    lines(dat2[[dataName]], col = "red")
    lines(dat3[[dataName]], col = "blue")

}
