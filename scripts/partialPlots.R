readData = function()
{
    dat1 <<- read.csv("./sim1_1.txt")
    dat2 <<- read.csv("./sim1_2.txt")
    dat3 <<- read.csv("./sim1_3.txt")
}


plotData = function(dataName, ...)
{
    plot(dat1[[dataName]], type = "l", main = dataName, ...)
    lines(dat2[[dataName]], col = "red")
    lines(dat3[[dataName]], col = "blue")

}
