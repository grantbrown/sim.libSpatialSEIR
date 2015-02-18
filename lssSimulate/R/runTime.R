runTimeSimKernel = function(dataRow){
    library(lssSimulate)
    genSeed=dataRow[1]
    nTpt = dataRow[2]
    nLoc = dataRow[3]
    simResults = generateMultiLocData(genSeed, nTpt = nTpt, population = rep(10000, nLoc))
    params = list(seed=genSeed + nTpt, outFileName = paste("rtsim_",genSeed+nTpt , "_", nLoc, ".txt", sep = ""),
         simResults=simResults)
    buildSmSampSimInstance(params)
    localModelObject$compartmentSamplingMode = 1
    localModelObject$useDecorrelation = 0
    localModelObject$performHybridStep = 0
    Time1 = system.time(localModelObject$simulate(5000))[3]
    localModelObject$compartmentSamplingMode = 17
    localModelObject$useDecorrelation = 0
    localModelObject$performHybridStep = 0
    Time2 = system.time(localModelObject$simulate(5000))[3]
    localModelObject$compartmentSamplingMode = 17
    localModelObject$useDecorrelation = 10
    localModelObject$performHybridStep = 10
    Time3 = system.time(localModelObject$simulate(5000))[3]
    return(c(Time1, Time2, Time3))
}

runSimulationIterTime = function(nTpt = c(seq(10,1000,10)), nLoc = 1,genSeed = 123123){
    seeds = genSeed + 100 + nTpt
    paramMat = cbind(seeds, nTpt, nLoc)
    cl = makeCluster(4)
    clusterExport(cl, c("generateMultiLocData", "buildSmSampSimInstance"))
    result = parApply(cl, paramMat, 1, runTimeSimKernel)
    stopCluster(cl)
}

