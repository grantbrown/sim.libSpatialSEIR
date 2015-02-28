spatialEstimationKernel = function(params, cl)
{                
    library(lssSimulate)
    genSeed=params$seed
    set.seed(genSeed+500)
    nTpt = params$nTpt
    population = params$population
    rho = params$rho;
    if (length(rho) != 3){
        stop("rho is of length 3 for this simulation")
    }
    DM1 = (1-diag(length(population)))/length(population)
    idx = sample(1:length(population))[1:4]
    DM2 = 0*DM1; DM2[idx[1], idx[2]] = 1; DM2[idx[2], idx[1]] = 1 
    DM3 = 0*DM1; DM3[idx[3], idx[4]] = 1; DM3[idx[4], idx[3]] = 1
    dmList = list(DM1, DM2, DM3)
    
    simResults = generateMultiLocData(genSeed, nTpt = nTpt, rho=rho, population = population, DM=dmList)

    fileNames = c(paste("sim2_1_", genSeed, ".txt", sep = ""),
		  paste("sim2_2_", genSeed, ".txt", sep = ""),
		  paste("sim2_3_", genSeed, ".txt", sep = ""))
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults))

    trueVals = parLapply(cl, paramsList, buildSmSampSimInstance)

    iterationParams = list(list(200000, 1000, 0.2, 0.05, 0.1),
                           list(200000, 1000, 0.2, 0.05, 0.1),
                           list(200000, 1000, 0.2, 0.05, 0.1))
    conv = FALSE
    while (!conv)
    {
        cat("Not converged, adding iterations...\n")
        parLapply(cl, iterationParams, additionalIterations)
        conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], maxVal=1.02)
    }

    #computeSim2Results(fileNames[1], fileNames[2], fileNames[3], simResults)

    dat1=read.csv(fileNames[1])
    dat2=read.csv(fileNames[2])
    dat3=read.csv(fileNames[3])
    list(dat1=dat1,dat2=dat2,dat3=dat3,simResults=simResults)
}



runSimulationSpatialEstimation = function(cellIterations = 50, 
                                          nTpt = 50, rho = c(0.1, 0.25, 0.5), 
                                          population = rep(10000, 10),
                                          genSeed = 123123,
                                          fitSeeds=c(812123,12301,5923)
                                          ){
    seeds = genSeed + 100*seq(1, cellIterations)
    params = lapply(seeds, function(x){list(seed = x, nTpt=nTpt, population=population, rho=rho)})
    main.cluster = makeCluster(2)
    clusterExport(main.cluster, c("fitSeeds", "buildSmSampSimInstance", 
				"simulationSmSampKernel"), envir = environment())
    outer.loop = function(paramVal){
	    fname = paste("./simLgSamp_results_0_", paramVal[["seed"]],".Rda.bz2", sep="")
	    if (file.exists(fname)){
	    	return(FALSE)
	    }
	    library(lssSimulate)
	    cl = makeCluster(3, outfile = "err.txt")
	    clusterExport(cl, c("fitSeeds", "buildSmSampSimInstance"), envir=environment()) 
	    f = function(paramVals)
	    {
    		library(lssSimulate)
	        spatialEstimationKernel(paramVals,cl)
	    }
	    simResults = f(paramVal) 
	    save(simResults, file=fname, 
		 compress="bzip2")
	    stopCluster(cl)
	    TRUE
   }
   parLapplyLB(main.cluster, params, outer.loop)	
   stopCluster(main.cluster)
}

