generateMultiLocData = function(seed, 
                                nTpt = 50,
                                 population=c(32,14,9,32,21,19,32), 
                                 beta_SE=c(-1.8650),
                                 rho = c(0.56),
                                 beta_RS=c(-100000),
                                 DM = NA)
{    
    set.seed(seed)
    E0 = rep(0, length(population))  
    I0 = c(1, rep(0, (length(population)-1))) 
    R0 = E0
    S0 = population - R0 - E0 - I0
    N = population

    nLoc = length(population)
    maxTpt = nTpt
    
    X = matrix(1, nrow = nLoc, ncol = 1)
    Z = NA

    X_RS = matrix(1, nrow = maxTpt)

    # Use overall spatial correlation parameter
    if (length(DM) == 1 && is.na(DM)){
	print("Creating Distance Matrix")
        distMatList = list((1-diag(length(N)))/length(N)) 
    }
    else{
        distMatList = DM
    }
    gamma_ei = 0.184
    gamma_ir = 0.142
    
    effectiveTransitionSampleSize = 1000
    timeIndex = 0:maxTpt

    out = generateData(seed + 1, 
                       nTpt,
                       S0,
                       E0,
                       I0,
                       R0,
                       timeIndex,
                       beta_SE,
                       beta_RS,
                       distMatList,
                       rho,
                       X,
                       Z,
                       X_RS,
                       gamma_ei,
                       gamma_ir,
                       effectiveTransitionSampleSize
                       )
    if (sum(out$I_star) < min(10, floor(sqrt(nrow(out$I_star)))) || (sum(apply(out$I_star, 2, sum) != 0) < 2 && ncol(out$I_star) != 1))
    {
	newSeed = seed + 1
        cat(paste("Epidemic died out too quickly, re-simulating with seed ", newSeed, "\n", sep = ""))
	return(generateMultiLocData(newSeed, 
                                nTpt,
                                    population,
                                    beta_SE,
                                    rho,
                                    beta_RS,
				    DM=DM))
    }
    else{
	return(out)
    }
}


simulationSmSampKernel = function(cl, genSeed, fitSeeds)
{                

    simResults = generateMultiLocData(genSeed)

    fileNames = c(paste("sim2_1", genSeed, ".txt", sep = ""),
		  paste("sim2_2", genSeed, ".txt", sep = ""),
		  paste("sim2_3", genSeed, ".txt", sep = ""))
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

buildSmSampSimInstance = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]

    set.seed(seed)

    #DataModel = buildDataModel(simResults$I_star, type = "overdispersion", params = c(500,100))
    DataModel = buildDataModel(simResults$I_star, type = "identity")


    priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
    ExposureModel = buildExposureModel_depricated(simResults$X, simResults$Z, 
                                       beta = c(priorBetaIntercept, rep(0, ((length(simResults$beta_SE))-1))), betaPriorPrecision = 1,
                                       nTpt=nrow(simResults$I_star))
    ReinfectionModel = buildReinfectionModel("SEIR")
    SamplingControl = buildSamplingControl(iterationStride=1000,
                                           sliceWidths = c(0.26,  # S_star
                                                           0.1,  # E_star
                                                           0.15, # I_star
                                                           0.22, # S0
                                                           0.24, # I0
                                                           0.8, # beta
                                                           0.2, # betaPrs
                                                           0.015, # rho
                                                           0.01, # gamma_ei
                                                           0.01, # gamma_ir
                                                           0.01 # phi
                                                          ))
    DistanceModel = buildDistanceModel(simResults$distMatList)
    TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-simResults$gamma_ei), 
                                                             1-exp(-simResults$gamma_ir), simResults$effectiveTransitionSampleSize,
                                                             simResults$effectiveTransitionSampleSize) 

    I0 = simResults$I0 
    E0 = I0*0
    S0 = simResults$N[1,] - I0 - E0
    InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                               S0 = S0, I0 = I0, E0 = E0, reinfection=FALSE, dataType = "I_star")
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,TransitionPriors,
                         InitContainer,SamplingControl)

    res$setRandomSeed(seed)
    for (i in 1:length(I0))
    {
        res$setTrace(i-1)
    }
    # Burn in tuning parameters
    for (i in 1:(200))
    {
        res$simulate(10)
        res$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    for (i in 1:(20))
    {
        res$simulate(100)
        res$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    res$compartmentSamplingMode = 17
    res$useDecorrelation = 10
    res$performHybridStep = 10

    # Store the model object in the global namespace of the node - can't pass these between sessions
    localModelObject <<- res
    return(list("model"=res,
                "fileName"=outFileName))    
}



runSimulationSmSamp = function(cellIterations = 50,
                               genSeed=123123, fitSeeds=c(812123,12301,5923))
{                     

    seeds = genSeed + 100*seq(1, cellIterations)
    main.cluster = makeCluster(2)
    clusterExport(main.cluster, c("fitSeeds", "cellIterations", "buildSmSampSimInstance", 
				"simulationSmSampKernel"), envir = environment())
    outer.loop = function(seedVal){
	    library(lssSimulate)
	    cl = makeCluster(3, outfile = "err.txt")
	    clusterExport(cl, c("buildSmSampSimInstance")) 
	    f = function(genSeed)
	    {
		library(lssSimulate)
		simulationSmSampKernel(cl, seedVal, fitSeeds)
	    }
	    simResults = f(seedVal) 
	    save(simResults, file=paste("./simSmSamp_results_0_", seedVal,".Rda.bz2", sep=""), 
		 compress="bzip2")
	    stopCluster(cl)
	    TRUE
   }
   parLapply(main.cluster, seeds, outer.loop)	
}


