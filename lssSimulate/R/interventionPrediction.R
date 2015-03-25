generateSingleLocData2_noInt = function(seed, 
                                        population=5363500, 
                                        ThrowAwayTpt=0,
                                        beta_SE=c(-1.4,-0.055),
                                        beta_RS=c(-100000))
{    
    set.seed(seed)
    data(Kikwit95) 
    E0 = 0  
    I0 = 1 
    R0 = 0
    S0 = population - R0 - E0 - I0

    maxTpt = 150
    nTpt = maxTpt - ThrowAwayTpt
    
    X = matrix(1)

    interventionDate = as.Date("05-09-1995", "%m-%d-%Y")
    hasIntervention = (Kikwit95$Date > interventionDate)
    Z = cbind(hasIntervention*(Kikwit95$Date - interventionDate))

    X_RS = matrix(1, nrow = maxTpt)



    distMatList = -1   
    rho = -1
    gamma_ei = 1/5
    gamma_ir = 1/7
    
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

    if (sum(out$I == 0) > 100 || max(out$I) < 15)
    {
        cat("Rejecting simulated data - epidemic died out too soon.\n")
        return(generateSingleLocData2_noInt(seed + rpois(1, 1000), ThrowAwayTpt=ThrowAwayTpt, beta_SE=beta_SE))
    }
    else
    {
        # disregard post-intervention information
        orig.result = out
        maxIdx = min(which(hasIntervention)) - 1
        out = list(S_star = out$S_star[1:maxIdx,,drop=FALSE],
                   E_star = out$E_star[1:maxIdx,,drop=FALSE],
                   I_star = out$I_star[1:maxIdx,,drop=FALSE],
                   R_star = out$R_star[1:maxIdx,,drop=FALSE],
                   S = out$S[1:maxIdx,,drop=FALSE],
                   E = out$E[1:maxIdx,,drop=FALSE],
                   I = out$I[1:maxIdx,,drop=FALSE],
                   R = out$R[1:maxIdx,,drop=FALSE],
                   S0=out$S0,
                   E0=out$E0,
                   I0=out$I0,
                   R0=out$R0,
                   N=out$N[1:maxIdx,,drop=FALSE],
                   X=out$X,
                   Z=NA,
                   X_RS=out$X_rs[1:maxIdx,,drop=FALSE],
                   beta_SE=out$beta_SE[1],
                   beta_RS=out$beta_RS,
                   gamma_ei=out$gamma_ei,
                   gamma_ir=out$gamma_ir,
                   distMatList=out$distMatList,
                   effectiveTransitionSampleSize=out$effectiveTransitionSampleSize,
                   rho=out$rho,
                   orig.result=orig.result 
                   )
    }
}



simulationIntPredKernel = function(cl, genSeed, fitSeeds, ThrowAwayTpt, beta_SE)
{
    simResults = generateSingleLocData2_noInt(genSeed, ThrowAwayTpt=ThrowAwayTpt, beta_SE=beta_SE)

    fileNames = c(paste("sim", genSeed, "_1.txt", sep=""),
		  paste("sim", genSeed, "_2.txt",sep="" ),
		  paste("sim", genSeed,"_3.txt", sep=""))
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults))

    trueVals = parLapply(cl, paramsList, buildSingleLocSimInstanceIntPred)

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

    dat1=read.csv(fileNames[1])
    dat2=read.csv(fileNames[2])
    dat3=read.csv(fileNames[3])
    list(dat1=dat1,dat2=dat2,dat3=dat3,simResults=simResults)
}

buildSingleLocSimInstanceIntPred = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]

    set.seed(seed)

    DataModel = buildDataModel(simResults$I_star, type = "overdispersion", phi=1)

    priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
    ExposureModel = buildExposureModel_depricated(simResults$X, simResults$Z, 
                                       beta = c(priorBetaIntercept, rep(0, ((length(simResults$beta_SE))-1))), betaPriorPrecision = 1,
                                       nTpt=nrow(simResults$I_star))
    
    ReinfectionModel = buildReinfectionModel("SEIR")
    SamplingControl = buildSamplingControl(iterationStride=1000,
                                           sliceWidths = c(0.3,  # S_star
                                                           0.3,  # E_star
                                                           0.3, # I_star
                                                           0.3, # S0
                                                           0.3, # I0
                                                           0.01, # beta
                                                           0.01, # betaPrs
                                                           0.01, # rho
                                                           0.01, # gamma_ei
                                                           0.01, # gamma_ir
                                                           0.01 # phi
                                                          ))
    DistanceModel = buildDistanceModel(list(matrix(0)))
    TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-simResults$gamma_ei), 
                                                             1-exp(-simResults$gamma_ir), simResults$effectiveTransitionSampleSize,
                                                             simResults$effectiveTransitionSampleSize) 

    I0 = max(simResults$I_star[1:10], simResults$I_star[2], 2)
    E0 = I0
    S0 = simResults$N[1] - I0 - E0
    InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                               S0 = S0, I0 = I0, E0 = E0, reinfection=FALSE, dataType = "I_star")
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,TransitionPriors,
                         InitContainer,SamplingControl)

    res$setRandomSeed(seed)
    res$setTrace(0)
    # Burn in tuning parameters
    for (i in 1:(500))
    {
        res$simulate(10)
        res$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    for (i in 1:(50))
    {
        res$simulate(100)
        res$updateSamplingParameters(0.2, 0.05, 0.01)
    }
    res$compartmentSamplingMode = 17
    res$useDecorrelation = 50
    res$performHybridStep = 50

    # Store the model object in the global namespace of the node - can't pass these between sessions
    localModelObject <<- res
    return(list("model"=res,
                "fileName"=outFileName))    
}


runSimulationIntPred = function(cellIterations = 50, ThrowAwayTpts=c(0,6,12,24),
                          genSeed=123123, fitSeeds=c(812123,12301,5923),
                          beta_SE=c(-1.4,-0.055))
{                     

    seeds = genSeed + 100*seq(1, cellIterations)
    main.cluster = makeCluster(2)
    # simulationIntPredKernel(cl, genSeed, fitSeeds, ThrowAwayTpt, beta_SE)

    clusterExport(main.cluster, c("fitSeeds", "cellIterations", "buildSingleLocSimInstanceIntPred", 
                                  "simulationIntPredKernel", "beta_SE"),
                  envir=environment())
    outer.loop = function(seedVal){
        library(lssSimulate)
        cl = makeCluster(3, outfile = "err.txt")
        clusterExport(cl, c("buildSingleLocSimInstanceIntPred"))
        f = function(genSeed){
            library(lssSimulate)
            simulationIntPredKernel(cl, seedVal, fitSeeds, 0, beta_SE)
        }
        simResults = f(seedVal)
        save(simResults, file = paste("./simInterventionPred_", seedVal, ".Rda.bz2", sep = ""),
             compress="bzip2")
        stopCluster(cl)
        TRUE
    }
    parLapply(main.cluster, seeds, outer.loop)
}
