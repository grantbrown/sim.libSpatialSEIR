generateSingleLocData.eaR0 = function(seed, 
                                  population=5363500, 
                                  ThrowAwayTpt=0,
                                  beta_SE=c(-1.5,-0.055),
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
    # use the real data here.
    out$I_star = matrix(Kikwit95$Count[1:nTpt], ncol=1)

    if (sum(out$I == 0) > 100)
    {
        cat("Rejecting simulated data - epidemic died out too soon.\n")
        return(generateSingleLocData.eaR0(seed + rpois(1, 1000), ThrowAwayTpt=ThrowAwayTpt))
    }
    out
}

hpdEst = function(vec)
{
    dens = density(vec)
    dens$x[which.max(dens$y)]
}

estimateR0 = function(params)
{
    extraIterations = params[[1]]
    batchSize = params[[2]]

    R0 = array(0, dim = c(nrow(localModelObject$I_star), 
                          ncol(localModelObject$I_star), extraIterations))
    empiricalR0 = R0
    for (i in 1:(extraIterations))
    {
        localModelObject$simulate(batchSize)
        for (j in (0:(nrow(localModelObject$I_star) - 1)))
        {
            R0[j,,i] = localModelObject$estimateR0(j)
            empiricalR0[j,,i] = apply(localModelObject$getIntegratedGenerationMatrix(j), 1, sum)
        }
    }

    R0Mean = apply(R0, 1:2, mean)
    R0LB = apply(R0, 1:2, quantile, probs = 0.05)
    R0UB = apply(R0, 1:2, quantile, probs = 0.95)

    empiricalR0Mean = apply(empiricalR0, 1:2, hpdEst)
    empiricalR0LB = apply(empiricalR0, 1:2, quantile, probs = 0.05)
    empiricalR0UB = apply(empiricalR0, 1:2, quantile, probs = 0.95)

    outList = list(R0 = list(mean=R0Mean,
                             LB=R0LB,
                             UB=R0UB),
                   empiricalR0 = list(mean = empiricalR0Mean,
                                      LB = empiricalR0LB,
                                      UB = empiricalR0UB)
                   )
    return(outList)
}

simulationKernel.eaR0 = function(cl, genSeed, fitSeeds, ThrowAwayTpt, underspecified)
{
    simResults = generateSingleLocData.eaR0(genSeed, ThrowAwayTpt=ThrowAwayTpt)

    fileNames = c("simR0_1.txt", "simR0_2.txt", "simR0_3.txt")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults = simResults, underspecified = underspecified),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults = simResults, underspecified = underspecified),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults = simResults, underspecified = underspecified))

    trueVals = parLapply(cl, paramsList, buildSingleLocSimInstance.eaR0)

    iterationParams = list(list(10000, 1000, 0.2, 0.05, 0.1),
                           list(10000, 1000, 0.2, 0.05, 0.1),
                           list(10000, 1000, 0.2, 0.05, 0.1))
    conv = FALSE
    while (!conv)
    {
        cat("Not converged, adding iterations...\n")
        parLapply(cl, iterationParams, additionalIterations)
        conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], maxVal = 1.02)
    }
    R0Params = list(list(1000,1000),
                    list(1000,1000),
                    list(1000,1000))


    R0Estimates = parLapply(cl, R0Params, estimateR0) 
    return(list(R0Estimates = R0Estimates, simResults = simResults))
}

buildSingleLocSimInstance.eaR0 = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]
    underspecified = params[[4]]

    set.seed(seed)

    I_star = matrix(simResults$I_star, ncol = 1)
    DataModel = buildDataModel(I_star, type = "overdispersion", params = c(10000,100))
    #DataModel = buildDataModel(simResults$I_star, type = "identity")


    priorBetaIntercept = log(mean(-log(1-(I_star/(simResults$N))))) 
    if (underspecified)
    {
        ExposureModel = buildExposureModel(simResults$X, Z=NA, 
                                           beta = c(priorBetaIntercept), betaPriorPrecision = 0.1,
                                           nTpt=nrow(I_star))
    }
    else
    {
        ExposureModel = buildExposureModel(simResults$X, Z=simResults$Z, 
                                           beta = c(priorBetaIntercept, 0), betaPriorPrecision = 0.1,
                                           nTpt=nrow(I_star))
    }
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
    DistanceModel = buildDistanceModel(list(matrix(0)))
    TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-simResults$gamma_ei), 
                                                             1-exp(-simResults$gamma_ir), simResults$effectiveTransitionSampleSize,
                                                             simResults$effectiveTransitionSampleSize) 

    I0 = max(I_star[1:10], I_star[2], 2)
    E0 = I0
    S0 = simResults$N[1] - I0 - E0
    InitContainer = buildInitialValueContainer(I_star, simResults$N, 
                                               S0 = S0, I0 = I0, E0 = E0, reinfection=FALSE, dataType = "I_star")
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,TransitionPriors,
                         InitContainer,SamplingControl)

    res$setRandomSeed(seed)
    res$setTrace(0)
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


runSimulation.eaR0 = function(ThrowAwayTpt=0, genSeed=123123, fitSeeds=c(812123,12301,5923), underspecified=FALSE)
{                     
    cl = makeCluster(3, outfile = "err.txt")
    print("Cluster Created")
    clusterExport(cl, c("buildSingleLocSimInstance.eaR0"))
    print("Variables Exported.") 
 
    simResults = simulationKernel.eaR0(cl, genSeed, fitSeeds,ThrowAwayTpt, underspecified)
    
    outData = simResults 
    specString = ifelse(underspecified, "underspec", "correct")
    save(outData, file=paste("./simR0_results_", specString, "_", ThrowAwayTpt, ".tmp", sep=""))

    print("Results obtained")
    stopCluster(cl)
    TRUE
}
