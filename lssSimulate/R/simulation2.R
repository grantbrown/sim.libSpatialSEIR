generateSingleLocData2 = function(seed, 
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
        return(generateSingleLocData2(seed + rpois(1, 1000), ThrowAwayTpt=ThrowAwayTpt, beta_SE=beta_SE))
    }
    out
}


simulation2Kernel = function(cl, genSeed, fitSeeds, ThrowAwayTpt, beta_SE)
{
    #TODO: Vary starting linear predictor parameters on each iteration 
    simResults = generateSingleLocData2(genSeed, ThrowAwayTpt=ThrowAwayTpt, beta_SE=beta_SE)

    fileNames = c("sim2_1.txt", "sim2_2.txt", "sim2_3.txt")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults))

    trueVals = parLapply(cl, paramsList, buildSingleLocSimInstance2)

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

buildSingleLocSimInstance2 = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]

    set.seed(seed)

    DataModel = buildDataModel(simResults$I_star, type = "overdispersion", params = c(500,100))
    #DataModel = buildDataModel(simResults$I_star, type = "identity")


    priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
    ExposureModel = buildExposureModel_depricated(simResults$X, simResults$Z, 
                                       beta = c(priorBetaIntercept, rep(0, ((length(simResults$beta_SE))-1))), betaPriorPrecision = 1)
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


computeSim2Results = function(fileName1, fileName2, fileName3, trueData)
{
    contains = function(trueVal, variableName, summarylist)
    {
        cat(paste("Checking for: ", variableName, "\n", sep = ""))
        sl = summarylist[[2]]
        variableIdx = which(rownames(sl) == variableName)
        if (length(variableIdx) == 0)
        {
            return(NA)
        }
        bounds = sl[variableIdx,][c(1,5)]
        trueVal < bounds[2] && trueVal > bounds[1]
    }

    getMean = function(variableName, summaryList)
    {
        sl = summaryList[[1]]
        variableIdx = which(rownames(sl) == variableName)
        if (length(variableIdx) == 0){
            return(NA)
        }
        sl[variableIdx, 1]
    }

    dat1 = read.csv(fileName1)
    dat2 = read.csv(fileName2)
    dat3 = read.csv(fileName3)

    dat1_full = dat1[floor(nrow(dat1)/2):nrow(dat1),]
    dat2_full = dat2[floor(nrow(dat2)/2):nrow(dat2),]
    dat3_full = dat3[floor(nrow(dat3)/2):nrow(dat3),]

    iteration = dat1$Iteration
    time = dat1$Time
    
    dat1 = dat1_full[,1:(ncol(dat1_full) - 2)]
    dat2 = dat2_full[,1:(ncol(dat2_full) - 2)]
    dat3 = dat3_full[,1:(ncol(dat3_full) - 2)]
    nonZeroVar = ((apply(dat1, 2, var) != 0) | (apply(dat2, 2, var) != 0) | (apply(dat3, 2, var) != 0))
    mcl = mcmc.list(as.mcmc(dat1[,nonZeroVar]),
                    as.mcmc(dat2[,nonZeroVar]),
                    as.mcmc(dat3[,nonZeroVar]))

    mcl.summary = summary(mcl)

    # Calculate parameter summaries
    trueBeta = trueData$beta_SE
    trueGamma_ei = trueData$gamma_ei
    trueGamma_ir = trueData$gamma_ir

    betaNames = paste("BetaP_SE_", seq(0, (length(trueBeta)-1)), sep = "")

    allCompareVariableNames = c(betaNames, "gamma_ei", "gamma_ir")
    allCompareVariableValues = c(trueBeta, trueGamma_ei, trueGamma_ir)

    resultsMatrix1 = matrix(0, nrow = length(allCompareVariableNames), ncol=5)
    colnames(resultsMatrix1) = c("CI_Contains", "CI_contains_0", "SignMatches","Bias", "BiasPct")
    rownames(resultsMatrix1) = allCompareVariableNames
    for (varNum in 1:length(allCompareVariableNames))
    {
        varName = allCompareVariableNames[varNum]
        varVal = allCompareVariableValues[varNum]
        bias = getMean(varName, mcl.summary) - varVal
        pctBias = bias/(abs(varVal))*100
        CIcontains = contains(varVal, varName, mcl.summary)
        CIcontainsZero = contains(0, varName, mcl.summary)
        estimateSignMatches = sign(varVal) ==  sign(getMean(varName, mcl.summary))
        resultsMatrix1[varNum,] = c(CIcontains, CIcontainsZero, estimateSignMatches, bias, pctBias)
    }

    # calculate estimation performance
    print("Process Compartments") 
    tpts = nrow(trueData$S) 
    Snames = paste("S_0_", seq(0,  (tpts-1)), sep = "")
    Enames = paste("E_0_", seq(0,  (tpts-1)), sep = "")
    Inames = paste("I_0_", seq(0,  (tpts-1)), sep = "")
    Rnames = paste("R_0_", seq(0,  (tpts-1)), sep = "")

    Smatrix = matrix(0, ncol = 3, nrow = length(Snames))
    rownames(Smatrix) = Snames
    colnames(Smatrix) = c("CI_Contains", "Bias", "BiasPct")

    Ematrix = matrix(0, ncol = 3, nrow = length(Enames))
    rownames(Ematrix) = Enames
    colnames(Ematrix) = c("CI_Contains", "Bias", "BiasPct")

    Imatrix = matrix(0, ncol = 3, nrow = length(Inames))
    rownames(Imatrix) = Inames
    colnames(Imatrix) = c("CI_Contains", "Bias", "BiasPct")

    Rmatrix = matrix(0, ncol = 3, nrow = length(Rnames))
    rownames(Rmatrix) = Rnames
    colnames(Rmatrix) = c("CI_Contains", "Bias", "BiasPct")

    processCompartment = function(variableName)
    {
        cat(paste("Processing: ", variableName, "\n", sep = ""))
        var = c(dat1_full[[variableName]], 
                dat2_full[[variableName]], 
                dat3_full[[variableName]]) 
        tmp = strsplit(variableName, "_")[[1]]

        variableType = tmp[1]
        variableIndex = as.numeric(tmp[length(tmp)]) + 1
        trueVal = trueData[[variableType]][variableIndex]
        
        meanVal = mean(var)
        bounds = quantile(var, c(0.025, 0.975))
        CIcontains = (trueVal > min(bounds) & trueVal < max(bounds))
        bias = meanVal - trueVal
        biasPct = bias/abs(trueVal)*100
        c(CIcontains, bias, biasPct) 
    }


    for (i in 1:length(Snames))
    {
        Smatrix[i,] = processCompartment(Snames[i])  
        Ematrix[i,] = processCompartment(Enames[i])  
        Imatrix[i,] = processCompartment(Inames[i])  
        Rmatrix[i,] = processCompartment(Rnames[i])  
    }
 
    return(list("params"=resultsMatrix1,
                "compartments"=list(Smatrix, Ematrix, Imatrix,Rmatrix),
                "iterations"=max(iteration), 
                "time"=max(time)))
}


runSimulation2 = function(cellIterations = 50, ThrowAwayTpts=c(0,6,12,24),
                          genSeed=123123, fitSeeds=c(812123,12301,5923),
                          beta_SE=c(-1.4,-0.055))
{                     
    cl = makeCluster(3, outfile = "err.txt")
    print("Cluster Created")
    clusterExport(cl, c("buildSingleLocSimInstance2"))
    print("Variables Exported.") 
 
    for (ThrowAwayTpt in ThrowAwayTpts)
    {
        f = function(genSeed)
        {
            simulation2Kernel(cl, genSeed, fitSeeds + genSeed,ThrowAwayTpt, beta_SE)
        }
        #simResults = lapply(genSeed + seq(1, cellIterations), f)
	itrSeeds = genSeed + seq(1, cellIterations)
	i = 1
	for (itrSeed in itrSeeds){
	    result = f(itrSeed)
	    save(result, file=paste("./sim2_results_", ThrowAwayTpt, "_", i, ".Rda.bz2", sep=""), 
        	 compress="bzip2")
	    i=i+1
	}
	
        #biasResults = simResults[[1]]$params
        #compartmentResults = simResults[[1]]$compartments
        #timeResult = simResults[[1]]$time
        #iterationResult = simResults[[1]]$iterations
        #if (length(simResults) > 1)
        #{
        #    for (i in 2:length(simResults))
        #    {
        #        biasResults = biasResults + simResults[[i]]$params
        #        timeResult = timeResult + simResults[[i]]$time
        #        iterationResult = iterationResult + simResults[[i]]$iterations
        #        for (j in 1:4)
        #        {
        #            compartmentResults[[j]] = compartmentResults[[j]] + (simResults[[i]]$compartments)[[j]]
        #        }
        #    }
        #}
        #biasResults = biasResults/length(simResults)
        #timeResult = timeResult/length(simResults)
        #iterationResult = iterationResult/length(simResults)
        #for (j in 1:4)
        #{
        #    compartmentResults[[j]] = compartmentResults[[j]]/length(simResults)
        #}

        #outData = list(biasResults=biasResults, timeResult=timeResult, iterationResult=iterationResult, throwAwayTpt=ThrowAwayTpt, compartmentResult=compartmentResults)
        #save(outData, file=paste("./sim2_results_", ThrowAwayTpt, ".tmp", sep=""))
    }
    print("Results obtained")
    stopCluster(cl)
    TRUE
}
