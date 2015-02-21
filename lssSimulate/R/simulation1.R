#generateSingleLocData = function(seed, population, NYears, TptPerYear, ThrowAwayTpt)
generateData = function(seed, 
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
{
    set.seed(seed)
    

    if (!is.na(Z) && nTpt > nrow(Z)/nrow(X))
    {
        stop(paste("Not enough information provided for ", nTpt, " time points.\n", sep = ""))
    }

    # Make a big matrix out of X and Z
    if (all(is.na(Z)))
    {
        bigX_SE = cbind(X[rep(c(1:nrow(X)), each=nTpt),])
    }
    else
    {

        Z = Z[rep(1:nTpt, nrow(X)),]
        bigX_SE = cbind(X[rep(c(1:nrow(X)), nTpt),], Z)
    }

    # Calculate temporal offsets. 
    uncumulate = function(x)
    {
        out = c(x[2:length(x)]-x[1:(length(x)-1)])
        ifelse(out >= 0, out, 0)
    }
    if (length(timeIndex) - nrow(X_RS) != 1)
    {
        stop("ERROR: timeIndex should be one longer than the number of time points in the study.")
    }
    offsets = uncumulate(timeIndex)[1:nTpt]

    # Calculate linear predictor for the intensity process. 
    eta_SE = exp(matrix((bigX_SE %*% beta_SE), ncol = nrow(X)))

    # Calculate linear predictor for reinfection process.
    X_RS = X_RS[1:nTpt,,drop=FALSE]
    eta_RS = exp(X_RS %*% beta_RS)

    # Calculate p_EI and p_IR
    p_EI = 1-exp(-gamma_ei*offsets) 
    p_IR = 1-exp(-gamma_ir*offsets)


    # p_RS
    p_RS = 1-exp(-offsets*(eta_RS))

    # Allocate compartments
    S_star = E_star = I_star = R_star = S = E = I = R = matrix(0, nrow = nTpt, ncol = nrow(X))

    # Declare N
    N = matrix(S0+E0+I0+R0, ncol = nrow(X), nrow = length(offsets), byrow=TRUE)

    # Run Simulation 
    offsetMatrix = matrix(offsets, nrow = length(offsets), ncol = nrow(X))
    S[1,] = S0
    E[1,] = E0
    I[1,] = I0
    R[1,] = R0

    p_SE = eta_SE
    for (i in 1:(nTpt))
    {
        # Calculate p_SE
        if (nrow(X) > 1)
        {
            # Case 1: spatial
            p_SE[i,] = p_SE[i,]*(I[i,]/N[i,]) 
            for (j in 1:length(distMatList))
            {
               p_SE[i,] = p_SE[i,] + rho[[j]]*distMatList[[j]] %*% (eta_SE[i,]*(I[i,]/N[i,]))
            }
            p_SE[i,] = 1-exp(-offsetMatrix[i,]*p_SE[i,])
        }
        else
        {
            # Case 2: single location
            p_SE[i,] = (p_SE[i,]*(I[i,]/N[i,]))
            p_SE[i,] = 1-exp(-offsetMatrix[i,]*p_SE[i,])
        }

        S_star[i,] = rbinom(rep(1, ncol(S)), R[i,], rep(p_RS[i], ncol(S)))
        E_star[i,] = rbinom(rep(1, ncol(S)), S[i,], p_SE[i,])
        I_star[i,] = rbinom(rep(1, ncol(S)), E[i,], rep(p_EI[i], ncol(S)))
        R_star[i,] = rbinom(rep(1, ncol(S)), I[i,], rep(p_IR[i], ncol(S)))

        if (i < nTpt)
        {
            S[i+1,] = S[i,] + S_star[i,] - E_star[i,] 
            E[i+1,] = E[i,] + E_star[i,] - I_star[i,]
            I[i+1,] = I[i,] + I_star[i,] - R_star[i,]
            R[i+1,] = R[i,] + R_star[i,] - S_star[i,]
        }
    }

    return(list(S_star=S_star,
                E_star=E_star,
                I_star=I_star,
                R_star=R_star,
                S=S,
                I=I,
                E=E,
                R=R,
                S0=S0,
                E0=E0,
                I0=I0,
                R0=R0,
                N =N, 
                X=X,
                Z=Z,
                X_RS=X_RS,
                beta_SE=beta_SE,
                beta_RS=beta_RS,
                gamma_ei=gamma_ei,
                gamma_ir=gamma_ir,
                distMatList=distMatList,
                p_SE=p_SE,
                p_EI=p_EI,
                p_IR=p_IR,
                p_RS=p_RS,
                effectiveTransitionSampleSize=effectiveTransitionSampleSize
                ))
}

generateSingleLocData = function(seed, 
                                 population, 
                                 NYears, 
                                 TptPerYear, 
                                 ThrowAwayTpt,
                                 beta_SE=c(0.5, 0.002, 0.5),
                                 beta_RS=c(-2.5))
{
    set.seed(seed)
    E0 = 0  
    I0 = floor(0.001*population)
    R0 = I0
    S0 = population - R0 - E0 - I0

    maxTpt = NYears*TptPerYear
    nTpt = maxTpt - ThrowAwayTpt
    
    X = matrix(1)
    #Z = matrix(sin(2*pi*((1:maxTpt)/TptPerYear)), ncol = 1) 
    Z = cbind(seq(1, maxTpt), sin(seq(1, maxTpt)/TptPerYear*2*pi))

    X_RS = matrix(dnorm(5*sin(seq(1, maxTpt)/TptPerYear*pi)), nrow = maxTpt)

    distMatList = -1   
    rho = -1
    gamma_ei = 2.3
    gamma_ir = 2.3
    
    effectiveTransitionSampleSize = 1000
    timeIndex = 0:maxTpt

    generateData(seed + 1, 
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

}


buildSingleLocSimInstance = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]

    set.seed(seed)

    #DataModel = buildDataModel(simResults$I_star, type = "overdispersion", params = c(10000,10000))
    DataModel = buildDataModel(simResults$I_star, type = "identity")


    priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
    ExposureModel = buildExposureModel_depricated(simResults$X, simResults$Z, 
                                       beta = c(priorBetaIntercept, rep(0, ((length(simResults$beta_SE))-1))), 
                                       betaPriorMean = rep(0, length(simResults$beta_SE)), 
                                       betaPriorPrecision = 0.01)
    ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = simResults$X_RS, 
                                             betaPrs = -c(0, rep(0,(length(simResults$beta_RS)-1))), 
                                             priorMean = 0,
                                             priorPrecision = 0.1)
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
    TransitionPriors$summary()

    I0 = max(simResults$I_star[1], simResults$I_star[2], 100)
    E0 = I0
    S0 = simResults$N[1] - I0 - E0
    InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                               S0 = S0,I0 = I0, E0 = E0, reinfection=TRUE,dataType="I_star")
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
    for (i in 1:(200))
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

additionalIterations = function(params)
{
    
    N = params[[1]]
    batchSize = params[[2]]
    targetRatio = params[[3]]
    targetWidth=params[[4]]
    proportionChange = params[[5]]
    for (i in 1:(N/batchSize))
    {
        localModelObject$simulate(batchSize)
        #localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
    }
}

simulation1Kernel = function(cl, genSeed, fitSeeds, population, NYears, TptPerYear, ThrowAwayTpt)
{
    #TODO: Vary starting linear predictor parameters on each iteration 
    simResults = generateSingleLocData(genSeed, population, NYears, TptPerYear, ThrowAwayTpt)

    fileNames = c("sim1_1.txt", "sim1_2.txt", "sim1_3.txt")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults)
                      )

    trueVals = parLapply(cl, paramsList, buildSingleLocSimInstance)

    iterationParams = list(list(200000, 1000, 0.2, 0.05, 0.1),
                           list(200000, 1000, 0.2, 0.05, 0.1),
                           list(200000, 1000, 0.2, 0.05, 0.1))
    conv = FALSE
    while (!conv)
    {
        cat("Not converged, adding iterations...\n")
        parLapply(cl, iterationParams, additionalIterations)
        conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], maxVal = 1.02)
        if (is.na(conv)){conv = FALSE;}
    }

    dat1 = read.csv(fileNames[1])
    dat2 = read.csv(fileNames[2])
    dat3 = read.csv(fileNames[3])
    list(dat1=dat1,dat2=dat2,dat3=dat3,simResults=simResults)
    #computeSim1Results(fileNames[1], fileNames[2], fileNames[3], simResults)
}

computeSim1Results = function(fileName1, fileName2, fileName3, trueData)
{
    contains = function(trueVal, variableName, summarylist)
    {
        cat(paste("Checking for: ", variableName, "\n", sep = ""))
        sl = summarylist[[2]]
        variableIdx = which(rownames(sl) == variableName)
        bounds = sl[variableIdx,][c(1,5)]
        trueVal < bounds[2] && trueVal > bounds[1]
    }

    getMean = function(variableName, summaryList)
    {
        sl = summaryList[[1]]
        variableIdx = which(rownames(sl) == variableName)
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
    
    dat1 = as.mcmc(dat1_full[,1:(ncol(dat1_full) - 2)])
    dat2 = as.mcmc(dat2_full[,1:(ncol(dat2_full) - 2)])
    dat3 = as.mcmc(dat3_full[,1:(ncol(dat3_full) - 2)])
    mcl = mcmc.list(dat1,dat2,dat3)

    mcl.summary = summary(mcl)

    # Calculate parameter summaries
    trueBeta = trueData$beta_SE
    trueBetaPrs = trueData$beta_RS
    trueGamma_ei = trueData$gamma_ei
    trueGamma_ir = trueData$gamma_ir

    betaNames = paste("BetaP_SE_", seq(0, (length(trueBeta)-1)), sep = "")
    betaPrsNames = paste("BetaP_RS_", seq(0, (length(trueBetaPrs)-1)), sep = "")

    allCompareVariableNames = c(betaNames, betaPrsNames, "gamma_ei", "gamma_ir")
    allCompareVariableValues = c(trueBeta, trueBetaPrs, trueGamma_ei, trueGamma_ir)

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
    tpts = nrow(trueData$X_RS) 
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

runSimulation1 = function(cellIterations = 50, ThrowAwayTpts=c(0,6,12,24),
                          genSeed=123123, fitSeeds=c(812123,12301,5923))
{                     
    cl = makeCluster(3, outfile = "err.txt")
    print("Cluster Created")
    clusterExport(cl, c("buildSingleLocSimInstance"))
    print("Variables Exported.") 
 
    for (ThrowAwayTpt in ThrowAwayTpts)
    {
        f = function(genSeed)
        {
            simulation1Kernel(cl, genSeed, fitSeeds + genSeed, 1000, 3, 12, ThrowAwayTpt)
        }
        itrSeeds = genSeed + seq(1, cellIterations)
        i = 1
        for (itrSeed in itrSeeds){
            result = f(itrSeed)
            save(result, file = paste("./sim1_results_", ThrowAwayTpt, "_", i, ".Rda.bz2", sep = ""), compress="bzip2")
            i = i+1
        }
    }
    print("Results obtained")
    stopCluster(cl)
    TRUE
}
