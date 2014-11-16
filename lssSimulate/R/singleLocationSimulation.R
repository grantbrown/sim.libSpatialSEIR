generateSingleLocData = function(seed, population, NYears, TptPerYear, ThrowAwayTpt)
{
    set.seed(seed)
    MaxTpt = NYears*TptPerYear
    N=population
    X = matrix(1, ncol = 1)
    Z = cbind(seq(1,NYears*TptPerYear), sin(seq(1,NYears*TptPerYear)/TptPerYear*2*pi))
    X_prs = cbind(1, 
                  sin((1:MaxTpt)/TptPerYear*2*pi), 
                  cos((6+1:MaxTpt)/TptPerYear*2*pi))

    trueBetaSEFixed = c(0.5)
    trueBetaSEVarying = c(0.002, .5)
    trueGamma = rep(0.0, MaxTpt)  

    eta_se = as.numeric((X %*% trueBetaSEFixed)) + (Z %*% trueBetaSEVarying)
    p_se = numeric(MaxTpt)
    gamma_ei = 2.3
    gamma_ir = 2.3

    p_ei = 1-exp(-gamma_ei)
    p_ir = 1-exp(-gamma_ir)

    trueBetaRS = c(-2.5, 1, -0.25) 
    eta_rs = X_prs %*% trueBetaRS
    p_rs = 1-exp(-exp(eta_rs))

    E0 = 0
    I0 = floor(0.001*N)
    R0 = floor(0.001*N) 
    S0 = N-E0-I0-R0

    S = matrix(0, nrow = 1, ncol = MaxTpt)
    E = matrix(0, nrow = 1, ncol = MaxTpt)
    I = matrix(0, nrow = 1, ncol = MaxTpt)
    R = matrix(0, nrow = 1, ncol = MaxTpt)
    S_star = matrix(0, nrow = 1, ncol = MaxTpt)
    E_star = matrix(0, nrow = 1, ncol = MaxTpt)
    I_star = matrix(0, nrow = 1, ncol = MaxTpt)
    R_star = matrix(0, nrow = 1, ncol = MaxTpt)

    for (i in 1:MaxTpt)
    {
        if (i == 1)
        {
            S[i] = S0
            E[i] = E0
            I[i] = I0
            R[i] = R0

            etaVal = -(I[i]/N)*exp(eta_se[i]) - trueGamma[i]
            p_se[i] = 1-exp(etaVal)

            S_star[i] = rbinom(1, R[i], p_rs[i])
            E_star[i] = rbinom(1, S[i], p_se[i])
            I_star[i] = rbinom(1, E[i], p_ei)
            R_star[i] = rbinom(1, I[i], p_ir)
        }
        else
        {
            S[i] = S[i-1] + S_star[i-1] - E_star[i-1]
            E[i] = E[i-1] + E_star[i-1] - I_star[i-1]
            I[i] = I[i-1] + I_star[i-1] - R_star[i-1]
            R[i] = R[i-1] + R_star[i-1] - S_star[i-1]

            etaVal = -(I[i]/N)*exp(eta_se[i]) - trueGamma[i]
            p_se[i] = 1-exp(etaVal)

            S_star[i] = rbinom(1, R[i], p_rs[i])
            E_star[i] = rbinom(1, S[i], p_se[i])
            I_star[i] = rbinom(1, E[i], p_ei)
            R_star[i] = rbinom(1, I[i], p_ir)
        }
    }

    if (ThrowAwayTpt != 0)
    {
        S0 = S[,ThrowAwayTpt+1]
        E0 = E[,ThrowAwayTpt+1]
        I0 = I[,ThrowAwayTpt+1]
        R0 = R[,ThrowAwayTpt+1]
    }

    S_star = S_star[,(ThrowAwayTpt + 1):ncol(S_star), drop = FALSE]
    E_star = E_star[,(ThrowAwayTpt + 1):ncol(E_star), drop = FALSE]
    I_star = I_star[,(ThrowAwayTpt + 1):ncol(I_star), drop = FALSE]
    R_star = R_star[,(ThrowAwayTpt + 1):ncol(R_star), drop = FALSE]

    S = S[,(ThrowAwayTpt + 1):ncol(S), drop = FALSE]
    E = E[,(ThrowAwayTpt + 1):ncol(E), drop = FALSE]
    I = I[,(ThrowAwayTpt + 1):ncol(I), drop = FALSE]
    R = R[,(ThrowAwayTpt + 1):ncol(R), drop = FALSE]

    Z = Z[(1+ThrowAwayTpt*nrow(S)):nrow(Z),]
    X_prs = X_prs[(ThrowAwayTpt + 1):nrow(X_prs),]

    # Transpose Everything to have TXP

    S0 = t(S0)
    E0 = t(E0)
    I0 = t(I0)
    R0 = t(R0)

    S_star = t(S_star)
    E_star = t(E_star)
    I_star = t(I_star)
    R_star = t(R_star)

    S = t(S)
    E = t(E)
    I = t(I)
    R = t(R)

    xDim = dim(X)
    zDim = dim(Z)
    xPrsDim = dim(X_prs)
    compMatDim = c(nrow(S), ncol(S))

    p_rs = p_rs[(ThrowAwayTpt +1):length(p_rs)]
    p_se = p_se[(ThrowAwayTpt +1):length(p_se)]
    eta_se = eta_se[(ThrowAwayTpt +1):length(eta_se)]
    beta = c(trueBetaSEFixed, trueBetaSEVarying) 
    betaPrs = trueBetaRS
    N = matrix(N, nrow = nrow(S), ncol = ncol(S))

    simResults = list("S_star"=S_star, "E_star" = E_star, "I_star" = I_star, "R_star" = R_star,
                "S"=S, "E"=E, "I"=I, "R"=R, "S0"=S0,"E0"=E0,"I0"=I0,"R0"=R0,"X"=X, "Z"=Z,"X_prs"=X_prs,
                "p_se"=p_se,"p_rs"=p_rs, "beta" = beta, "betaPrs"=betaPrs, "N"=N, "gamma_ei"=gamma_ei, 
                "gamma_ir"=gamma_ir)


}

buildSingleLocSimInstance = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]

    set.seed(seed)

    DataModel = buildDataModel(simResults$I_star, type = "overdispersion", params = c(10,0.1))

    ExposureModel = buildExposureModel(simResults$X, simResults$Z, 
                                       beta = c(2, rep(0, ((length(simResults$beta))-1))), betaPriorPrecision = 0.1)
    ReinfectionModel = buildReinfectionModel("SEIRS", X_prs = simResults$X_prs, 
                                             betaPrs = -c(4, rep(0,(length(simResults$betaPrs)-1))), 
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
    TransitionPriors = buildTransitionPriorsManually(2300,1000,2300,1000)

    I0 = max(simResults$I_star[1], simResults$I_star[2], 100)
    E0 = I0
    S0 = simResults$N[1] - I0 - E0
    InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                               S0 = S0,, I0 = I0, E0 = E0)
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,DistanceModel,TransitionPriors,
                         InitContainer,SamplingControl)

    res$setRandomSeed(seed+1)
    res$simulate(10000)
    res$compartmentSamplingMode = 16
    res$useDecorrelation = 100
    res$performHybridStep = 100

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
        localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
    }
}

checkConvergence = function(fileName1,fileName2,fileName3,maxVal=1.1,useUpper=FALSE,verbose=TRUE)
{
    dat1 = read.csv(fileName1)
    dat2 = read.csv(fileName2)
    dat3 = read.csv(fileName3)
    
    iteration = dat1$Iteration

    dat1 = as.mcmc(dat1[,1:(ncol(dat1) - 2)])
    dat2 = as.mcmc(dat2[,1:(ncol(dat2) - 2)])
    dat3 = as.mcmc(dat3[,1:(ncol(dat3) - 2)])
    mcl = mcmc.list(dat1,dat2,dat3)
    diag = gelman.diag(mcl, multivariate=FALSE)
    if (useUpper)
    {
        criterion = max(diag[[1]][,2])
        maxParam = rownames(diag[[1]])[which.max(diag[[1]][,2])]
    }
    else
    {
        criterion = max(diag[[1]][,1])
        maxParam = rownames(diag[[1]])[which.max(diag[[1]][,1])]
    }
    if (verbose)
    {
        cat(paste("Convergence Criterion: ", round(criterion, 2), 
                  ", param: ", maxParam, "\n", sep = ""))
    }
    criterion < maxVal
}



simulation1Kernel = function(cl, genSeed, fitSeeds, population, NYears, TptPerYear, ThrowAwayTpt)
{
    #simResults = generateSingleLocData = function(seed, population, NYears, TptPerYear, ThrowAwayTpt)
    simResults = generateSingleLocData(genSeed, population, NYears, TptPerYear, ThrowAwayTpt)

    fileNames = c("sim1_1.txt", "sim1_2.txt", "sim1_3.txt")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], simResults),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], simResults),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], simResults)
                      )

    trueVals = parLapply(cl, paramsList, buildSingleLocSimInstance)

    iterationParams = list(list(20000, 100, 0.15, 0.05, 0.1),
                           list(20000, 100, 0.15, 0.05, 0.1),
                           list(20000, 100, 0.15, 0.05, 0.1))
    conv = FALSE
    while (!conv)
    {
        cat("Not converged, adding iterations...\n")
        parLapply(cl, iterationParams, additionalIterations)
        conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3])
    }

    computeSim1Results(fileNames[1], fileNames[2], fileNames[3], simResults)
}

computeSim1Results = function(fileName1, fileName2, fileName3, trueData)
{
    contains = function(trueVal, variableName, summarylist)
    {
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

    dat1 = dat1[floor(nrow(dat1)/2):nrow(dat1),]
    dat2 = dat2[floor(nrow(dat2)/2):nrow(dat2),]
    dat3 = dat3[floor(nrow(dat3)/2):nrow(dat3),]


    iteration = dat1$Iteration
    time = dat1$Time
    

    dat1 = as.mcmc(dat1[,1:(ncol(dat1) - 2)])
    dat2 = as.mcmc(dat2[,1:(ncol(dat2) - 2)])
    dat3 = as.mcmc(dat3[,1:(ncol(dat3) - 2)])
    mcl = mcmc.list(dat1,dat2,dat3)

    mcl.summary = summary(mcl)

    trueBeta = trueData$beta
    trueBetaPrs = trueData$betaPrs
    trueGamma_ei = trueData$gamma_ei
    trueGamma_ir = trueData$gamma_ir

    betaNames = paste("BetaP_SE_", seq(0, (length(trueBeta)-1)), sep = "")
    betaPrsNames = paste("BetaP_RS_", seq(0, (length(trueBetaPrs)-1)), sep = "")

    allCompareVariableNames = c(betaNames, betaPrsNames, "gamma_ei", "gamma_ir")
    allCompareVariableValues = c(trueBeta, trueBetaPrs, trueGamma_ei, trueGamma_ir)

    outMatrix = matrix(0, nrow = length(allCompareVariableNames), ncol=3)
    colnames(outMatrix) = c("CI_Contains", "Bias", "BiasPct")
    rownames(outMatrix) = allCompareVariableNames
    for (varNum in 1:length(allCompareVariableNames))
    {
        varName = allCompareVariableNames[varNum]
        varVal = allCompareVariableValues[varNum]
        bias = getMean(varName, mcl.summary) - varVal
        outMatrix[varNum,] = c(contains(varVal, varName, mcl.summary), bias, bias/(abs(varVal))*100)
    }
    return(list("data"=outMatrix,"iterations"=max(iteration), "time"=max(time)))
}

runSimulation1 = function(cellIterations = 50, ThrowAwayTpts=c(0, 6, 12, 18, 24),
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
            simulation1Kernel(cl, genSeed, fitSeeds, 100000, 3, 12, ThrowAwayTpt)
        }
        simResults = lapply(genSeed + seq(1, cellIterations), f)
        biasResults = simResults[[1]]$data
        timeResult = simResults[[1]]$time
        iterationResult = simResults[[1]]$iterations
        for (i in 2:length(simResults))
        {
            biasResults = biasResults + simResults[[i]]$data
            timeResult = timeResult + simResults[[i]]$time
            iterationResult = iterationResult + simResults[[i]]$iterations
        }
        biasResults = biasResults/length(simResults)
        timeResult = timeResult/length(simResults)
        iterationResult = iterationResult/length(simResults)
        outData = cbind(biasResults, timeResult, iterationResult, ThrowAwayTpt)
        colnames(outData) = c("Coverage", "AvgBias", "AvgBiasPct", "AvgTime","AvgIterations", "ThrowAwayTpt")
        write.csv(outData, file=paste("./sim1_results_", ThrowAwayTpt, ".txt",sep=""))
    }
    print("Results obtained")
    stopCluster(cl)
    TRUE
}
