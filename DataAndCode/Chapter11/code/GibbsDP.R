for = (iter in 1:(MaxIter-1)) { 
    PostDirichletParms[iter,] = TempDirichletParms 
    ProjectedFailTime = rep(NA, nCensor) 
    for = (j in 1:nCensor) { 
#  Grab the remaining time for the censored observation. 
        StartIndx = findInterval(CensorTime[j], TimeGrid) 
        PotentialFailTimes = TimeGrid[(StartIndx+1):nGrid] 
#  Use multinomial sampling over the remaining time grid to 
#  generate a future failure time for the censored observation.
#  The multinomial probabilities are proportional to the 
#  parameters from vector TempDirichletParms corresponding 
#  to the part of the time grid to the right of the censored time.  
        ProjectedFailTime[j] =
            sample(PotentialFailTimes, 1, prob=Probs[iter, (StartIndx+1):nGrid]) 
    } 
#  Add the imputed failure times to the current Dirichlet  
#  distribution parameter vector. 
    MatchIndx = match(ProjectedFailTime, TimeGrid) 
    oMatchIndx = order(MatchIndx) 
    TempTable = table(MatchIndx) 
    for (i in 1:length(TempTable)) { 
        loc = as.integer(names(TempTable)[i]) 
        PostDirichletParms[iter, loc] = TempDirichletParms[loc] + TempTable[i] 
    } 
#  Generate a sample from the posterior Dirichlet distribution 
#  (with the imputed failure times at this iteration) 
#  via random gammas.
    Probs[iter+1,] = rdirichlet(1, PostDirichletParms[iter,]) 
}
