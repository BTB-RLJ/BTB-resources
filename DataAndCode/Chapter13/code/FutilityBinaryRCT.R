for (iSim in 1:nSim) { 
    for (iGroup in 1:nGroup) { 
#		Generate data for the new cohort randomized to each of the two treatments 
        Control[iGroup, 2] = rbinom(1, n0PerGroup[iGroup], p0) 
        Trt[iGroup, 2] = rbinom(1, n1PerGroup[iGroup], p1) 
#		Parameters for the posterior beta distribution 
        TotalControlResponses = sum(Control[1:iGroup, 2]) 
        TotalControlNoResponse = sum(Control[1:iGroup, 1]) - TotalControlResponses 
        aControlPost = aControl + TotalControlRespond  
        bControlPost = bControl + TotalControlNoRespond 
        TotalTrtResponses = sum(Trt[1:iGroup, 2]) 
        TotalTrtNoResponse = sum(Trt[1:iGroup, 1]) - TotalTrtResponses 
        aTrtPost = aTrt + TotalTrtResponses; bTrtPost = bTrt + TotalTrtNoResponses 
#		Generate remaining obs from predictive beta-binomial dist'n 
        n0Remaining = n0 - n0Entered; n1Remaining = n1 - n1Entered 
        Pred.pValue = rep(NA, nPred) 
        if (iGroup < nGroup) { 
            x0Pred = RandomBetaBinomial(nPred, n0Remaining, aControlPost, bControlPost) 
            PredTotalControlResponses = TotalControlResponses + x0Pred 
            PredTotalControlNoResponse = TotalControlNoResponse + 
                (n0Remaining - x0Pred) 
            x1Pred = RandomBetaBinomial(nPred, n1Remaining, aTrtPost, bTrtPost) 
            PredTotalTrtResponses = TotalTrtResponses + x1Pred 
            PredTotalTrtNoResponse = TotalTrtNoResponse + (n1Remaining - x1Pred) 
#		Generate future obs'ns and combine with current obs'ns for final test.  
            Pred.pValue = rep(NA, nPred) 
            for (iPred in 1:nPred) { 
                Tab = matrix(c(PredTotalTrtRespond[iPred], PredTotalTrtNoRespond[iPred], 
                                        PredTotalControlRespond[iPred], PredTotalControlNoRespond[iPred]),  
                                        ncol=2, byrow=T) 
                 Pred.pValue[iPred] = fisher.test(Tab, alternative="greater")$p.value  
            } 
            PredProbReject[iSim, iGroup] =  sum(Pred.pValue < pValueToReject) / length(Pred.pValue) 
            if (PredProbReject[iSim, iGroup] < CutOff) { 
                Stop[iSim] = 1.0 
                break 
            } 
        } 
    } 
    SampleSize[iSim,] = c(n0temp, n1temp) 
    if ((iGroup == nGroup) & (Stop[iSim]==0)) { 
        TabFinal = matrix(c(TotalTrtRespond, TotalTrtNoRespond, 
                                        TotalControlRespond, TotalControlNoRespond),
					ncol=2, byrow=T)
        pValue[iSim] = fisher.test(TabFinal, alternative="g")$p.value 
        SampleSize[iSim,] = c(n0temp, n1temp) 
    }
}
