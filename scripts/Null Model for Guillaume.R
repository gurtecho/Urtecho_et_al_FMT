
set.seed(123)
ParallelSim9 <- function(BINARYMAT,NUMBERRANDOMCOMM,CORES){
  
  # initialize the parallel infrastructure
  library(doParallel)
  registerDoParallel(CORES)
  numberiterations <- NUMBERRANDOMCOMM
  
  # do the simulation in parallel, with X cores each performing an independent simulation of Y shuffles
  # the .combine=rbind means the resulting data will be rbinded together
  simulatedresults <- foreach(i = icount(numberiterations), .combine=rbind,.packages=c('reshape2',"data.table")) %dopar% {
    # get the initial otu matrix
    testmat2 <- BINARYMAT
    
    # perform the sim9 shuffling X times
    for(j in 1:25000){
      testmat2 <- EcoSimR::sim9_single(testmat2)
    }
    
    # melt it and only save instances that contain an actual count, label the iteration
    testres <- as.data.table(reshape2::melt(testmat2))[value > 0,]
    testres$Iteration <- i
    print(paste("Finished Iteration",i))
    return(testres)
  }
  names(simulatedresults) <- c("Particle","OTU","Count","Iteration")
  
  return(simulatedresults)
}

# perform the shuffling of otus in parallel
ParallelLabelShuffle <- function(BINARYMAT,NUMBERRANDOMCOMM,CORES){
  
  # initialize the parallel infrastructure
  library(doParallel)
  registerDoParallel(CORES)
  numberiterations <- NUMBERRANDOMCOMM
  
  # shuffle the otu labels in parallel, and bind them together
  shuffledresults <- foreach(i = icount(numberiterations), .combine=rbind,.packages=c('reshape2',"data.table")) %dopar% {
    shuffledmatrix <- BINARYMAT
    
    print("Shuffling the OTU labels")
    colnames(shuffledmatrix) <- sample(colnames(shuffledmatrix))
    
    testres <- as.data.table(reshape2::melt(shuffledmatrix))[value > 0,]
    testres$Iteration <- i
    print(paste("Finished Iteration",i))
    return(testres)
  }
  names(shuffledresults) <- c("Particle","OTU","Count","Iteration")
  
  return(shuffledresults)
}


DistributionDistanceComparison <- function(ORIGINALMATRIX,SIMRESULTS,DISTANCEMETRIC,DONOR,RANDOMIZATION){
  
  # calculate the otu associations in each simulatied community in parallel
  distanceiteration <- foreach(i = unique(SIMRESULTS$Iteration), .combine=rbind,.packages=c('reshape2',"parallelDist","data.table")) %dopar% {
    
    # cast data and put into appropriate form
    casted <- dcast.data.table(SIMRESULTS[Iteration == i],OTU~Particle,value.var = "Count")
    setnafill(casted,fill=0, cols=2:ncol(casted))
    castmat <- as.matrix(casted[,2:ncol(casted)])
    row.names(castmat) <- casted$OTU
    
    # perform the distance calculation and remove duplicate values
    distres <- as.matrix(parallelDist::parallelDist(castmat,DISTANCEMETRIC))
    
    distres[upper.tri(distres)] <- NA
    diag(distres) <- NA
    
    # make it a datatable and give it the iteration number
    distmelt <- as.data.table(reshape2::melt(distres,na.rm=T))
    distmelt$Iteration <- i
    
    return(distmelt)
  }
  
  # calculate distances for the original community
  distres <- as.matrix(parallelDist::parallelDist(t(ORIGINALMATRIX),DISTANCEMETRIC))
  
  # melt the original distance and merge it with the simulated communities
  distmelt <- as.data.table(reshape2::melt(distres,na.rm=T))
  names(distmelt)[3] <- "Original"
  combinedmatrix <- merge(distanceiteration,distmelt,by=c("Var1","Var2"))
  

  
  # calculated means/sds for each otu pair, then calculate zscores and p values
  combinedsummary <- combinedmatrix[,.(meandist=mean(value),
                                       distsd=sd(value),
                                       median=median(value),
                                       Original=unique(Original)),
                                    by=.(Var1,Var2)]
  
  combinedsummary[,Zscore := ((meandist-Original)/distsd)]
  combinedsummary[,zpvalue := 2*pnorm(abs(Zscore), lower.tail=FALSE)]
  
  combinedsummary$zpadjusted <- signif(p.adjust(combinedsummary$zpvalue,"fdr"), 3)
  


  otuorder <- names(colSums(ORIGINALMATRIX))[order(colSums(ORIGINALMATRIX),decreasing = T)]
  combinedsummary$Var1 <- factor(combinedsummary$Var1,levels=otuorder)
  combinedsummary$Var2 <- factor(combinedsummary$Var2,levels=otuorder)
  
  combinedsummary2 <- combinedsummary
  names(combinedsummary2)[1:2] <- c("Var2","Var1") 
  
  combinedplot <- rbind(combinedsummary,combinedsummary2)
  
  
  # get the significant pairs
  combinedplot$Significant <- combinedplot$zpadjusted < .05
  
  
  # make sure there is no weirdness with the parallel stuff
  closeAllConnections() 
  
  return(combinedsummary)
}





