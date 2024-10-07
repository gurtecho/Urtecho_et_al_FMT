Compute_Spearman_Correlation_And_Melt <- function(OTUTABLE,DATANAME,READCOUNT){
  print("Calculating")
  Correlations <- Hmisc::rcorr(OTUTABLE,type="spearman")
  #melt into pairwise data
  print("Melting")
  MeltedCorrelations <- reshape2::melt(Correlations$r)
  names(MeltedCorrelations) <- c("OTU1","OTU2","Correlation")
  MeltedSignificance <- reshape2::melt(Correlations$P)
  names(MeltedSignificance) <- c("OTU1","OTU2","Significance")
  AllMelted <- merge(MeltedCorrelations,MeltedSignificance)
  #label DF with
  print("Labeling Sample")
  AllMelted$Sample <- paste(DATANAME)
  AllMelted$RarefactionThreshold <- paste(READCOUNT)
  #adjust p-values
  AllMelted$SpearmanAdjustedP <- p.adjust(AllMelted$Significance,"bonferroni")
  AllMelted <- subset(AllMelted,Correlation < 1)
  #Returning
  print("Done")
  return(AllMelted)
}