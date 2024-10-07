library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }   
  }
}


data.merge <- read.table("./assembly.COG.txt", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.merge) <- c("assembly", "CDS", "COG", "evalue")
data.geneNum.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly) %>% 
                            summarise(CDScount = length(unique(CDS))),
                            stringsAsFactors = F)
rownames(data.geneNum.merge.stat) <- data.geneNum.merge.stat$assembly
data.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly, COG) %>% 
                            summarise(count = n()),
                            stringsAsFactors = F)
data.merge.stat.dcast <- dcast(data.merge.stat, assembly ~ COG, fill = 0)
rownames(data.merge.stat.dcast) <- data.merge.stat.dcast$assembly
data.merge.stat.dcast <- data.merge.stat.dcast[,-1]
data.merge.stat.dcast.normalize <- sweep(data.merge.stat.dcast, 1, data.geneNum.merge.stat[rownames(data.merge.stat.dcast),"CDScount"], "/")


data.supplier.info <- read.table("./assembly.supplier.tsv", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.supplier.info) <- c("assembly", "supplier")
rownames(data.supplier.info) <- data.supplier.info$assembly
data.tax.info <- read.table("./assembly.taxonomy_parsed.tsv", stringsAsFactors = F, header = T, sep = "\t")
rownames(data.tax.info) <- data.tax.info$assembly
data.abundance.info <- read.table("./assembly.abundanceMean.stat", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.abundance.info) <- c("assembly", "abundance")
rownames(data.abundance.info) <- data.abundance.info$assembly
assembly.metadata.df <- data.frame(assembly = rownames(data.supplier.info),
								supplier = data.supplier.info[rownames(data.supplier.info), "supplier"],
								abundance = data.abundance.info[rownames(data.supplier.info), "abundance"],
								order = data.tax.info[rownames(data.supplier.info),"order"],
								family = data.tax.info[rownames(data.supplier.info),"family"],
								row.names = rownames(data.supplier.info),
								stringsAsFactors = F)

### MDS for all Bacteroidales ###
sample.list <- assembly.metadata.df[which(assembly.metadata.df$order == "Bacteroidales"),"assembly"]
data.matrix.used <- data.merge.stat.dcast.normalize[sample.list,]

data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            stringsAsFactors = F)
pdf("./figures/01.COG_analysis.MDS_Bacteroidales.pdf", height = 5, width = 7)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = family, shape = supplier)) +
    geom_point(alpha = 0.9, size = 1.2) +
    scale_shape_manual(values = c(16, 16, 8, 16)) +
    scale_color_manual(values = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_classic() +
    labs(x = "MDS1 (COG categories)", y = "MDS2 (COG categories)", title = "Order: Bacteroidales")
dev.off()
##################################


### MDS for all Muribaculaceae ###
sample.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"]
data.matrix.used <- data.merge.stat.dcast.normalize[sample.list,]

data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            stringsAsFactors = F)
pdf("./figures/01.COG_analysis.MDS_Muribaculaceae.pdf", height = 5, width = 7)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = supplier)) +
    geom_point(alpha = 0.9, size = 1.2) +
    scale_shape_manual(values = c(16, 16, 8, 16)) +
    scale_color_manual(values = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_classic() +
    labs(x = "MDS1 (COG categories)", y = "MDS2 (COG categories)", title = "Family: Muribaculaceae")
dev.off()
##################################


### MDS for all assembly ###
sample.list <- assembly.metadata.df[,"assembly"]
data.matrix.used <- data.merge.stat.dcast.normalize[sample.list,]

data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            stringsAsFactors = F)
pdf("./figures/01.COG_analysis.MDS_all.pdf", height = 5, width = 10)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = order, shape = supplier)) +
    geom_point(alpha = 0.9, size = 1.2) +
    scale_shape_manual(values = c(16, 16, 8, 16)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_classic() +
    labs(x = "MDS1 (COG categories)", y = "MDS2 (COG categories)", title = "All assembly")
dev.off()

pdf("./figures/01.COG_analysis.MDS_all_colorChange.pdf", height = 5, width = 8)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = supplier)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_shape_manual(values = c(16, 16, 8, 16)) +
	scale_color_manual(values = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_classic() +
    labs(x = "MDS1 (COG categories)", y = "MDS2 (COG categories)", title = "All assembly")
dev.off()
#############################


sample.null.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae" & assembly.metadata.df$supplier != "AZC7E"),"assembly"]
sample.exp.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae" & assembly.metadata.df$supplier == "AZC7E"),"assembly"]

#sample.null.list <- assembly.metadata.df[which(assembly.metadata.df$family != "Muribaculaceae" & assembly.metadata.df$order == "Bacteroidales"),"assembly"]
#sample.exp.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"]

data.dataframe.used <- data.merge

data.COG.notes <- read.table("./archive/cognames2003-2014.tab", stringsAsFactors = F, header = F, sep = "\t", quote = "") 
colnames(data.COG.notes) <- c("COG", "term", "note")
rownames(data.COG.notes) <- data.COG.notes$COG
COG.unique.list <- as.character(unique(data.merge$COG))

COG.info <- data.frame(row.names = COG.unique.list,
					COG = COG.unique.list,
					annotation = data.COG.notes[COG.unique.list, "note"])

### individual COG analysis ###
COG.enrich.info <- data.frame(COG_term = rep(" ", nrow(COG.info)),
                            COG_notes = rep(" ", nrow(COG.info)),
                            null_species = rep(0, nrow(COG.info)),
                            exp_species = rep(0, nrow(COG.info)),
                            null_COG = rep(0, nrow(COG.info)),
                            exp_COG = rep(0, nrow(COG.info)),
                            null_gene = rep(0, nrow(COG.info)),
                            exp_gene = rep(0, nrow(COG.info)),
                            null_COG_gene = rep(0, nrow(COG.info)),
                            exp_COG_gene = rep(0, nrow(COG.info)),
	                      	pvalue_species = rep(0, nrow(COG.info)),
                            odds_species = rep(0, nrow(COG.info)),
							pvalue_gene = rep(0, nrow(COG.info)),
                            odds_gene = rep(0, nrow(COG.info)),
                            stringsAsFactors = F)

for (i in 1:nrow(COG.info)){
    tmp.COG.term <- as.character(COG.info[i, "COG"])
    tmp.COG.note <- as.character(COG.info[i, "annotation"])
    tmp.df <- data.dataframe.used[which(data.dataframe.used$COG == tmp.COG.term),]
    totalSample <- unique(tmp.df$assembly)
    null_species = length(sample.null.list)
    exp_species = length(sample.exp.list)
    null_COG = length(totalSample[which(totalSample %in% sample.null.list)])
    exp_COG = length(totalSample[which(totalSample %in% sample.exp.list)])
	null_gene = sum(data.geneNum.merge.stat[sample.null.list, "CDScount"])
	exp_gene = sum(data.geneNum.merge.stat[sample.exp.list, "CDScount"])
	null_COG_gene = length(unique(tmp.df[which(tmp.df$assembly %in% sample.null.list),"CDS"]))
	exp_COG_gene = length(unique(tmp.df[which(tmp.df$assembly %in% sample.exp.list),"CDS"]))
    testor_species <- rbind(c(exp_COG, exp_species - exp_COG), c(null_COG, null_species - null_COG))
    model_species <- fisher.test(testor_species, alternative="two.sided")
	testor_gene <- rbind(c(exp_COG_gene, exp_gene - exp_COG_gene), c(null_COG_gene, null_gene - null_COG_gene))
    model_gene <- fisher.test(testor_gene, alternative="two.sided")
    COG.enrich.info[i,] <- c(tmp.COG.term, tmp.COG.note, null_species, exp_species, null_COG, exp_COG, 
							null_gene, exp_gene, null_COG_gene, exp_COG_gene,
							model_species$p.value, model_species$estimate,
							model_gene$p.value, model_gene$estimate)
}
COG.enrich.info$pvalue_species <- as.numeric(COG.enrich.info$pvalue_species)
COG.enrich.info$odds_species <- as.numeric(COG.enrich.info$odds_species)
COG.enrich.info$pvalue_gene <- as.numeric(COG.enrich.info$pvalue_gene)
COG.enrich.info$odds_gene <- as.numeric(COG.enrich.info$odds_gene)
COG.enrich.info$null_COG <- as.numeric(COG.enrich.info$null_COG)
COG.enrich.info$exp_COG <- as.numeric(COG.enrich.info$exp_COG)
COG.enrich.info$null_species <- as.numeric(COG.enrich.info$null_species)
COG.enrich.info$exp_species <- as.numeric(COG.enrich.info$exp_species)
COG.enrich.info$null_COG_gene <- as.numeric(COG.enrich.info$null_COG_gene)
COG.enrich.info$exp_COG_gene <- as.numeric(COG.enrich.info$exp_COG_gene)
COG.enrich.info$null_gene <- as.numeric(COG.enrich.info$null_gene)
COG.enrich.info$exp_gene <- as.numeric(COG.enrich.info$exp_gene)
COG.enrich.info$padj_species <- p.adjust(COG.enrich.info$pvalue_species, method = "BH")
COG.enrich.info$padj_gene <- p.adjust(COG.enrich.info$pvalue_gene, method = "BH")

Pvalue_threshold <- 0.05
COG.enrich.info$sig_gene <- "Not Sig"
COG.enrich.info[which(COG.enrich.info$padj_gene < Pvalue_threshold), "sig_gene"] <- "Sig"
COG.enrich.info$sig_species <- "Not Sig"
COG.enrich.info[which(COG.enrich.info$padj_species < Pvalue_threshold), "sig_species"] <- "Sig"

plotVector <- vector("list", length = 2)

p1 <- ggplot(COG.enrich.info, aes(x = log2(odds_gene), y = -log10(padj_gene))) +
	geom_point(aes(color = sig_gene), alpha = 0.3, size = 1) +
	geom_text_repel(data = COG.enrich.info[which(COG.enrich.info$sig_gene == "Sig"),], aes(label = paste(COG_term, ": ", COG_notes, sep = "")), size = 1.5, color = "#a50f15") +
	scale_color_manual(values=c("black", "red")) +
	labs(title = "COG enrichment on gene (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
	theme(legend.position="none") +
	xlim(-5.8, 6.2) + ylim(0, 5) + 
	theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(COG.enrich.info, aes(x = log2(odds_species), y = -log10(padj_species))) +
	geom_point(aes(color = sig_species), alpha = 0.3, size = 1) +
#	geom_text_repel(data = COG.enrich.info[which(COG.enrich.info$sig_species == "Sig"),], aes(label = paste(COG_term, ": ", COG_notes, sep = "")), size = 1.5, color = "#a50f15") +
	scale_color_manual(values=c("black", "red")) +
	labs(title = "COG enrichment on assembly (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
	theme(legend.position="none") +
	xlim(-5.8, 6.2) + ylim(0, 5) + 
	theme(plot.title = element_text(hjust = 0.5))
plotVector[[1]] <- p1
plotVector[[2]] <- p2

pdf("./figures/01.COG_analysis.individualTest.pdf", width = 10, height = 4) 
multiplot(plotlist = plotVector, cols = 2)
dev.off()


##############################



data.dataframe.used <- data.merge

data.COG.notes <- read.table("./archive/cognames2003-2014.tab", stringsAsFactors = F, header = F, sep = "\t", quote = "")
colnames(data.COG.notes) <- c("COG", "term", "note")
rownames(data.COG.notes) <- data.COG.notes$COG
COG.unique.list <- as.character(unique(data.merge$COG))

COG.info <- data.frame(row.names = COG.unique.list,
                    COG = COG.unique.list,
                    annotation = data.COG.notes[COG.unique.list, "note"])


assembly.cluster <- read.table("./assembly.CAZy_GH.cluster.tsv", stringsAsFactors = F, sep = "\t", header = F)
colnames(assembly.cluster) <- c("assembly", "cluster")
rownames(assembly.cluster) <- assembly.cluster$assembly

data.cluster.meta <- data.frame(row.names = assembly.cluster$assembly,
                                assembly = assembly.cluster$assembly,
                                cluster = assembly.cluster$cluster,
                                supplier = assembly.metadata.df[assembly.cluster$assembly,"supplier"],
                                abundance = assembly.metadata.df[assembly.cluster$assembly,"abundance"],
                                stringsAsFactors = F)

cluster.list <- c("CAZy_GH_1", "CAZy_GH_2", "CAZy_GH_3")
plotVector <- vector("list", length = 6)
for (j in 1:3){
    tmpCluster <- cluster.list[j]
    sample.null.list <- data.cluster.meta[data.cluster.meta$cluster != tmpCluster, "assembly"]
    sample.exp.list <- data.cluster.meta[data.cluster.meta$cluster == tmpCluster, "assembly"]
	COG.enrich.info <- data.frame(COG_term = rep(" ", nrow(COG.info)),
                            COG_notes = rep(" ", nrow(COG.info)),
                            null_species = rep(0, nrow(COG.info)),
                            exp_species = rep(0, nrow(COG.info)),
                            null_COG = rep(0, nrow(COG.info)),
                            exp_COG = rep(0, nrow(COG.info)),
                            null_gene = rep(0, nrow(COG.info)),
                            exp_gene = rep(0, nrow(COG.info)),
                            null_COG_gene = rep(0, nrow(COG.info)),
                            exp_COG_gene = rep(0, nrow(COG.info)),
	                      	pvalue_species = rep(0, nrow(COG.info)),
                            odds_species = rep(0, nrow(COG.info)),
							pvalue_gene = rep(0, nrow(COG.info)),
                            odds_gene = rep(0, nrow(COG.info)),
                            stringsAsFactors = F)
	for (i in 1:nrow(COG.info)){
    	tmp.COG.term <- as.character(COG.info[i, "COG"])
    	tmp.COG.note <- as.character(COG.info[i, "annotation"])
    	tmp.df <- data.dataframe.used[which(data.dataframe.used$COG == tmp.COG.term),]
    	totalSample <- unique(tmp.df$assembly)
    	null_species = length(sample.null.list)
    	exp_species = length(sample.exp.list)
    	null_COG = length(totalSample[which(totalSample %in% sample.null.list)])
    	exp_COG = length(totalSample[which(totalSample %in% sample.exp.list)])
		null_gene = sum(data.geneNum.merge.stat[sample.null.list, "CDScount"])
		exp_gene = sum(data.geneNum.merge.stat[sample.exp.list, "CDScount"])
		null_COG_gene = length(unique(tmp.df[which(tmp.df$assembly %in% sample.null.list),"CDS"]))
		exp_COG_gene = length(unique(tmp.df[which(tmp.df$assembly %in% sample.exp.list),"CDS"]))
    	testor_species <- rbind(c(exp_COG, exp_species - exp_COG), c(null_COG, null_species - null_COG))
    	model_species <- fisher.test(testor_species, alternative="greater")
		testor_gene <- rbind(c(exp_COG_gene, exp_gene - exp_COG_gene), c(null_COG_gene, null_gene - null_COG_gene))
    	model_gene <- fisher.test(testor_gene, alternative="greater")
    	COG.enrich.info[i,] <- c(tmp.COG.term, tmp.COG.note, null_species, exp_species, null_COG, exp_COG, 
							null_gene, exp_gene, null_COG_gene, exp_COG_gene,
							model_species$p.value, model_species$estimate,
							model_gene$p.value, model_gene$estimate)
	}
	COG.enrich.info$pvalue_species <- as.numeric(COG.enrich.info$pvalue_species)
	COG.enrich.info$odds_species <- as.numeric(COG.enrich.info$odds_species)
	COG.enrich.info$pvalue_gene <- as.numeric(COG.enrich.info$pvalue_gene)
	COG.enrich.info$odds_gene <- as.numeric(COG.enrich.info$odds_gene)
	COG.enrich.info$padj_species <- p.adjust(COG.enrich.info$pvalue_species, method = "BH")
	COG.enrich.info$padj_gene <- p.adjust(COG.enrich.info$pvalue_gene, method = "BH")
	Pvalue_threshold <- 0.01
	COG.enrich.info$sig_gene <- "Not Sig"
	COG.enrich.info[which(COG.enrich.info$padj_gene < Pvalue_threshold), "sig_gene"] <- "Sig"
	COG.enrich.info$sig_species <- "Not Sig"
	COG.enrich.info[which(COG.enrich.info$padj_species < Pvalue_threshold), "sig_species"] <- "Sig"

	p1 <- ggplot(COG.enrich.info, aes(x = log2(odds_gene), y = -log10(padj_gene))) +
	geom_point(aes(color = sig_gene), alpha = 0.3, size = 1) +
	geom_text_repel(data = COG.enrich.info[which(COG.enrich.info$sig_gene == "Sig"),], aes(label = paste(COG_term, ": ", COG_notes, sep = "")), size = 1.5, color = "#a50f15") +
	scale_color_manual(values=c("black", "red")) +
	labs(title = "COG enrichment on gene (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
	theme(legend.position="none") +
	xlim(-5.8, 6.2) + ylim(0, 5) + 
	theme(plot.title = element_text(hjust = 0.5))
	p2 <- ggplot(COG.enrich.info, aes(x = log2(odds_species), y = -log10(padj_species))) +
	geom_point(aes(color = sig_species), alpha = 0.3, size = 1) +
	geom_text_repel(data = COG.enrich.info[which(COG.enrich.info$sig_species == "Sig"),], aes(label = paste(COG_term, ": ", COG_notes, sep = "")), size = 1.5, color = "#a50f15") +
	scale_color_manual(values=c("black", "red")) +
	labs(title = "COG enrichment on assembly (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
	theme(legend.position="none") +
	xlim(-5.8, 6.2) + ylim(0, 5) + 
	theme(plot.title = element_text(hjust = 0.5))
    plotVector[[(j - 1) * 2 + 1]] <- p1
    plotVector[[(j - 1) * 2 + 2]] <- p2
}
pdf("./figures/01.COG_analysis.cluster.individualTest.pdf", width = 15, height = 8)
multiplot(plotlist = plotVector, cols = 3)
dev.off()





