library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)
library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggdendro)
library(ape)
library(ggforce)
library(RColorBrewer)
library(vegan)
library(ggrepel)

ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}



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

data.merge <- read.table("./assembly.KEGG.txt", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.merge) <- c("assembly", "CDS", "KEGG", "evalue")
data.merge <- data.merge[which(!is.na(data.merge$KEGG)),]
data.geneNum.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly) %>% 
                            summarise(CDScount = length(unique(CDS))),
                            stringsAsFactors = F)
rownames(data.geneNum.merge.stat) <- data.geneNum.merge.stat$assembly
data.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly, KEGG) %>% 
                            summarise(count = n()),
                            stringsAsFactors = F)
data.merge.stat.dcast <- dcast(data.merge.stat, assembly ~ KEGG, fill = 0)
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
								phylum = data.tax.info[rownames(data.supplier.info),"phylum"],
								order = data.tax.info[rownames(data.supplier.info),"order"],
								family = data.tax.info[rownames(data.supplier.info),"family"],
								row.names = rownames(data.supplier.info),
								stringsAsFactors = F)



pathway.annotation <- read.table("./archive/KEGG_pathway.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(pathway.annotation) <- c("A_term", "A_title", "B_term", "B_title", "Pathway_term", "Pathway_title", "KO", "gene", "annotation")
pathway.info <- data.frame(pathway.annotation %>%
						group_by(Pathway_term, Pathway_title) %>%
						summarize(KO_num = n()), stringsAsFactors = F)
assembly.pathway.matrix <- matrix(rep(0, nrow(data.merge.stat.dcast) * nrow(pathway.info)), nrow = nrow(data.merge.stat.dcast))
colnames(assembly.pathway.matrix) <- pathway.info$Pathway_term
rownames(assembly.pathway.matrix) <- rownames(data.merge.stat.dcast)
for (i in 1:nrow(pathway.info)){
	tmp.pathway.term <- pathway.info[i, "Pathway_term"]
	tmp.pathway.notes <- pathway.info[i, "Pathway_title"]
	tmp.pathway.KO <- pathway.annotation[which(pathway.annotation$Pathway_term == tmp.pathway.term),"KO"]
	tmp.pathway.KO.used <- intersect(tmp.pathway.KO, colnames(data.merge.stat.dcast))
	if (length(tmp.pathway.KO.used) == 0){
		tmp.vector <- rep(0, nrow(assembly.pathway.matrix))
	}else{
		if (length(tmp.pathway.KO.used) == 1){
			tmp.vector <- data.merge.stat.dcast[,tmp.pathway.KO.used]
		}else{
			tmp.matrix <- data.merge.stat.dcast[,tmp.pathway.KO.used]
			tmp.vector <- rowSums(tmp.matrix)
		}
	}
	assembly.pathway.matrix[,tmp.pathway.term] <- tmp.vector
}

module.annotation <- read.table("./archive/KEGG_module.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(module.annotation) <- c("A_title", "B_title", "C_title", "Module_term", "Module_title", "KO", "gene", "annotation")
module.info <- data.frame(module.annotation %>%
                        group_by(Module_term, Module_title) %>%
                        summarize(KO_num = n()), stringsAsFactors = F)
assembly.module.matrix <- matrix(rep(0, nrow(data.merge.stat.dcast) * nrow(module.info)), nrow = nrow(data.merge.stat.dcast))
colnames(assembly.module.matrix) <- module.info$Module_term
rownames(assembly.module.matrix) <- rownames(data.merge.stat.dcast)
for (i in 1:nrow(module.info)){
	tmp.module.term <- module.info[i, "Module_term"]
	tmp.module.notes <- module.info[i, "Module_title"]
	tmp.module.KO <- module.annotation[which(module.annotation$Module_term == tmp.module.term),"KO"]
	tmp.module.KO.used <- intersect(tmp.module.KO, colnames(data.merge.stat.dcast))
	if (length(tmp.module.KO.used) == 0){
		tmp.vector <- rep(0, nrow(assembly.module.matrix))
	}else{
		if (length(tmp.module.KO.used) == 1){
			tmp.vector <- data.merge.stat.dcast[,tmp.module.KO.used]
		}else{
			tmp.matrix <- data.merge.stat.dcast[,tmp.module.KO.used]
			tmp.vector <- rowSums(tmp.matrix)
		}
	}
	assembly.module.matrix[,tmp.module.term] <- tmp.vector
}


### MDS for all assembly ###
sample.list <- intersect(assembly.metadata.df[,"assembly"], rownames(assembly.pathway.matrix))

data.matrix.used <- assembly.pathway.matrix[sample.list,]
data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							phylum = assembly.metadata.df[rownames(data.merge.mds$points),"phylum"],
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            stringsAsFactors = F)
family.order <- read.table("./family.order", header = F, sep = "\t", stringsAsFactors = F)$V1
data.mds.result$family.show <- data.mds.result$family
data.mds.result[which(!(data.mds.result$family.show %in% family.order)), "family.show"] <- "Others"

pdf("./figures/03.KEGG_analysis.MDS_pathway_all.pdf", height = 3.2, width = 5.5)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = factor(family.show, levels = family.order), shape = supplier)) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_shape_manual(values = c(1, 3, 16, 8)) +
    theme_bw() +
	scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(12)) +
    labs(x = "MDS1 (KEGG pathway)", y = "MDS2 (KEGG pathway)", color = "family") +
	theme(axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) 
dev.off()



### MDS for all Muribaculaceae ###
sample.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"]

data.matrix.used <- assembly.pathway.matrix[sample.list,]
data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							phylum = assembly.metadata.df[rownames(data.merge.mds$points),"phylum"],
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            stringsAsFactors = F)

family.order <- read.table("./family.order", header = F, sep = "\t", stringsAsFactors = F)$V1
data.mds.result$family.show <- data.mds.result$family
data.mds.result[which(!(data.mds.result$family.show %in% family.order)), "family.show"] <- "Others"

pdf("./figures/03.KEGG_analysis.MDS_pathway_Muribaculaceae.pdf", height = 2.7, width = 4)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = supplier, shape = supplier)) +
    geom_point(alpha = 0.9, size = 1.3) +
	scale_shape_manual(values = c(1, 3, 16, 8)) +
	theme_bw() +
    scale_color_manual(values = colorRampPalette(brewer.pal(4, "Dark2"))(4)) +
    labs(x = "MDS1 (KEGG pathway)", y = "MDS2 (KEGG pathway)") +
	theme(axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) 
dev.off()



pathway.annotation <- read.table("./archive/KEGG_pathway.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(pathway.annotation) <- c("A_term", "A_title", "B_term", "B_title", "Pathway_term", "Pathway_title", "KO", "gene", "annotation")
pathway.annotation <- pathway.annotation[which(pathway.annotation$A_title == "Metabolism"), ]

pathway.info <- data.frame(pathway.annotation %>%
						group_by(Pathway_term, Pathway_title) %>%
						summarize(KO_num = n()), stringsAsFactors = F)


module.annotation <- read.table("./archive/KEGG_module.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(module.annotation) <- c("A_title", "B_title", "C_title", "Module_term", "Module_title", "KO", "gene", "annotation")
module.info <- data.frame(module.annotation %>%
						group_by(Module_term, Module_title) %>%
						summarize(KO_num = n()), stringsAsFactors = F)

KO.annotation <- read.table("./archive/KEGG_KO.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(KO.annotation) <- c("A_term", "A_title", "B_term", "B_title", "Pathway_term", "Pathway_title", "KO", "gene", "annotation")
KO.info <- data.frame(KO.annotation %>%
					group_by(KO, annotation) %>%
					summarize(duplicates = n()),
					stringsAsFactors = F)
rownames(KO.info) <- KO.info$KO



sample.used <- assembly.metadata.df[, "assembly"]
data.dataframe.used <- data.merge

pathway.enrich.info <- data.frame()

for (i in 1:nrow(pathway.info)){
	print(i)
	tmp.pathway.term <- pathway.info[i, "Pathway_term"]
	tmp.pathway.notes <- pathway.info[i, "Pathway_title"]
	tmp.pathway.KO <- pathway.annotation[which(pathway.annotation$Pathway_term == tmp.pathway.term),"KO"]
	tmp.pathway.allGene <- data.dataframe.used[which(data.dataframe.used$KEGG %in% tmp.pathway.KO),]
	if (nrow(tmp.pathway.allGene) >= 1){
		tmp.pathway.allGene.stat <- data.frame(tmp.pathway.allGene %>%
									group_by(assembly) %>%
									summarise(count = n()),
									stringsAsFactors = F)
		tmp.pathway.allGene.stat$ratio <- tmp.pathway.allGene.stat$count / length(tmp.pathway.KO)
		tmp.pathway.allGene.stat$pathway.term <- tmp.pathway.term
		tmp.pathway.allGene.stat$pathway.note <- tmp.pathway.notes
		pathway.enrich.info <- rbind(pathway.enrich.info, tmp.pathway.allGene.stat)
	}
}

pathway.enrich.info$supplier <- assembly.metadata.df[as.character(pathway.enrich.info$assembly), "supplier"]
pathway.enrich.info$abundance <- assembly.metadata.df[as.character(pathway.enrich.info$assembly), "abundance"]

pathway.enrich.supplierStat <- data.frame(pathway.enrich.info %>%
								group_by(supplier, pathway.term, pathway.note) %>%
								summarise(abundance = sum(ratio * abundance),
										count = sum(ratio)),
								stringsAsFactors = F)
pathway.enrich.supplierStat$pathway.show <- paste(pathway.enrich.supplierStat$pathway.term, pathway.enrich.supplierStat$pathway.note, sep = ", ")
pathway.enrich.supplierStat.dcast <- dcast(pathway.enrich.supplierStat[,c("supplier", "pathway.show", "abundance")], pathway.show ~ supplier, fill = 0)
rownames(pathway.enrich.supplierStat.dcast) <- pathway.enrich.supplierStat.dcast$pathway.show
pathway.enrich.supplierStat.dcast <- pathway.enrich.supplierStat.dcast[,-1]
pathway.enrich.supplierStat.melt <- melt(as.matrix(pathway.enrich.supplierStat.dcast))
colnames(pathway.enrich.supplierStat.melt) <- c("pathway.show", "supplier", "abundance")

pathway.enrich.supplierStat.mean <- rowMeans(as.matrix(pathway.enrich.supplierStat.dcast))
pathway.enrich.supplierStat.sd <- rowSds(as.matrix(pathway.enrich.supplierStat.dcast))
pathway.enrich.supplierStat.zscore <- pathway.enrich.supplierStat.sd / pathway.enrich.supplierStat.mean
pathway.enrich.supplierStat.zscore <- sort(pathway.enrich.supplierStat.zscore, decreasing = T)

pathway.toshow <- names(pathway.enrich.supplierStat.zscore[1:25])

hc <- hclust(dist(log10(pathway.enrich.supplierStat.dcast[pathway.toshow,] + 0.0001)))
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
py <- ggdend(ddata$segments) + coord_flip()
row.ord <- order.dendrogram(dhc)
row.name <- rownames(pathway.enrich.supplierStat.dcast[pathway.toshow,])[row.ord]


pdf("./figures/00.KEGG_pathway_abundance.all.pdf", width = 8.8, height = 3.5)
ggplot(pathway.enrich.supplierStat.melt[which(pathway.enrich.supplierStat.melt$pathway.show %in% pathway.toshow),],
    aes(x = factor(pathway.show, levels = row.name),
		y = supplier)) +
    geom_tile(aes(fill = log10(abundance + 0.0001)), colour = "black", size = 0.3) + 
    scale_fill_gradient2(limit = c(-4, 0), midpoint = -2, mid = "#DBDDE0", low = "#647CBD", high = "#E0624E") +
	theme(axis.text.x = element_text(color = "black", size = 5.5, hjust = 1, angle = 45),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(x = "Samples", y = "Otu-level taxonomy", fill = "log10(Relative abundance)")
dev.off()


assembly.metadata.supplier.stat <- data.frame(assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),] %>%
								group_by(supplier) %>%
								summarise(count = n()),
								stringsAsFactors = F)



sample.used <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"]
data.dataframe.used <- data.merge

pathway.enrich.info <- data.frame()

for (i in 1:nrow(pathway.info)){
	print(i)
	tmp.pathway.term <- pathway.info[i, "Pathway_term"]
	tmp.pathway.notes <- pathway.info[i, "Pathway_title"]
	tmp.pathway.KO <- pathway.annotation[which(pathway.annotation$Pathway_term == tmp.pathway.term),"KO"]
	tmp.pathway.allGene <- data.dataframe.used[which(data.dataframe.used$KEGG %in% tmp.pathway.KO &
										data.dataframe.used$assembly %in% sample.used),]
	if (nrow(tmp.pathway.allGene) >= 1){
		tmp.pathway.allGene.stat <- data.frame(tmp.pathway.allGene %>%
									group_by(assembly) %>%
									summarise(count = n()),
									stringsAsFactors = F)
		tmp.pathway.allGene.stat$ratio <- tmp.pathway.allGene.stat$count / length(tmp.pathway.KO)
		tmp.pathway.allGene.stat$pathway.term <- tmp.pathway.term
		tmp.pathway.allGene.stat$pathway.note <- tmp.pathway.notes
		pathway.enrich.info <- rbind(pathway.enrich.info, tmp.pathway.allGene.stat)
	}
}

pathway.enrich.info$supplier <- assembly.metadata.df[as.character(pathway.enrich.info$assembly), "supplier"]
pathway.enrich.info$abundance <- assembly.metadata.df[as.character(pathway.enrich.info$assembly), "abundance"]

pathway.enrich.supplierStat <- data.frame(pathway.enrich.info %>%
								group_by(supplier, pathway.term, pathway.note) %>%
								summarise(abundance = sum(ratio * abundance),
										count = sum(ratio)),
								stringsAsFactors = F)

pathway.enrich.supplierStat$pathway.show <- paste(pathway.enrich.supplierStat$pathway.term, pathway.enrich.supplierStat$pathway.note, sep = ", ")
pathway.enrich.supplierStat.dcast <- dcast(pathway.enrich.supplierStat[,c("supplier", "pathway.show", "abundance")], pathway.show ~ supplier, fill = 0)
rownames(pathway.enrich.supplierStat.dcast) <- pathway.enrich.supplierStat.dcast$pathway.show
pathway.enrich.supplierStat.dcast <- pathway.enrich.supplierStat.dcast[,-1]
pathway.enrich.supplierStat.melt <- melt(as.matrix(pathway.enrich.supplierStat.dcast))
colnames(pathway.enrich.supplierStat.melt) <- c("pathway.show", "supplier", "abundance")

pathway.enrich.supplierStat.mean <- rowMeans(as.matrix(pathway.enrich.supplierStat.dcast))
pathway.enrich.supplierStat.sd <- rowSds(as.matrix(pathway.enrich.supplierStat.dcast))
pathway.enrich.supplierStat.zscore <- pathway.enrich.supplierStat.sd / pathway.enrich.supplierStat.mean
pathway.enrich.supplierStat.zscore <- sort(pathway.enrich.supplierStat.zscore, decreasing = T)

pathway.toshow <- names(pathway.enrich.supplierStat.zscore[1:25])

hc <- hclust(dist(log10(pathway.enrich.supplierStat.dcast[pathway.toshow,] + 0.0001)))
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
py <- ggdend(ddata$segments) + coord_flip()
row.ord <- order.dendrogram(dhc)
row.name <- rownames(pathway.enrich.supplierStat.dcast[pathway.toshow,])[row.ord]

pdf("./figures/00.KEGG_pathway_abundance.Muribaculaceae.pdf", width = 8.8, height = 3.8)
ggplot(pathway.enrich.supplierStat.melt[which(pathway.enrich.supplierStat.melt$pathway.show %in% pathway.toshow),],
    aes(x = factor(pathway.show, levels = row.name),
		y = supplier)) +
    geom_tile(aes(fill = log10(abundance + 0.0001)), colour = "black", size = 0.3) + 
    scale_fill_gradient2(limit = c(-4, 0), midpoint = -2, mid = "#DBDDE0", low = "#647CBD", high = "#E0624E") +
    theme(axis.text.x = element_text(color = "black", size = 5.5, hjust = 1, angle = 45),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(x = "Samples", y = "Otu-level taxonomy", fill = "log10(Relative abundance)")
dev.off()







assembly.cluster <- read.table("./assembly.CAZy_GH.cluster.tsv", stringsAsFactors = F, sep = "\t", header = F)
colnames(assembly.cluster) <- c("assembly", "cluster")
rownames(assembly.cluster) <- assembly.cluster$assembly

data.cluster.meta <- data.frame(row.names = assembly.cluster$assembly,
                                assembly = assembly.cluster$assembly,
                                cluster = assembly.cluster$cluster,
                                supplier = assembly.metadata.df[assembly.cluster$assembly,"supplier"],
                                abundance = assembly.metadata.df[assembly.cluster$assembly,"abundance"],
                                stringsAsFactors = F)

pathway.annotation <- read.table("./archive/KEGG_pathway.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(pathway.annotation) <- c("A_term", "A_title", "B_term", "B_title", "Pathway_term", "Pathway_title", "KO", "gene", "annotation")
pathway.info <- data.frame(pathway.annotation %>%
                        group_by(Pathway_term, Pathway_title) %>%
                        summarize(KO_num = n()), stringsAsFactors = F)
module.annotation <- read.table("./archive/KEGG_module.txt", sep = "\t", stringsAsFactors = F, quote = "")
colnames(module.annotation) <- c("A_title", "B_title", "C_title", "Module_term", "Module_title", "KO", "gene", "annotation")
module.info <- data.frame(module.annotation %>%
                        group_by(Module_term, Module_title) %>%
                        summarize(KO_num = n()), stringsAsFactors = F)

data.dataframe.used <- data.merge
cluster.list <- c("CAZy_GH_1", "CAZy_GH_2", "CAZy_GH_3")
plotVector <- vector("list", length = 6)
for (j in 1:3){
    tmpCluster <- cluster.list[j]
    sample.null.list <- data.cluster.meta[data.cluster.meta$cluster != tmpCluster, "assembly"]
    sample.exp.list <- data.cluster.meta[data.cluster.meta$cluster == tmpCluster, "assembly"]

	pathway.enrich.info <- data.frame(Pathway_term = rep(" ", nrow(pathway.info)),
								Pathway_notes = rep(" ", nrow(pathway.info)),
								null_gene = rep(0, nrow(pathway.info)),
								exp_gene = rep(0, nrow(pathway.info)),
								null_KO = rep(0, nrow(pathway.info)),
								exp_KO = rep(0, nrow(pathway.info)),
								pvalue = rep(0, nrow(pathway.info)),
								odds = rep(0, nrow(pathway.info)),
								stringsAsFactors = F)
	for (i in 1:nrow(pathway.info)){
		tmp.pathway.term <- pathway.info[i, "Pathway_term"]
		tmp.pathway.notes <- pathway.info[i, "Pathway_title"]
		tmp.pathway.KO <- pathway.annotation[which(pathway.annotation$Pathway_term == tmp.pathway.term),"KO"]
		tmp.pathway.allGene <- data.dataframe.used[which(data.dataframe.used$KEGG %in% tmp.pathway.KO),]
		null_KO <- length(unique(tmp.pathway.allGene[which(tmp.pathway.allGene$assembly %in% sample.null.list),"CDS"]))
		exp_KO <- length(unique(tmp.pathway.allGene[which(tmp.pathway.allGene$assembly %in% sample.exp.list),"CDS"]))
		null_gene <- sum(data.geneNum.merge.stat[sample.null.list,"CDScount"])
		exp_gene <- sum(data.geneNum.merge.stat[sample.exp.list,"CDScount"])
		testor <- rbind(c(exp_KO, exp_gene - exp_KO), c(null_KO, null_gene - null_KO))
		model <- fisher.test(testor, alternative = "greater")
		pathway.enrich.info[i,] <- c(tmp.pathway.term, tmp.pathway.notes, null_gene, exp_gene, null_KO, exp_KO, model$p.value, model$estimate)
	}
	pathway.enrich.info$pvalue <- as.numeric(pathway.enrich.info$pvalue)
	pathway.enrich.info$odds <- as.numeric(pathway.enrich.info$odds)
	pathway.enrich.info$padj <- p.adjust(pathway.enrich.info$pvalue, method = "BH")
	module.enrich.info <- data.frame(Module_term = rep(" ", nrow(module.info)),
								Module_notes = rep(" ", nrow(module.info)),
								null_gene = rep(0, nrow(module.info)),
								exp_gene = rep(0, nrow(module.info)),
								null_KO = rep(0, nrow(module.info)),
								exp_KO = rep(0, nrow(module.info)),
								pvalue = rep(0, nrow(module.info)),
								odds = rep(0, nrow(module.info)),
								stringsAsFactors = F)
	for (i in 1:nrow(module.info)){
		tmp.module.term <- module.info[i, "Module_term"]
		tmp.module.notes <- module.info[i, "Module_title"]
		tmp.module.KO <- module.annotation[which(module.annotation$Module_term == tmp.module.term),"KO"]
		tmp.module.allGene <- data.dataframe.used[which(data.dataframe.used$KEGG %in% tmp.module.KO),]
		null_KO <- length(unique(tmp.module.allGene[which(tmp.module.allGene$assembly %in% sample.null.list),"CDS"]))
		exp_KO <- length(unique(tmp.module.allGene[which(tmp.module.allGene$assembly %in% sample.exp.list),"CDS"]))
		null_gene <- sum(data.geneNum.merge.stat[sample.null.list,"CDScount"])
		exp_gene <- sum(data.geneNum.merge.stat[sample.exp.list,"CDScount"])
		testor <- rbind(c(exp_KO, exp_gene - exp_KO), c(null_KO, null_gene - null_KO))
		model <- fisher.test(testor, alternative = "greater")
		module.enrich.info[i,] <- c(tmp.module.term, tmp.module.notes, null_gene, exp_gene, null_KO, exp_KO, model$p.value, model$estimate)
	}
	module.enrich.info$pvalue <- as.numeric(module.enrich.info$pvalue)
	module.enrich.info$odds <- as.numeric(module.enrich.info$odds)
	module.enrich.info$padj <- p.adjust(module.enrich.info$pvalue, method = "BH")
	Pvalue_threshold <- 0.05
	pathway.enrich.info$sig <- "Not Sig"
	pathway.enrich.info[which(pathway.enrich.info$padj < Pvalue_threshold), "sig"] <- "Sig"
	module.enrich.info$sig <- "Not Sig"
	module.enrich.info[which(module.enrich.info$padj < Pvalue_threshold), "sig"] <- "Sig"

	p1 <- ggplot(pathway.enrich.info, aes(x = log2(odds), y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.3, size = 1) +
    geom_text_repel(data = pathway.enrich.info[which(pathway.enrich.info$sig == "Sig"),], aes(label = paste(Pathway_term, ": ", Pathway_notes, sep = "")), size = 3, color = "#a50f15") +
    scale_color_manual(values=c("black", "red")) +
    labs(title = "KEGG pathway enrichment on gene (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))

	p2 <- ggplot(module.enrich.info, aes(x = log2(odds), y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.3, size = 1) +
    geom_text_repel(data = module.enrich.info[which(module.enrich.info$sig == "Sig"),], aes(label = paste(Module_term, ": ", Module_notes, sep = "")), size = 3, color = "#a50f15") +
    scale_color_manual(values=c("black", "red")) +
    labs(title = "KEGG module enrichment on gene (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))

    plotVector[[(j - 1) * 2 + 1]] <- p1
    plotVector[[(j - 1) * 2 + 2]] <- p2
}

pdf("./figures/03.KEGG_analysis.cluster.individualTest.pdf", width = 15, height = 8)
multiplot(plotlist = plotVector, cols = 3)
dev.off()





