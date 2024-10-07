library(ggplot2)
library(reshape)
library(dplyr)
library(ggtree)
library(RColorBrewer)

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
								class = data.tax.info[rownames(data.supplier.info),"class"],
                                order = data.tax.info[rownames(data.supplier.info),"order"],
                                family = data.tax.info[rownames(data.supplier.info),"family"],
                                row.names = rownames(data.supplier.info),
                                stringsAsFactors = F)

assembly.metadata.df.stat <- data.frame(assembly.metadata.df %>%
										group_by(supplier, family) %>%
										summarise(assembly_count = n()),
										stringsAsFactors = F)
assembly.metadata.df.stat <- assembly.metadata.df.stat[which(assembly.metadata.df.stat$assembly_count >= 1),]
assembly.metadata.df.stat[which(is.na(assembly.metadata.df.stat$family)), "family"] <- "Others"

assembly.metadata.familysort <- data.frame(assembly.metadata.df.stat %>%
											group_by(family) %>%
											summarise(count = sum(assembly_count)),
											stringsAsFactors = F)
assembly.metadata.familysort <- assembly.metadata.familysort[order(assembly.metadata.familysort$count, decreasing = T),]
rownames(assembly.metadata.familysort) <- assembly.metadata.familysort$family

tax.to.show <- 12
assembly.metadata.familysort$family.show <- assembly.metadata.familysort$family
assembly.metadata.familysort[(tax.to.show):nrow(assembly.metadata.familysort), "family.show"] <- "Others"

family.order <- unique(assembly.metadata.familysort$family.show)

assembly.metadata.df.stat$family.show <- assembly.metadata.familysort[as.character(assembly.metadata.df.stat$family), "family.show"]
assembly.metadata.df.showStat <- data.frame(assembly.metadata.df.stat %>%
									group_by(supplier, family.show) %>%
									summarise(assembly_count = sum(assembly_count)),
									stringsAsFactors = F)

write.table(family.order, "./family.order", row.names = F, col.names = F, quote = F)

pdf("./figures/00.taxa_count.pdf", width = 7, height = 2.7)
ggplot(assembly.metadata.df.showStat, aes(x = supplier, y = assembly_count,
			fill = factor(family.show, levels = family.order))) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.9) + 
    theme_bw() +
	scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(length(family.order))) +
    labs(x = "Suppliers", y = "# of assemblies", fill = "Family-level taxonomy") +
	theme(axis.text.x = element_text(color = "black", size = 9),
		axis.text.y = element_text(color = "black", size = 9))
dev.off()








tax.to.show <- 12
assembly.metadata.familysort$family.show <- assembly.metadata.familysort$family
assembly.metadata.familysort[(tax.to.show):nrow(assembly.metadata.familysort), "family.show"] <- "Others"

family.order <- unique(assembly.metadata.familysort$family.show)

assembly.metadata.df.stat$family.show <- assembly.metadata.familysort[as.character(assembly.metadata.df.stat$family), "family.show"]
ggplot(assembly.metadata.df.stat, aes(x = supplier, y = assembly_count, fill = family.show)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    theme_classic() +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8),"Set1"))(length(family.order)) +                            
    labs(x = "Suppliers", y = "# of assemblies") +
    theme(axis.text.x = element_text(color="black", size = 8))




### abundance for AZC3J ###
data.mapping <- read.table("./output_11_reAssemblyMapping.stat.txt", stringsAsFactors = F, header = T, sep = "\t")
rownames(data.mapping) <- data.mapping$sampleID
mouse.info <- read.table("./sample.info", header = F, stringsAsFactors = F, sep = "\t")
colnames(mouse.info) <- c("sample", "supplier", "mouse", "day")
rownames(mouse.info) <- mouse.info$sample
data.mapping$mouse <- mouse.info[data.mapping$sampleID, "mouse"]
data.mapping$day <- mouse.info[data.mapping$sampleID, "day"]




ggplot(data.mapping, aes(x = supplier, y = totalPairQC)) + 
	geom_violin(fill="gray", trim = F) +
	geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.3) + 
	theme_classic() + 
	scale_y_log10(limits = c(100, 100000000), breaks = scales::trans_breaks("log10", function(x) 10^x),
				labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
	labs(x = "Suppliers", y = "# of reads after QC") + 
	theme(axis.text.x = element_text(color = "black", size = 12))



ggplot(data.mapping, aes(x = supplier, y = totalMappedReads / totalPairNoSpikein / 2 * 100)) + 
	geom_violin(fill="gray", trim = F) + 
	geom_boxplot(position = position_dodge(0.9), width = 0.15) + 
	theme_classic() + ylim(0, 100) + 
	labs(x = "Suppliers", y = "Overall alignment rate(%)") + 
	theme(axis.text.x = element_text(color = "black", size = 12))


data.weight <- read.table("./archive/sample_weight.tsv", stringsAsFactors = F, header = T, sep = "\t")
rownames(data.weight) <- data.weight$SampleID

data.assembly <- read.table("./output_10_reAssemblyDuplicates.stat", stringsAsFactors = F, header = T, sep = "\t")
rownames(data.assembly) <- data.assembly$assemblyID
checkM.info <- read.table("./output_09_reAssemblyAnnotation/merge.checkM.stat.txt", stringsAsFactors = F, header = F, sep = "\t", quote = "")
rownames(checkM.info) <- checkM.info$V1
data.assembly$Completeness <- checkM.info[data.assembly$assemblyID, "V12"]
data.assembly$Contamination <- checkM.info[data.assembly$assemblyID, "V13"]
data.assembly$Heterogeneity <- checkM.info[data.assembly$assemblyID, "V14"]

data.assembly.stat <- data.frame(data.assembly %>%
								group_by(supplier) %>%
								summarise(assembly_count = n()),
								stringsAsFactors = F)
rownames(data.assembly.stat) <- data.assembly.stat$supplier
data.assembly$assembly_count <- data.assembly.stat[data.assembly$supplier, "assembly_count"]

ggplot(data.assembly, aes(x = paste(supplier, "(", assembly_count, ")", sep = ""), y = Heterogeneity)) + 
	geom_violin(fill="gray", trim = t) +
	geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.3) + 
	theme_classic() + ylim(0, 100) + 
	labs(x = "Suppliers", y = "Strain heterogeneity") + 
	theme(axis.text.x = element_text(color = "black", size = 8))

ggplot(data.assembly, aes(x = paste(supplier, "(", assembly_count, ")", sep = ""), y = largestContig)) + 
	geom_violin(fill="gray", trim = F) +
	geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.3) + 
	theme_classic() + 
	scale_y_log10(limits = c(1, 10000000), breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
	labs(x = "Suppliers", y = "length of largest contigs") + 
	theme(axis.text.x = element_text(color = "black", size = 8))











supplier <- "AZC5T"
data.assembly_AZC3J <- read.table(paste("./archive/mappingStat.", supplier, ".stat", sep = ""), stringsAsFactors = F, header = T, sep = "\t")
rownames(data.assembly_AZC3J) <- data.assembly_AZC3J$sampleID
data.assembly_AZC3J <- data.assembly_AZC3J[,-1]
data.assembly_AZC3J.relative <- sweep(data.assembly_AZC3J, 1, data.mapping[rownames(data.assembly_AZC3J), "totalPairNoSpikein"] * 2, "/")
data.assembly_AZC3J.reAb <- sweep(data.assembly_AZC3J.relative, 2, data.assembly[colnames(data.assembly_AZC3J.relative),"totalLength"] / 1000000, "/")
write.table(t(data.assembly_AZC3J.reAb), paste("assembly.abundance_relative.", supplier, ".stat", sep = ""), row.names = T, quote = F, sep = "\t")


data.assembly_AZC3J.spikeIn <- sweep(data.assembly_AZC3J, 1, data.weight[rownames(data.assembly_AZC3J), "weight"] * data.mapping[rownames(data.assembly_AZC3J),"spikeInMappedReads"], "/")
data.assembly_AZC3J.abAb <- sweep(data.assembly_AZC3J.spikeIn, 2, data.assembly[colnames(data.assembly_AZC3J.spikeIn),"totalLength"] / 1000000, "/")
write.table(t(data.assembly_AZC3J.abAb), paste("assembly.abundance.", supplier, ".stat", sep = ""), row.names = T, quote = F, sep = "\t")



data.assembly_AZC3J.abAb.perA <- colMeans(data.assembly_AZC3J.abAb)
write.table(data.frame(assembly = names(data.assembly_AZC3J.abAb.perA), AbMean = data.assembly_AZC3J.abAb.perA), "assembly.abundance.AZC7E_mean.stat", row.names = F, quote = F, sep = "\t")

sample.AZC3J.list <- read.table("./sample_day.AZC3J.list", stringsAsFactors = F)$V1
data.assembly_AZC3J.abAb.melt <- melt(as.matrix(data.assembly_AZC3J.absoluteAbundance))
ggplot(data.assembly_AZC3J.abAb.melt, aes(factor(Var1, levels = sample.AZC3J.list), Var2)) + 
	geom_tile(aes(fill = value), colour = "white") + 
	scale_fill_gradient(limits = c(0, 4), low = "white", high = "firebrick") +
	theme(axis.text.x = element_text(color = "black", size = 5, hjust = 1, angle = 45),
		axis.text.y = element_text(color = "black", size = 3)) +
	labs(x = "Samples", y = "Bins", fill = "Absolute abundance")
############################



### ANI heatmap for Muribaculaceae ###
sample.Muribaculaceae.list <- read.table("./assembly_Muribaculaceae.list", stringsAsFactors = F, header = F)$V1
data.ani <- read.table("./fastANI.merge.txt", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.ani) <- c("assembly1", "assembly2", "ANI")
data.Muribaculaceae.ani <- data.ani[which(data.ani$assembly1 %in% sample.Muribaculaceae.list & 
								data.ani$assembly2 %in% sample.Muribaculaceae.list),]

ggplot(data = data.Muribaculaceae.ani, aes(
        x = factor(assembly1, levels = sample.Muribaculaceae.list),
        y = factor(assembly2, levels = sample.Muribaculaceae.list))) +
    geom_tile(aes(fill = ANI), colour = "white") + 
    scale_fill_gradient(low = "white", high = "firebrick") +
    theme(axis.text.x = element_text(color = "black", size = 4, hjust = 1, angle = 45),
        axis.text.y = element_text(color = "black", size = 5)) +
    labs(y = "Assembly2", x = "Assembly1", fill = "ANI")
######################################


### abundance bar plot for Bacteroidales ###
sample.Bacteroidales.list <- read.table("./assembly_Bacteroidales.list", stringsAsFactors = F, header = F)$V1
data.assembly.abundance <- read.table("./assembly.abundanceMean.stat", stringsAsFactors = F, header = F)
colnames(data.assembly.abundance) <- c("assembly","abundance")
rownames(data.assembly.abundance) <- data.assembly.abundance$assembly
data.supplier <- read.table("./analysis_00_CAZy/Bacteroidales.supplier.txt", stringsAsFactors = F, header = F)
colnames(data.supplier) <- c("assembly", "supplier")
rownames(data.supplier) <- data.supplier$assembly
data.family <- read.table("./analysis_00_CAZy/Bacteroidales.family.txt", stringsAsFactors = F, header = F)
colnames(data.family) <- c("assembly", "family")
rownames(data.family) <- data.family$assembly
data.Bacteroidales.metadata <- data.frame(assembly = sample.Bacteroidales.list, 
										abundance = data.assembly.abundance[sample.Bacteroidales.list, "abundance"],
										supplier = data.supplier[sample.Bacteroidales.list, "supplier"],
										family = data.family[sample.Bacteroidales.list, "family"],
										stringsAsFactors = F)
data.Bacteroidales.metadata.stat <- data.frame(data.Bacteroidales.metadata %>%
											group_by(supplier, family) %>%
											summarise(abundance = sum(abundance)),
											stringsAsFactors = F)

family.list <- c("Muribaculaceae", "Bacteroidaceae", "Rikenellaceae", "Marinifilaceae", "UBA932", "UBA933")
ggplot(data = data.Bacteroidales.metadata.stat, aes(x = supplier, y = abundance, fill = factor(family, levels = rev(family.list)))) +
	geom_bar(stat = "identity", width = 0.7) + 
	theme_classic() +
	scale_fill_manual(values = rev(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))) + 
	theme(axis.text.x = element_text(color = "black", size = 10),
		axis.text.y = element_text(color = "black", size = 8)) +
	labs(y = "Absolute abundance", x = "Supplier", fill = "family", title = "Order: Bacteroidales")
#############################################




### tree analysis for Muribaculaceae ###
data.tree <- read.newick("./merge.align.newick")
data.assembly.abundance <- read.table("./assembly.abundanceMean.stat", stringsAsFactors = F, header = F)
colnames(data.assembly.abundance) <- c("assembly","abundance")
rownames(data.assembly.abundance) <- data.assembly.abundance$assembly
data.family <- read.table("./analysis_00_CAZy/Bacteroidales.family.txt", stringsAsFactors = F, header = F)
colnames(data.family) <- c("assembly", "family")
rownames(data.family) <- data.family$assembly
data.supplier <- read.table("./analysis_00_CAZy/Bacteroidales.supplier.txt", stringsAsFactors = F, header = F)
colnames(data.supplier) <- c("assembly", "supplier")
rownames(data.supplier) <- data.supplier$assembly
p <- ggtree(data.tree, size = 0.5) +
		geom_tiplab(size = 2) +
		theme_tree2()

data.abundance.panel <- data.frame(id = data.tree$tip.label, 
								abundance = data.assembly.abundance[data.tree$tip.label,"abundance"],
								stringsAsFactors = F)
data.supplier.panel <- data.frame(id = data.tree$tip.label,
								supplier = data.supplier[data.tree$tip.label,"supplier"],
								stringsAsFactors = F)

data.supplier.matrix <- data.frame(row.names = data.supplier.panel$id, supplier = data.supplier.panel[,c("supplier")])

p2 <- facet_plot(p, panel = "absolute abundance", data = data.abundance.panel, geom=geom_segment,
				aes(x = 0, xend = abundance, y = y, yend = y), size = 3, color = "steelblue")


gheatmap(p, data.supplier.matrix, offset = 0.1, width = 0.1, font.size = 3) +
	scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"))
#########################################













