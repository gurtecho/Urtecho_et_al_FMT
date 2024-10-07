library(tidyverse)
library(ggrepel)
library(reshape2)
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
library(cowplot)
library(ggtree)
library(rstatix)
library(dunn.test)

setwd("/Users/guillaumeurtecho/Dropbox/fmt/Urtecho_et_al_FMT/experiments/20220507_AZmice_metagenome_vender_zip")

ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}


data.merge <- read.table("./assembly.CAZy_GH.txt", stringsAsFactors = F, header = F, sep = "\t")
colnames(data.merge) <- c("assembly", "CDS", "CAZy_GH", "evalue")
data.geneNum.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly) %>% 
                            summarise(CDScount = length(unique(CDS))),
                            stringsAsFactors = F)
rownames(data.geneNum.merge.stat) <- data.geneNum.merge.stat$assembly
data.merge.stat <- data.frame(data.merge %>% 
                            group_by(assembly, CAZy_GH) %>% 
                            summarise(count = n()),
                            stringsAsFactors = F)
data.merge.stat.dcast <- dcast(data.merge.stat, assembly ~ CAZy_GH, fill = 0)
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

assembly.cluster <- read.table("./assembly.CAZy_GH.cluster.tsv", stringsAsFactors = F, sep = "\t", header = F)
colnames(assembly.cluster) <- c("assembly", "cluster")
rownames(assembly.cluster) <- assembly.cluster$assembly

data.cluster.meta <- data.frame(row.names = assembly.cluster$assembly,
								assembly = assembly.cluster$assembly,
                                cluster = assembly.cluster$cluster,
                                supplier = assembly.metadata.df[assembly.cluster$assembly,"supplier"],
                                abundance = assembly.metadata.df[assembly.cluster$assembly,"abundance"],
                                stringsAsFactors = F)



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
family.order <- read.table("./family.order", header = F, sep = "\t", stringsAsFactors = F)$V1
data.mds.result$family.show <- data.mds.result$family
data.mds.result[which(!(data.mds.result$family.show %in% family.order)), "family.show"] <- "Others"

#pdf("./figures/02.CAZy_GH_analysis.MDS_all.pdf", height = 3.2, width = 5.5)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = factor(family.show, levels = family.order), shape = supplier)) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_shape_manual(values = c(1, 3, 16, 8)) +
    theme_bw() +
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(12)) +
    labs(x = "MDS-1 (CAZymes GH)", y = "MDS-2 (CAZymes GH)", color = "family") +
    theme(axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9))

#pdf("./figures/02.CAZy_GH_analysis.MDS_all.pdf", height = 3.2, width = 5.5)

dev.off() 





### MDS for all Muribaculaceae ###
sample.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"] 
data.matrix.used <- data.merge.stat.dcast.normalize[sample.list,] #%>% select(GH.toShow)

data.merge.cor <- cor(t(data.matrix.used), method = "spearman")
data.merge.dist = sqrt(2 - 2 * data.merge.cor)
data.merge.mds <- cmdscale(data.merge.dist, eig = TRUE, add = TRUE, k = 2)
data.mds.result <- data.frame(assembly = rownames(data.merge.mds$points),
							MDS1 = data.merge.mds$points[,1], MDS2 = data.merge.mds$points[,2],
                            supplier = assembly.metadata.df[rownames(data.merge.mds$points),"supplier"], 
							order = assembly.metadata.df[rownames(data.merge.mds$points),"order"],
                            family = assembly.metadata.df[rownames(data.merge.mds$points),"family"],
                            abundance = assembly.metadata.df[rownames(data.merge.mds$points),"abundance"],    
                            cluster = data.cluster.meta[rownames(data.merge.mds$points),"cluster"],
							stringsAsFactors = F)


family.order <- read.table("./family.order", header = F, sep = "\t", stringsAsFactors = F)$V1
data.mds.result$family.show <- data.mds.result$family
data.mds.result[which(!(data.mds.result$family.show %in% family.order)), "family.show"] <- "Others"


#pdf("./figures/03.CAZy_GH_analysis.MDS_Muribaculaceae.pdf", height = 2.7, width = 4)
ggplot(data.mds.result, aes(x = MDS1, y = MDS2, color = supplier, shape = supplier)) +
    geom_point(alpha = 0.9, size = 1.3) +
    scale_shape_manual(values = c(1, 3, 16, 8)) +
    theme_bw() +
    scale_color_manual(values = colorRampPalette(brewer.pal(4, "Dark2"))(4)) +
    labs(x = "MDS1 (KEGG pathway)", y = "MDS2 (KEGG pathway)") +
    theme(axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9))
dev.off()

# Perform PERMANOVA test, th
permanova_result <- vegan::adonis2(data.merge.dist ~ supplier, data = data.mds.result)
print(permanova_result)

#statistics on this
data_long <- data.merge.stat.dcast.normalize %>%
  rownames_to_column(var = "assembly") %>%
  left_join(assembly.metadata.df %>% select(assembly, supplier), by = "assembly") %>%
  pivot_longer(cols = starts_with("GH"), 
               names_to = "GH", 
               values_to = "Hits") %>%
  mutate(supplier = factor(supplier)) 


# ANOVA or Kruskal-Wallis test
# Replace anova_test with kruskal_test if assumptions for ANOVA are not met
anova_results <- data_long %>%
  group_by(GH) %>%
  anova_test(Hits ~ supplier) %>%
  get_anova_table()

anova_results

# If ANOVA is significant, run post-hoc Tukey's HSD
# For Kruskal-Wallis, replace tukey_hsd with dunn_test
post_hoc_results <- data_long %>%
  group_by(GH) %>%
  tukey_hsd(Hits ~ supplier) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

post_hoc_results 

# Filter for comparisons involving AZC7E and then identify where AZC7E has significantly higher estimates
significant_GH_AZC7E <- post_hoc_results %>%
  filter((group1 == "AZC7E" | group2 == "AZC7E") & p.adj.signif != "ns") %>%
  mutate(significant_in_AZC7E = ifelse(group1 == "AZC7E" & estimate > 0 | group2 == "AZC7E" & estimate < 0, TRUE, FALSE)) %>%
  filter(significant_in_AZC7E == TRUE) %>%
  select(GH, group1, group2, estimate, p.adj)


#############################

sample.used <- assembly.metadata.df[, "assembly"]

used.merge.stat.dcast <- data.merge.stat.dcast[sample.used,]
used.merge.stat.melt <- melt(as.matrix(used.merge.stat.dcast))
colnames(used.merge.stat.melt) <- c("assembly", "GHterm", "count")
used.merge.stat.melt$supplier <- assembly.metadata.df[as.character(used.merge.stat.melt$assembly), "supplier"]
used.merge.stat.melt$abundance <- assembly.metadata.df[as.character(used.merge.stat.melt$assembly), "abundance"]

used.merge.stat.melt$count.binary <- 0
used.merge.stat.melt[which(used.merge.stat.melt$count > 0), "count.binary"] <- 1

used.merge.stat.supplierStat <- data.frame(used.merge.stat.melt %>%
								group_by(supplier, GHterm) %>%
								summarise(ab = sum(count.binary * abundance)),
								stringsAsFactors = F)

used.merge.stat.supplierStat.dcast <- dcast(used.merge.stat.supplierStat, GHterm ~ supplier, fill = 0)
rownames(used.merge.stat.supplierStat.dcast) <- used.merge.stat.supplierStat.dcast$GHterm
used.merge.stat.supplierStat.dcast <- used.merge.stat.supplierStat.dcast[,-1]
used.merge.stat.supplierStat.melt <- melt(as.matrix(used.merge.stat.supplierStat.dcast))
colnames(used.merge.stat.supplierStat.melt) <- c("GHterm", "supplier", "ab")

used.merge.stat.supplierStat.mean <- rowMeans(as.matrix(used.merge.stat.supplierStat.dcast))
used.merge.stat.supplierStat.sd <- rowSds(as.matrix(used.merge.stat.supplierStat.dcast))
used.merge.stat.supplierStat.zscore <- used.merge.stat.supplierStat.sd / used.merge.stat.supplierStat.mean
used.merge.stat.supplierStat.zscore <- sort(used.merge.stat.supplierStat.zscore, decreasing = T)

GH.toShow <- names(used.merge.stat.supplierStat.zscore[1:25])

hc <- hclust(dist(log10(used.merge.stat.supplierStat.dcast[GH.toShow,] + 0.0001)))
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
py <- ggdend(ddata$segments) + coord_flip()
row.ord <- order.dendrogram(dhc)
row.name <- rownames(used.merge.stat.supplierStat.dcast[GH.toShow,])[row.ord]


#pdf("./figures/00.CAZy_abundance.all.pdf", width = 7.8, height = 2)
ggplot(used.merge.stat.supplierStat.melt[which(used.merge.stat.supplierStat.melt$GHterm %in% GH.toShow),],
    aes(x = factor(GHterm, levels = row.name),
        y = supplier)) +
    geom_tile(aes(fill = log10(ab + 0.0001)), colour = "black", size = 0.3) +
    scale_fill_gradient2(limit = c(-4, 1), midpoint = -1.5, mid = "#DBDDE0", low = "#647CBD", high = "#E0624E") +
    theme(axis.text.x = element_text(color = "black", size = 5.5, hjust = 1, angle = 45),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(x = "Samples", y = "Otu-level taxonomy", fill = "log10(Relative abundance)")
dev.off()



sample.used <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae"),"assembly"]

used.merge.stat.dcast <- data.merge.stat.dcast[sample.used,]
used.merge.stat.melt <- melt(as.matrix(used.merge.stat.dcast))
colnames(used.merge.stat.melt) <- c("assembly", "GHterm", "count")
used.merge.stat.melt$supplier <- assembly.metadata.df[as.character(used.merge.stat.melt$assembly), "supplier"]
used.merge.stat.melt$abundance <- assembly.metadata.df[as.character(used.merge.stat.melt$assembly), "abundance"]

used.merge.stat.melt$count.binary <- 0
used.merge.stat.melt[which(used.merge.stat.melt$count > 0), "count.binary"] <- 1

used.merge.stat.supplierStat <- data.frame(used.merge.stat.melt %>%
								group_by(supplier, GHterm) %>%
								summarise(ab = sum(count.binary * abundance)),
								stringsAsFactors = F)

used.merge.stat.supplierStat.dcast <- dcast(used.merge.stat.supplierStat, GHterm ~ supplier, fill = 0)
rownames(used.merge.stat.supplierStat.dcast) <- used.merge.stat.supplierStat.dcast$GHterm
used.merge.stat.supplierStat.dcast <- used.merge.stat.supplierStat.dcast[,-1]
used.merge.stat.supplierStat.melt <- melt(as.matrix(used.merge.stat.supplierStat.dcast))
colnames(used.merge.stat.supplierStat.melt) <- c("GHterm", "supplier", "ab")

used.merge.stat.supplierStat.mean <- rowMeans(as.matrix(used.merge.stat.supplierStat.dcast))
used.merge.stat.supplierStat.sd <- rowSds(as.matrix(used.merge.stat.supplierStat.dcast))
used.merge.stat.supplierStat.zscore <- used.merge.stat.supplierStat.sd / used.merge.stat.supplierStat.mean
used.merge.stat.supplierStat.zscore <- sort(used.merge.stat.supplierStat.zscore, decreasing = T)

GH.toShow <- names(used.merge.stat.supplierStat.zscore[1:50])

hc <- hclust(dist(log10(used.merge.stat.supplierStat.dcast[GH.toShow,] + 0.0001)))
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
py <- ggdend(ddata$segments) + coord_flip()
row.ord <- order.dendrogram(dhc)
row.name <- rownames(used.merge.stat.supplierStat.dcast[GH.toShow,])[row.ord]


#pdf("./figures/00.CAZy_abundance.Muribaculaceae.pdf", width = 10, height = 2)
ggplot(used.merge.stat.supplierStat.melt[which(used.merge.stat.supplierStat.melt$GHterm %in% GH.toShow),],
    aes(x = factor(GHterm, levels = row.name),
        y = supplier)) +
    geom_tile(aes(fill = log10(ab + 0.0001)), colour = "black", size = 0.3) +
    scale_fill_gradient2(limit = c(-4, 1), midpoint = -1.5, mid = "#DBDDE0", low = "#647CBD", high = "#E0624E") +
  theme_cowplot() +
    theme(axis.text.x = element_text(color = "black", size = 5.5, hjust = 1, angle = 45),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(x = "Samples", y = "Otu-level taxonomy", fill = "log10(Relative abundance)") 
dev.off()



##################GU here#####################

#Can we do anova on these?

#just ones shown, slightly significant! Yay.
anova_result <- aov(ab ~ supplier, data = used.merge.stat.supplierStat.melt[which(used.merge.stat.supplierStat.melt$GHterm %in% GH.toShow),])
summary(anova_result)

tukey_result <- TukeyHSD(anova_result)
summary(tukey_result)

#all? This is more significant but eh
anova_result <- aov(ab ~ supplier, data = used.merge.stat.supplierStat.melt)
summary(anova_result) 


#Now let's see if Env MAGs have more GH families than others on average

GH_counts <- data.merge %>% separate(assembly, into = c('supplier', 'MAGno'), remove = F) %>% group_by(assembly) %>% 
#  filter(., grepl("GH",CAZy_GH)) %>%
  mutate(num_hits = n()) %>%
  ungroup %>% 
  select(assembly, supplier, num_hits) %>% distinct() %>%
  left_join(., data.tax.info) 

plot <- GH_counts %>% filter(family == 'Muribaculaceae') %>%
  ggplot(., aes(x=supplier, y=num_hits)) + 
      stat_boxplot(geom = 'errorbar', width = .15) +
      geom_boxplot(outlier.shape = NA, width = .35, fill = 'white') +
      geom_jitter(aes(color = family), position = position_jitter(width = 0.2)) +
    #  scale_color_brewer(palette = "Set1") +
      theme_cowplot() +
      geom_signif(comparisons = list(c("AZC3J", "AZC7E")), 
                  map_signif_level = FALSE,
                  test = "wilcox.test") + 
      labs(x = '', y = 'Num GH hits') +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14), 
            legend.title = element_blank())


plotly::ggplotly(plot)

#How many unique families in each cohort?

data.merge %>%
  separate(assembly, into = c('supplier', 'MAGno'), remove = F) %>% 
  left_join(., data.tax.info) %>% filter(order == 'Bacteroidales') %>%
  select(supplier, CAZy_GH) %>% 
  distinct() %>% #write.csv("venn.GH.csv", quote = F, row.names = F)
  group_by(supplier) %>% 
  mutate(num_GH = n()) %>%
  ungroup %>% 
  select(supplier, num_GH) %>% distinct() 

#How many unique CDS in each cohort?

data.merge %>%
  separate(assembly, into = c('supplier', 'MAGno'), remove = F) %>% 
  left_join(., data.tax.info) %>% filter(order == 'Bacteroidales') %>%
  select(supplier, CDS) %>% group_by(supplier, CDS) %>% mutate(num_duplicates = n()) %>% ungroup %>%
  distinct() %>% #write.csv("venn.GH.csv", quote = F, is row.names = F)
  group_by(supplier) %>% 
  mutate(num_GH = n()) %>%
  ungroup %>% 
  select(supplier, num_GH) %>% distinct() 

#how many duplicated vs non-duplicated CAZymes?

temp <- data.merge %>%
  separate(assembly, into = c('supplier', 'MAGno'), remove = F) %>% 
  left_join(., data.tax.info) %>% filter(order == 'Bacteroidales') %>%
  select(supplier, CDS) %>% group_by(supplier, CDS) %>% mutate(num_duplicates = n()) %>% ungroup %>%
  distinct() %>% mutate()

 

### CAZy_GH enrichment analysis ###
sample.null.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae" & assembly.metadata.df$supplier != "AZC7E"),"assembly"]
sample.exp.list <- assembly.metadata.df[which(assembly.metadata.df$family == "Muribaculaceae" & assembly.metadata.df$supplier == "AZC7E"),"assembly"]

GH.enrich.info <- data.frame(GH = rep(" ", ncol(data.merge.stat.dcast)),
						null_species = rep(0, ncol(data.merge.stat.dcast)),
						exp_species = rep(0, ncol(data.merge.stat.dcast)),
						null_GH = rep(0, ncol(data.merge.stat.dcast)),
						exp_GH = rep(0, ncol(data.merge.stat.dcast)),
						null_gene = rep(0, ncol(data.merge.stat.dcast)),
						exp_gene = rep(0, ncol(data.merge.stat.dcast)),
						null_GH_gene = rep(0, ncol(data.merge.stat.dcast)),
						exp_GH_gene = rep(0, ncol(data.merge.stat.dcast)),
						pvalue_species = rep(0, ncol(data.merge.stat.dcast)),
						odds_species = rep(0, ncol(data.merge.stat.dcast)),
						pvalue_gene = rep(0, ncol(data.merge.stat.dcast)),
                        odds_gene = rep(0, ncol(data.merge.stat.dcast)),
						stringsAsFactors = F)                                          
for (i in 1:ncol(data.merge.stat.dcast)){
	tmpGH = colnames(data.merge.stat.dcast)[i]
   	null_species = length(sample.null.list)
    exp_species = length(sample.exp.list)
	total_assembly = rownames(data.merge.stat.dcast)[which(data.merge.stat.dcast[,tmpGH] > 0)]
	null_GH = length(total_assembly[which(total_assembly %in% sample.null.list)])
	exp_GH = length(total_assembly[which(total_assembly %in% sample.exp.list)])
	null_gene = sum(data.geneNum.merge.stat[sample.null.list, "CDScount"])
	exp_gene = sum(data.geneNum.merge.stat[sample.exp.list, "CDScount"])
	null_GH_gene = sum(data.merge.stat.dcast[sample.null.list, tmpGH])
	exp_GH_gene = sum(data.merge.stat.dcast[sample.exp.list, tmpGH])
	testor_species <- rbind(c(exp_GH, exp_species - exp_GH), c(null_GH, null_species - null_GH))
    model_species <- fisher.test(testor_species, alternative="two.sided")
    testor_gene <- rbind(c(exp_GH_gene, exp_gene - exp_GH_gene), c(null_GH_gene, null_gene - null_GH_gene))
    model_gene <- fisher.test(testor_gene, alternative="two.sided") 
	GH.enrich.info[i,] <- c(tmpGH, null_species, exp_species, null_GH, exp_GH, 
						null_gene, exp_gene, null_GH_gene, exp_GH_gene, 
						model_species$p.value, model_species$estimate,
                       	model_gene$p.value, model_gene$estimate)
}
GH.enrich.info$pvalue_species <- as.numeric(GH.enrich.info$pvalue_species)
GH.enrich.info$odds_species <- as.numeric(GH.enrich.info$odds_species)
GH.enrich.info$pvalue_gene <- as.numeric(GH.enrich.info$pvalue_gene)
GH.enrich.info$odds_gene <- as.numeric(GH.enrich.info$odds_gene)
GH.enrich.info$null_GH <- as.numeric(GH.enrich.info$null_GH)
GH.enrich.info$exp_GH <- as.numeric(GH.enrich.info$exp_GH)
GH.enrich.info$null_species <- as.numeric(GH.enrich.info$null_species)
GH.enrich.info$exp_species <- as.numeric(GH.enrich.info$exp_species)
GH.enrich.info$null_GH_gene <- as.numeric(GH.enrich.info$null_GH_gene)
GH.enrich.info$exp_GH_gene <- as.numeric(GH.enrich.info$exp_GH_gene)
GH.enrich.info$null_gene <- as.numeric(GH.enrich.info$null_gene)
GH.enrich.info$exp_gene <- as.numeric(GH.enrich.info$exp_gene)
GH.enrich.info$padj_species <- p.adjust(GH.enrich.info$pvalue_species, method = "BH")
GH.enrich.info$padj_gene <- p.adjust(GH.enrich.info$pvalue_gene, method = "BH")

Pvalue_threshold <- 0.05
GH.enrich.info$sig_gene <- "Not Sig"
GH.enrich.info[which(GH.enrich.info$padj_gene < Pvalue_threshold), "sig_gene"] <- "Sig"
GH.enrich.info$sig_species <- "Not Sig"
GH.enrich.info[which(GH.enrich.info$padj_species < Pvalue_threshold), "sig_species"] <- "Sig"

plotVector <- vector("list", length = 2)

p1 <- ggplot(GH.enrich.info, aes(x = log2(odds_gene), y = -log10(padj_gene))) +
    geom_point(aes(color = sig_gene), alpha = 0.3, size = 1) +
    geom_text_repel(data = GH.enrich.info[which(GH.enrich.info$padj_gene != 1),], aes(label = GH), size = 3, color = "#a50f15") +
    scale_color_manual(values=c("black", "red")) +
    labs(title = "GH enrichment on gene (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(GH.enrich.info, aes(x = log2(odds_species), y = -log10(padj_species))) +
    geom_point(aes(color = sig_species), alpha = 0.3, size = 1) +
    scale_color_manual(values=c("black", "red")) +
    labs(title = "GH enrichment on assembly (Envigo vs. Others)", x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))
plotVector[[1]] <- p1
plotVector[[2]] <- p2

#pdf("./figures/02.CAZy_GH_analysis.individualTest.pdf", width = 10, height = 4)
plot_grid(plotlist = plotVector, ncols = 2)
dev.off()



####################################


### GH cluster analysis ###

cluster.list <- c("CAZy_GH_1", "CAZy_GH_2", "CAZy_GH_3")
data.cluster.meta.stat <- data.frame(data.cluster.meta %>%
									group_by(cluster, supplier) %>%
									summarise(assemblyCount = n(), abundanceSum = sum(abundance)),
									stringsAsFactors = F)
data.cluster.meta.stat.supplier <- data.frame(data.cluster.meta %>%
									group_by(supplier) %>%
									summarise(abundanceSum = sum(abundance)),
									stringsAsFactors = F)
rownames(data.cluster.meta.stat.supplier) <- data.cluster.meta.stat.supplier$supplier
data.cluster.meta.stat$cluster <- factor(data.cluster.meta.stat$cluster, levels = cluster.list)
data.cluster.meta.stat$totalAb <- data.cluster.meta.stat.supplier[data.cluster.meta.stat$supplier,"abundanceSum"]
data.cluster.meta.stat$abundancePercent <- data.cluster.meta.stat$abundanceSum / data.cluster.meta.stat$totalAb
data.cluster.meta.stat <- plyr::ddply(data.cluster.meta.stat, "supplier",
                   						transform, label_ypos_count = cumsum(assemblyCount))
data.cluster.meta.stat <- plyr::ddply(data.cluster.meta.stat, "supplier",
                                        transform, label_ypos_abundance = cumsum(abundanceSum))

#pdf("./figures/02.CAZy_GH_analysis.cluster.species.pdf", width = 5, height = 5)
ggplot(data = data.cluster.meta.stat, aes(x = supplier, y = assemblyCount, fill = factor(cluster, levels = cluster.list))) +
    geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE)) +
	geom_text(aes(y = label_ypos_count, label = assemblyCount), vjust = 1.3, color = "black", size = 2.5) +
    theme_classic() +
    scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
    theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(y = "Number of species/MAGs", x = "Supplier", fill = "GH cluster", title = "Family: Muribaculaceae")
dev.off()


#pdf("./figures/02.CAZy_GH_analysis.cluster.abundance.pdf", width = 5, height = 5)
ggplot(data = data.cluster.meta.stat, aes(x = supplier, y = abundanceSum, fill = factor(cluster, levels = cluster.list))) +
    geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE)) + 
	geom_text(aes(y = label_ypos_abundance, label = paste(round(abundancePercent, 4) * 100, "%", sep = "")), vjust = 1.3, color = "black", size = 2.5) +
    theme_classic() +
    scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) + 
    theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 8)) +
    labs(y = "Accumulated absolute abundance", x = "Supplier", fill = "GH cluster", title = "Family: Muribaculaceae")
dev.off()


data.tree <- treeio::read.newick("./assembly.tree.newick")
p <- ggtree(data.tree, size = 0.5) +
        geom_tiplab(size = 2) +
        theme_tree2()

data.cluster.panel <- data.frame(id = data.tree$tip.label,
                                cluster = data.cluster.meta[data.tree$tip.label,"cluster"],
                                stringsAsFactors = F)
data.cluster.matrix <- data.frame(row.names = data.cluster.panel$id, cluster = data.cluster.panel[,c("cluster")])

#pdf("./figures/02.CAZy_GH_analysis.cluster.tree.pdf", width = 6, height = 7)
ggtree::gheatmap(p, data.cluster.matrix, offset = 0.1, width = 0.1, font.size = 3) +
    scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb"))
dev.off()
###########################















cluster.list <- c("CAZy_GH_1", "CAZy_GH_2", "CAZy_GH_3")
plotVector <- vector("list", length = 6)
for (j in 1:3){
	tmpCluster <- cluster.list[j]
	sample.null.list <- data.cluster.meta[data.cluster.meta$cluster != tmpCluster, "assembly"]
	sample.exp.list <- data.cluster.meta[data.cluster.meta$cluster == tmpCluster, "assembly"]

	GH.enrich.info <- data.frame(GH = rep(" ", ncol(data.merge.stat.dcast)),
						null_species = rep(0, ncol(data.merge.stat.dcast)),
						exp_species = rep(0, ncol(data.merge.stat.dcast)),
						null_GH = rep(0, ncol(data.merge.stat.dcast)),
						exp_GH = rep(0, ncol(data.merge.stat.dcast)),
						null_gene = rep(0, ncol(data.merge.stat.dcast)),
						exp_gene = rep(0, ncol(data.merge.stat.dcast)),
						null_GH_gene = rep(0, ncol(data.merge.stat.dcast)),
						exp_GH_gene = rep(0, ncol(data.merge.stat.dcast)),
						pvalue_species = rep(0, ncol(data.merge.stat.dcast)),
						odds_species = rep(0, ncol(data.merge.stat.dcast)),
						pvalue_gene = rep(0, ncol(data.merge.stat.dcast)),
                        odds_gene = rep(0, ncol(data.merge.stat.dcast)),
						stringsAsFactors = F)                                          
	for (i in 1:ncol(data.merge.stat.dcast)){
		tmpGH = colnames(data.merge.stat.dcast)[i]
   		null_species = length(sample.null.list)
    	exp_species = length(sample.exp.list)
		total_assembly = rownames(data.merge.stat.dcast)[which(data.merge.stat.dcast[,tmpGH] > 0)]
		null_GH = length(total_assembly[which(total_assembly %in% sample.null.list)])
		exp_GH = length(total_assembly[which(total_assembly %in% sample.exp.list)])
		null_gene = sum(data.geneNum.merge.stat[sample.null.list, "CDScount"])
		exp_gene = sum(data.geneNum.merge.stat[sample.exp.list, "CDScount"])
		null_GH_gene = sum(data.merge.stat.dcast[sample.null.list, tmpGH])
		exp_GH_gene = sum(data.merge.stat.dcast[sample.exp.list, tmpGH])
		testor_species <- rbind(c(exp_GH, exp_species - exp_GH), c(null_GH, null_species - null_GH))
    	model_species <- fisher.test(testor_species, alternative="greater")
    	testor_gene <- rbind(c(exp_GH_gene, exp_gene - exp_GH_gene), c(null_GH_gene, null_gene - null_GH_gene))
    	model_gene <- fisher.test(testor_gene, alternative="greater") 
		GH.enrich.info[i,] <- c(tmpGH, null_species, exp_species, null_GH, exp_GH, 
						null_gene, exp_gene, null_GH_gene, exp_GH_gene, 
						model_species$p.value, model_species$estimate,
                       	model_gene$p.value, model_gene$estimate)
	}
	GH.enrich.info$pvalue_species <- as.numeric(GH.enrich.info$pvalue_species)
	GH.enrich.info$odds_species <- as.numeric(GH.enrich.info$odds_species)
	GH.enrich.info$pvalue_gene <- as.numeric(GH.enrich.info$pvalue_gene)
	GH.enrich.info$odds_gene <- as.numeric(GH.enrich.info$odds_gene)
	GH.enrich.info$padj_species <- p.adjust(GH.enrich.info$pvalue_species, method = "BH")
	GH.enrich.info$padj_gene <- p.adjust(GH.enrich.info$pvalue_gene, method = "BH")

	Pvalue_threshold <- 0.05
	GH.enrich.info$sig_gene <- "Not Sig"
	GH.enrich.info[which(GH.enrich.info$padj_gene < Pvalue_threshold), "sig_gene"] <- "Sig"
	GH.enrich.info$sig_species <- "Not Sig"
	GH.enrich.info[which(GH.enrich.info$padj_species < Pvalue_threshold), "sig_species"] <- "Sig"

	p1 <- ggplot(GH.enrich.info, aes(x = log2(odds_gene), y = -log10(padj_gene))) +
    geom_point(aes(color = sig_gene), alpha = 0.3, size = 1) +
    geom_text_repel(data = GH.enrich.info[which(GH.enrich.info$sig_gene == "Sig"),], aes(label = GH), size = 3, color = "#a50f15") +
    scale_color_manual(values=c("black", "red")) +
    labs(title = paste("GH enrichment on gene (", tmpCluster ," vs. Others)", sep = ""), x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))
	p2 <- ggplot(GH.enrich.info, aes(x = log2(odds_species), y = -log10(padj_species))) +
    geom_point(aes(color = sig_species), alpha = 0.3, size = 1) +
	geom_text_repel(data = GH.enrich.info[which(GH.enrich.info$sig_species == "Sig"),], aes(label = GH), size = 3, color = "#a50f15") +
    scale_color_manual(values=c("black", "red")) +
    labs(title = paste("GH enrichment on assembly (", tmpCluster," vs. Others)", sep = ""), x = "log2OddsRatio", y = "-log10(adjusted P-value)") +
    theme(legend.position="none") +
    xlim(-5.8, 6.2) + ylim(0, 5) +
    theme(plot.title = element_text(hjust = 0.5))
	plotVector[[(j - 1) * 2 + 1]] <- p1
	plotVector[[(j - 1) * 2 + 2]] <- p2
}

#pdf("./figures/02.CAZy_GH_analysis.cluster.individualTest.pdf", width = 15, height = 8)
plot_grid(plotlist = plotVector, ncols = 3)
dev.off()
############################


assembly.cluster <- read.table("./assembly.CAZy_GH.cluster.tsv", stringsAsFactors = F, sep = "\t", header = F)
colnames(assembly.cluster) <- c("assembly", "cluster")
rownames(assembly.cluster) <- assembly.cluster$assembly

data.cluster.meta <- data.frame(row.names = assembly.cluster$assembly,
                                assembly = assembly.cluster$assembly,
                                cluster = assembly.cluster$cluster,
                                supplier = assembly.metadata.df[assembly.cluster$assembly,"supplier"],
                                abundance = assembly.metadata.df[assembly.cluster$assembly,"abundance"],
                                stringsAsFactors = F)

data.cluster.meta$geneWithGH <- data.geneNum.merge.stat[data.cluster.meta$assembly,"CDScount"]
data.cluster.meta$totalGHannotated <- rowSums(data.merge.stat.dcast[data.cluster.meta$assembly,])
data.cluster.meta$uniqueGHannotated <- rowSums(as.matrix(data.merge.stat.dcast[data.cluster.meta$assembly,] > 0) + 0)

#pdf("./figures/02.CAZy_GH_analysis.cluster.CDScount.pdf", width = 6, height = 4)
ggplot(data.cluster.meta, aes(x = supplier, y = geneWithGH, fill = cluster)) + 
	geom_boxplot(position = position_dodge(0.9)) + 
	theme_classic() + ylim(0, 500) + 
	scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
	labs(x = "Suppliers", y = "CDS with GH annotation") + 
    theme(axis.text.x = element_text(color="black", size = 8))
dev.off()

#pdf("./figures/02.CAZy_GH_analysis.cluster.totalGHcount.pdf", width = 6, height = 4)
ggplot(data.cluster.meta, aes(x = supplier, y = totalGHannotated, fill = cluster)) + 
	geom_boxplot(position = position_dodge(0.9)) + 
	theme_classic() + ylim(0, 600) + 
	scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
	labs(x = "Suppliers", y = "Number of total GH annotation") + 
    theme(axis.text.x = element_text(color="black", size = 8))
dev.off()

#pdf("./figures/02.CAZy_GH_analysis.cluster.uniqueGHcount.pdf", width = 6, height = 4)
ggplot(data.cluster.meta, aes(x = supplier, y = uniqueGHannotated, fill = cluster)) + 
	geom_boxplot(position = position_dodge(0.9)) + 
	theme_classic() + ylim(20, 80) + 
	scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
	labs(x = "Suppliers", y = "Number of total GH annotation") + 
    theme(axis.text.x = element_text(color="black", size = 8))
dev.off()

supplier.list <- c("AZC3J", "AZC5T", "AZC7E", "AZC9C")
cluster.list <- c("CAZy_GH_1", "CAZy_GH_2", "CAZy_GH_3")
data.GHcovered.info <- data.frame(supplier = rep("", length(supplier.list) * length(cluster.list)),
								cluster = rep("", length(supplier.list) * length(cluster.list)),
								GHcovered = rep("", length(supplier.list) * length(cluster.list)),
								stringsAsFactors = F)
for (i in 1:length(supplier.list)){
	for (j in 1:length(cluster.list)){
		tmpSupplier <- supplier.list[i]
		tmpCluster <- cluster.list[j]
		tmpAssembly.list <- data.cluster.meta[which(data.cluster.meta$supplier == tmpSupplier & data.cluster.meta$cluster == tmpCluster),"assembly"]
		tmpvalue <- sum(as.matrix(colSums(data.merge.stat.dcast[tmpAssembly.list,]) > 0) + 0)
		data.GHcovered.info[(i - 1) * length(cluster.list) + j,] <- c(tmpSupplier,  tmpCluster, tmpvalue)
	}
}
data.GHcovered.info$GHcovered <- as.numeric(data.GHcovered.info$GHcovered)
#pdf("./figures/02.CAZy_GH_analysis.cluster.uniqueGHcovered.pdf", width = 6, height = 3.6)
ggplot(data.GHcovered.info, aes(x = supplier, y = GHcovered - 25, fill = cluster)) +
	geom_bar(stat="identity", color="black", position =position_dodge(), width = 0.8) + 
    theme_classic() + ylim(0, 100) +
    scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
    labs(x = "Suppliers", y = "Number of unique GHs in MAGs") +
	theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))	
dev.off()









