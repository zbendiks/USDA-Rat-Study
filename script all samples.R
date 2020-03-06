#Zachary Bendiks
#Feb 2020
#USDA 1 RNA-Seq analysis - all samples

#annotations downloaded from MG-RAST

rm(list = ls())

#packages 
library(edgeR)
library(ggplot2)
library(vegan)
library(reshape2)
library(gplots)
library(ggdendro)
library(dendextend)
library(magrittr)
#1. Reading in and formatting RefSeq table####

USDA_table_refseq <- read.table("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/RefSeq all samples.tsv",
                                stringsAsFactors = FALSE,
                                header = TRUE,
                                sep = "\t")
#clean up names 
names(USDA_table_refseq) <- gsub("Rat.Cecal.Metatranscriptome.", "", names(USDA_table_refseq))
names(USDA_table_refseq)<- gsub( "X68_out.extendedFrags.fastq", "68", names(USDA_table_refseq))
#removing non-bacteria reads
USDA_table_refseq_bacteria <- USDA_table_refseq[USDA_table_refseq$domain == "Bacteria",]
#reordering column names so CON, WG, and RS samples are together
CON_sample_names <- c("68", "26", "58", "90")
WG_sample_names <- c("16", "24", "32", "40", "48")
RS_sample_names <- c("12", "20", "60","18", "76", "84")

USDA_table_refseq_bacteria <- USDA_table_refseq_bacteria[, c(names(USDA_table_refseq_bacteria)[1:6], 
                                                             CON_sample_names,
                                                             WG_sample_names,
                                                             RS_sample_names)]

#2. edgeR normalization, differential expression of RefSeq annotations at genus level, and volcano plots####
group = factor(c(rep("CON", 4), rep("WG", 5), rep("RS", 6)),
                    levels = c("CON", "WG", "RS"))

#genus level
USDA_table_refseq_genus <- USDA_table_refseq_bacteria[,-c(1:5)]
rownames(USDA_table_refseq_genus) <- USDA_table_refseq_bacteria[,"genus"]
USDA_table_refseq_genus <- USDA_table_refseq_genus[,-1]
colSums(USDA_table_refseq_genus)
#View(USDA_table_refseq_genus)
rnaseqMatrix_refseq_genus = round(USDA_table_refseq_genus) 
rnaseqMatrix_refseq_genus = rnaseqMatrix_refseq_genus[rowSums(cpm(rnaseqMatrix_refseq_genus) > 1) >= 2,]
nrow(rnaseqMatrix_refseq_genus) #582 genera after filtering out low abundance features

#set up DGE 
exp_study_refseq_genus = DGEList(counts=rnaseqMatrix_refseq_genus, group=group)
#calculate TMM normalization factors for library size adjustment
exp_study_refseq_genus = calcNormFactors(exp_study_refseq_genus, method = "TMM")
genus_refseq_TMM <- cpm.DGEList(exp_study_refseq_genus)
#calculate the average and genewise dispersion values, then plot 
exp_study_refseq_genus = estimateCommonDisp(exp_study_refseq_genus)
exp_study_refseq_genus = estimateTagwiseDisp(exp_study_refseq_genus)
plotBCV(exp_study_refseq_genus) 
#slight trend but not vert robust, won't fit a linear model to it 

#set up model matrix using the group variable and use to re-estimate dispersion
design <- model.matrix(~group)
colnames(design) <- c("(Intercept)", "WG", "RS")
exp_study_refseq_genus <- estimateDisp(exp_study_refseq_genus, design)

#make contrasts to compare each exp groups to the control
my.contrasts <- makeContrasts(CON_WG = CON-WG,
                              CON_RS = CON-RS,
                              levels = design)

#perform quasi-likelihood F-test which teh edgeR devs recommend for bulk RNA-Seq analysis 
fit <- glmQLFit(exp_study_refseq_genus, design)
#comparing the CON and WG groups with coef = 2
qlf_CON_WG <- glmQLFTest(fit, coef = 2) 
qlf_CON_WG_toptags <- topTags(qlf_CON_WG, n = nrow(qlf_CON_WG$table))
nrow(qlf_CON_WG_toptags$table) #582
sum(qlf_CON_WG_toptags$table$FDR < 0.05) #0 of 582 DE 
sum(qlf_CON_WG_toptags$table$FDR > 0.05) #582 n.d.
qlf_CON_WG_toptags$table[rownames(qlf_CON_WG_toptags$table) == "Lactobacillus",]
#comparing the CON and RS groups with coef = 3
qlf_CON_RS <- glmQLFTest(fit, coef = 3) 
qlf_CON_RS_toptags <- topTags(qlf_CON_RS, n = nrow(qlf_CON_RS$table))

nrow(qlf_CON_RS_toptags) #582
sum(qlf_CON_RS_toptags$table$FDR < 0.05) #234 of 582 DE 
sum(qlf_CON_RS_toptags$table$FDR < 0.05 & qlf_CON_RS_toptags$table$logFC > 0) #120 of 234  upregulated
sum(qlf_CON_RS_toptags$table$FDR < 0.05 & qlf_CON_RS_toptags$table$logFC < 0) #114 of 234  upregulated
sum(qlf_CON_RS_toptags$table$FDR > 0.05) #348 n.d.
qlf_CON_RS_toptags$table[rownames(qlf_CON_RS_toptags$table) == "Lactobacillus",]
qlf_CON_RS_toptags$table[rownames(qlf_CON_RS_toptags$table) == "Bacteroides",]
qlf_CON_RS_toptags$table[rownames(qlf_CON_RS_toptags$table) == "Parabacteroides",]
qlf_CON_RS_toptags$table[rownames(qlf_CON_RS_toptags$table) == "Prevotella",]
qlf_CON_RS_toptags$table[rownames(qlf_CON_RS_toptags$table) == "Porphyromonas",]
qlf_CON_RS_toptags$table[grep("Desulf" , rownames(qlf_CON_RS_toptags$table)),]

#save edgeR results to files 
setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
write.csv(qlf_CON_WG$table, file = "CON WG RefSeq DE genera.csv", row.names = TRUE)
write.csv(qlf_CON_RS$table, file = "CON RS RefSeq DE genera.csv", row.names = TRUE)

#adding full taxa names to RS table 
#subset refseq table to just include the taxa included for DE analysis
USDA_table_refseq_bacteria_toptags <- USDA_table_refseq_bacteria[USDA_table_refseq_bacteria$genus %in% rownames(qlf_CON_RS_toptags$table),]
nrow(USDA_table_refseq_bacteria_toptags) #582 functions as expected
#order the genus annotations alphabetically in both objects 
USDA_table_refseq_bacteria_toptags <- USDA_table_refseq_bacteria_toptags[order(USDA_table_refseq_bacteria_toptags$genus),]
USDA_result_table_refseq <- qlf_CON_RS_toptags$table[order(rownames(qlf_CON_RS_toptags$table)),]
nrow(USDA_result_table_refseq) #582 functions
sum(rownames(USDA_result_table_refseq) == USDA_table_refseq_bacteria_toptags$genus)
#add the full strings to result table, write to CSV
USDA_result_table_refseq_full_string <- cbind(USDA_table_refseq_bacteria_toptags[,1:5],
                                              rownames(USDA_result_table_refseq),
                                              USDA_result_table_refseq)
View(USDA_result_table_refseq_full_string)
names(USDA_result_table_refseq_full_string)[6] <- "genus"
write.csv(USDA_result_table_refseq_full_string, file = "CON RS RefSeq DE genera.csv", row.names = FALSE)


#volcano plots of WG and RS comparisons
setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("CON WG RefSeq Volcano.tiff", width = 2, height = 1.5, units = "in", res = 1200)
par(mgp=c(0,0.1,0), mar=c(2,1.5,0.1,0.1))
par(oma=c(0,0,0,0))

plot(qlf_CON_WG_toptags$table$logFC, 
     -1*log10(qlf_CON_WG_toptags$table$FDR), 
     xlab=NA, 
     ylab=NA,
     pch=16,
     xaxt = "none", yaxt = "none",
     col=ifelse(qlf_CON_WG_toptags$table$FDR<=0.1 & qlf_CON_WG_toptags$table$logFC < 0, "gray50",
                ifelse(qlf_CON_WG_toptags$table$FDR<=0.1 & qlf_CON_WG_toptags$table$logFC > 0, "darkorange", "black")),
     cex = 0.5,
     ylim = c(0,1.65))

abline(-1*log10(0.05),0, lty = 2, lwd = 1, col = "grey40") #add line at FDR = 0.05 
abline(-1*log10(0.1),0, lty = 2, lwd = 1, col = "grey40") #add line at FDR = 0.1 
abline(v = 0, lty = 2, lwd = 1, col = "grey40") #add line at logFC = 0

legend(x = -0.75, 
       y =  1.85, 
       legend = "   CON \nn=0 (1)", 
       text.col = "gray50",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)
legend(x = 0.25, 
       y = 1.85, 
       legend = "WG \nn=0", 
       text.col = "darkorange",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)
legend(x = -2.15, 
       y = 0.45, 
       legend = "N.D. \nn=582", 
       text.col = "black",
       bty = "n",
       cex = 0.60,
       xjust = 0.5)
legend(x = 1.65, 
       y = -1*log10(0.05) + 0.2, 
       legend = "p = 0.05", 
       text.col = "red",
       bty = "n",
       cex = 0.4,
       xjust = 0.5,
       text.font = 3)
legend(x = 1.65, 
       y = -1*log10(0.1) + 0.2, 
       legend = "p = 0.1  ", 
       text.col = "red",
       bty = "n",
       cex = 0.4,
       xjust = 0.5,
       text.font = 3)

axis(side = 1, cex.axis = 0.6, tck=-.025)
axis(side = 2, cex.axis = 0.6, tck=-.025)

mtext(side = 1, line = 0.75, "log(Fold Change)", cex = 0.75) #line chances label distance from plot
mtext(side = 2, line = 0.75, "-log(FDR)", cex = 0.75)

dev.off()


tiff("CON RS RefSeq Volcano.tiff", width = 2, height = 1.5, units = "in", res = 1200)
par(mgp=c(0,0.1,0), mar=c(2,1.5,0.1,0.1))
par(oma=c(0,0,0,0))

plot(qlf_CON_RS_toptags$table$logFC, 
     -1*log10(qlf_CON_RS_toptags$table$FDR), 
     xlab=NA, 
     ylab=NA,
     pch=16,
     xaxt = "none", yaxt = "none",
     col=ifelse(qlf_CON_RS_toptags$table$FDR<=0.05 & qlf_CON_RS_toptags$table$logFC < 0, "gray50",
                ifelse(qlf_CON_RS_toptags$table$FDR<=0.05 & qlf_CON_RS_toptags$table$logFC > 0, "blue", "black")),
     cex = 0.5,
     ylim = c(0,7.25),
     xlim = c(-4.25, 5.25))


abline(-1*log10(0.05),0, lty = 2, lwd = 1, col = "grey40") #add line at FDR = 0.05 
abline(v = 0, lty = 2, lwd = 1, col = "grey40") #add line at logFC = 0

legend(x = -1.65, 
       y =  8.1, 
       legend = "  CON \nn=114", 
       text.col = "gray50",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

legend(x = 0.75, 
       y = 8.1, 
       legend = "RS \nn=120", 
       text.col = "blue",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

legend(x = -3.15, 
       y = 1.5, 
       legend = "N.D. n=348", 
       text.col = "black",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

axis(side = 1, cex.axis = 0.6, tck=-.025)
axis(side = 2, cex.axis = 0.6, tck=-.025)

mtext(side = 1, line = 0.75, "log(Fold Change)", cex = 0.75)
mtext(side = 2, line = 0.75, "-log(FDR)", cex = 0.75)

dev.off()


#3. Condense RefSeq annotations at order level, normalize, and plot####
#order level - condense rows by taxa, normalize, and make figure
length(unique(USDA_table_refseq_bacteria$genus)) #593
length(unique(USDA_table_refseq_bacteria$family)) #243
length(unique(USDA_table_refseq_bacteria$order)) #110

USDA_table_refseq_order <- USDA_table_refseq_bacteria[,-c(1,2,3,5,6)]
USDA_table_refseq_order <- aggregate(. ~ USDA_table_refseq_order$order, 
          data = USDA_table_refseq_order[2:ncol(USDA_table_refseq_order)], 
          FUN = "sum")
rownames(USDA_table_refseq_order) <- USDA_table_refseq_order$`USDA_table_refseq_order$order`
USDA_table_refseq_order <- USDA_table_refseq_order[,-1]
#nrow(USDA_table_refseq_order) #110 as expected

#edgeR normalization of order level counts 
rnaseqMatrix_refseq_order = round(USDA_table_refseq_order) 
rnaseqMatrix_refseq_order = rnaseqMatrix_refseq_order[rowSums(cpm(rnaseqMatrix_refseq_order) > 1) >= 2,]
nrow(rnaseqMatrix_refseq_order) #110 still

#set up DGE 
exp_study_refseq_order = DGEList(counts=rnaseqMatrix_refseq_order, group=group)
#calculate TMM normalization factors for library size adjustment
exp_study_refseq_order = calcNormFactors(exp_study_refseq_order, method = "TMM")
order_refseq_TMM <- cpm.DGEList(exp_study_refseq_order)

refseq_order_sums <- data.frame(rowSums(order_refseq_TMM))
refseq_order_sums_ordered <- order(refseq_order_sums,
                           refseq_order_sums$rowSums.order_refseq_TMM.)
refseq_order_sums <- refseq_order_sums[refseq_order_sums_ordered, ,drop = FALSE]

#setting 100,000 TMM counts as cutoff, aggregating all other orders within 'Other'. Will give 9 taxa + aggregated group 
high_abun_orders <- row.names(
  refseq_order_sums[refseq_order_sums$rowSums.order_refseq_TMM. > 100000, , drop = FALSE]
  )
low_abun_orders <- row.names(
  refseq_order_sums[refseq_order_sums$rowSums.order_refseq_TMM. < 100000, , drop = FALSE]
)

refseq_order_condensed <- order_refseq_TMM[row.names(order_refseq_TMM) %in% high_abun_orders, , drop = FALSE]
refseq_order_condensed <- refseq_order_condensed[order(rowSums(refseq_order_condensed), decreasing = TRUE), ] 

Other <- order_refseq_TMM[row.names(order_refseq_TMM) %in% low_abun_orders, , drop = FALSE]
Other <- colSums(Other)

refseq_order_condensed <- rbind(refseq_order_condensed, Other)
refseq_order_condensed_melt <- melt(refseq_order_condensed)
names(refseq_order_condensed_melt) <- c("RefSeq Order", "Sample", "TMM")

#set orders as factors, specify order
refseq_order_condensed_melt$`RefSeq Order` <- factor(refseq_order_condensed_melt$`RefSeq Order`,
       levels = c("Clostridiales",
                  "Bacteroidales",
                  "Erysipelotrichales",
                  "Lactobacillales",
                  "Bacillales",
                  "Bifidobacteriales",
                  "Coriobacteriales",
                  "Thermoanaerobacterales",
                  "Burkholderiales",
                  "Other"))

refseq_order_condensed_melt$Sample <- factor(refseq_order_condensed_melt$Sample,
                                                     levels = unique(refseq_order_condensed_melt$Sample))


Refseq_order_all_samples_TMM_abundance <- ggplot(data=refseq_order_condensed_melt, 
       aes(x=Sample, y=TMM, fill=`RefSeq Order`)) +
  geom_rect(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 4.5, fill = "gray", alpha = 0.5)+
  geom_rect(ymin = -Inf, ymax = Inf, xmin = 4.5, xmax = 9.5, fill = "orange", alpha = 0.5)+
  geom_rect(ymin = -Inf, ymax = Inf, xmin = 9.5, xmax = Inf, fill = "lightblue", alpha = 0.5)+
  geom_bar(stat="identity", color = "black", width = 1) +
  scale_fill_manual(values = rainbow(10))+
  theme_bw()+
  annotate("text", x = 2.5, y = 2000000, label = "bold(CON)", parse = TRUE, size = 2.5)+
  annotate("text", x = 7, y = 2000000, label = "bold(WG)", parse = TRUE, size = 2.5)+
  annotate("text", x = 12.5, y = 2000000, label = "bold(RS)", parse = TRUE, size = 2.5)+
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.key.height = unit(0.75, "line"),
        legend.key.width = unit(1, "line"))
  
setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("USDA_Refseq_order_all_samples_TMM_abundance.tiff", width = 4.5, height = 3.0, units = "in", compression = "lzw", res = 1200)
Refseq_order_all_samples_TMM_abundance
dev.off()

#4. Condense RefSeq annotations at family level, normalize, and plot####
#family level - condense rows by taxa, normalize, and make figure
length(unique(USDA_table_refseq_bacteria$genus)) #593
length(unique(USDA_table_refseq_bacteria$family)) #243
length(unique(USDA_table_refseq_bacteria$order)) #110

USDA_table_refseq_family <- USDA_table_refseq_bacteria[,-c(1,2,3,4,6)]
USDA_table_refseq_family <- aggregate(. ~ USDA_table_refseq_family$family, 
                                     data = USDA_table_refseq_family[2:ncol(USDA_table_refseq_family)], 
                                     FUN = "sum")
rownames(USDA_table_refseq_family) <- USDA_table_refseq_family$`USDA_table_refseq_family$family`
USDA_table_refseq_family <- USDA_table_refseq_family[,-1]
#nrow(USDA_table_refseq_family) #243 as expected

#edgeR normalization of family level counts 
rnaseqMatrix_refseq_family = round(USDA_table_refseq_family) 
rnaseqMatrix_refseq_family = rnaseqMatrix_refseq_family[rowSums(cpm(rnaseqMatrix_refseq_family) > 1) >= 2,]
nrow(rnaseqMatrix_refseq_family) #240

#set up DGE 
exp_study_refseq_family = DGEList(counts=rnaseqMatrix_refseq_family, group=group)
#calculate TMM normalization factors for library size adjustment
exp_study_refseq_family = calcNormFactors(exp_study_refseq_family, method = "TMM")
family_refseq_TMM <- cpm.DGEList(exp_study_refseq_family)

refseq_family_sums <- data.frame(rowSums(family_refseq_TMM))
refseq_family_sums_ordered <- order(refseq_family_sums,
                                   refseq_family_sums$rowSums.family_refseq_TMM.)
refseq_family_sums <- refseq_family_sums[refseq_family_sums_ordered, ,drop = FALSE]

#setting 300,000 TMM counts as cutoff, aggregating all other familys within 'Other'. Will give 9 taxa + aggregated group 
high_abun_familys <- row.names(
     refseq_family_sums[refseq_family_sums$rowSums.family_refseq_TMM. > 300000, , drop = FALSE]
)
low_abun_familys <- row.names(
     refseq_family_sums[refseq_family_sums$rowSums.family_refseq_TMM. < 300000, , drop = FALSE]
)

refseq_family_condensed <- family_refseq_TMM[row.names(family_refseq_TMM) %in% high_abun_familys, , drop = FALSE]
refseq_family_condensed <- refseq_family_condensed[order(rowSums(refseq_family_condensed), decreasing = TRUE), ] 

Other <- family_refseq_TMM[row.names(family_refseq_TMM) %in% low_abun_familys, , drop = FALSE]
Other <- colSums(Other)

refseq_family_condensed <- rbind(refseq_family_condensed, Other)
refseq_family_condensed_melt <- melt(refseq_family_condensed)
names(refseq_family_condensed_melt) <- c("RefSeq family", "Sample", "TMM")

#set familys as factors, specify family
refseq_family_condensed_melt$`RefSeq family` <- factor(refseq_family_condensed_melt$`RefSeq family`,
                                                     levels = unique(refseq_family_condensed_melt$`RefSeq family`))

refseq_family_condensed_melt$Sample <- factor(refseq_family_condensed_melt$Sample,
                                             levels = unique(refseq_family_condensed_melt$Sample))

col_palette <- rainbow(11)
col_palette[5] <- "chartreuse3"
col_palette[11] <- "black"

Refseq_family_all_samples_TMM_abundance <- ggplot(data=refseq_family_condensed_melt, 
                                                 aes(x=Sample, y=TMM, fill=`RefSeq family`)) +
     geom_rect(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 4.5, fill = "gray", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 4.5, xmax = 9.5, fill = "orange", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 9.5, xmax = Inf, fill = "lightblue", alpha = 0.5)+
     geom_bar(stat="identity", width = 1) +
     scale_fill_manual(values = col_palette)+
     theme_bw()+
     annotate("text", x = 2.5, y = 2500000, label = "bold(CON)", parse = TRUE, size = 3.5)+
     annotate("text", x = 7, y = 2500000, label = "bold(WG)", parse = TRUE, size = 3.5)+
     annotate("text", x = 12.5, y = 2500000, label = "bold(RS)", parse = TRUE, size = 3.5)+
     theme(axis.text.y = element_text(colour = "black"),
           axis.text.x = element_blank(),
           legend.text = element_text(size=8),
           legend.title = element_text(size=10),
           legend.key.height = unit(0.75, "line"),
           legend.key.width = unit(1, "line"))

setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("USDA_Refseq_family_all_samples_TMM_abundance.tiff", width = 4.5, height = 3.0, units = "in", compression = "lzw", res = 1200)
Refseq_family_all_samples_TMM_abundance
dev.off()

#remaking family plot but TMM proportions
#View(refseq_family_condensed)
refseq_family_condensed_proportion <- prop.table(refseq_family_condensed, 2)
refseq_family_condensed_proportion_melt <- melt(refseq_family_condensed_proportion)
names(refseq_family_condensed_proportion_melt) <- c("RefSeq family", "Sample", "TMM Proportion")

#set familys as factors, specify family
refseq_family_condensed_proportion_melt$`RefSeq family` <- factor(refseq_family_condensed_proportion_melt$`RefSeq family`,
                                                       levels = unique(refseq_family_condensed_proportion_melt$`RefSeq family`))

refseq_family_condensed_proportion_melt$Sample <- factor(refseq_family_condensed_proportion_melt$Sample,
                                              levels = unique(refseq_family_condensed_proportion_melt$Sample))

col_palette <- rainbow(11)
col_palette[5] <- "chartreuse3"
col_palette[11] <- "black"


Refseq_family_all_samples_TMM_proportion <- ggplot(data=refseq_family_condensed_proportion_melt, 
                                                  aes(x=Sample, y=`TMM Proportion`, fill=`RefSeq family`)) +
     geom_rect(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 4.5, fill = "gray", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 4.5, xmax = 9.5, fill = "orange", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 9.5, xmax = Inf, fill = "lightblue", alpha = 0.5)+
     geom_bar(stat="identity", width = 1) +
     scale_fill_manual(values = col_palette)+
     theme_bw()+
     annotate("text", x = 2.5, y = 1.05, label = "bold(CON)", parse = TRUE, size = 3.5)+
     annotate("text", x = 7, y = 1.05, label = "bold(WG)", parse = TRUE, size = 3.5)+
     annotate("text", x = 12.5, y = 1.05, label = "bold(RS)", parse = TRUE, size = 3.5)+
     theme(axis.text.y = element_text(colour = "black"),
           axis.text.x = element_blank(),
           legend.text = element_text(size=8),
           legend.title = element_text(size=10),
           legend.key.height = unit(0.75, "line"),
           legend.key.width = unit(1, "line"))

setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("USDA_Refseq_family_all_samples_TMM_proportion.tiff", width = 4.5, height = 3.0, units = "in", compression = "lzw", res = 1200)
Refseq_family_all_samples_TMM_proportion
dev.off()

#10. Plot log FC of highly expressed significant RefSeq genera ####
#sort by full string object by logFC
USDA_result_table_refseq_string <- USDA_result_table_refseq_full_string[order(USDA_result_table_refseq_full_string$logFC, decreasing = TRUE),]
#extract significant features
USDA_result_table_refseq_string_sig <- USDA_result_table_refseq_string[USDA_result_table_refseq_string$FDR < 0.05,]
#sort from highest to least normalized counts
View(head(USDA_result_table_refseq_string_sig[order(USDA_result_table_refseq_string_sig$logCPM, decreasing = TRUE),], n = 40))
#extract the 40 most abundant refseq genera
USDA_result_table_refseq_string_sig_top40 <- head(USDA_result_table_refseq_string_sig[order(USDA_result_table_refseq_string_sig$logCPM, decreasing = TRUE),], n = 40)

#make a column with taxa names as factors so ggplot plots taxa in the specified order (sorted by logFC)
USDA_result_table_refseq_string_sig_top40 <- USDA_result_table_refseq_string_sig_top40[order(USDA_result_table_refseq_string_sig_top40$logFC, decreasing = TRUE),]
USDA_result_table_refseq_string_sig_top40$genus <- factor(USDA_result_table_refseq_string_sig_top40$genus,
                                                              levels = unique(USDA_result_table_refseq_string_sig_top40$genus))

View(USDA_result_table_refseq_string_sig_top40)
unique(USDA_result_table_refseq_string_sig_top40$family) #25 unique family
unique(USDA_result_table_refseq_string_sig_top40$order) #13 unique orders

ggplot(data = USDA_result_table_refseq_string_sig_top40, aes(x = genus, y =logFC)) +
     geom_bar(stat = "identity", aes(fill = order))+
     #  theme(axis.text.y=element_text(color = result_table_fam_sig_only$color))+
     coord_flip()

#create a plot of the normalized counts for each gene (object already ordered properly so just plot the values)
counts_plot_by_refseq <- ggplot(data = USDA_result_table_refseq_string_sig_top40,
                             aes(x = logCPM, y = genus)) + 
     geom_point(size = 2.5, aes(color = order), show.legend = FALSE)+
     geom_segment(size = 1.25, aes(x = 9,
                                   xend = logCPM,
                                   y = genus, 
                                   yend = genus,
                                   color = order), show.legend = FALSE)+
     theme_classic() +
     scale_y_discrete(limits=rev(levels(USDA_result_table_refseq_string_sig_top40$genus))) +
     xlab("log\n(TMM)")+ylab("")+xlim(c(9,18))+
     theme(axis.text.y = element_blank(),
           axis.text.x = element_text(color = "black"),
           plot.margin = unit(c(0.1,0.1,0.175,0.0), "cm")) # , , bottom, 

refseq_barchart_function_logFC_top40 <- ggplot(data = USDA_result_table_refseq_string_sig_top40, aes(x = genus, y =logFC)) +
     geom_bar(stat = "identity", aes(fill = order), color = "black")+
     theme(axis.text.y=element_text(color = "black", size = 8), 
           axis.text.x=element_text(color = "black"),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_line(color = "gray80"),
           panel.background = element_rect(fill = "white", colour = "black"),
           legend.key.size = unit(0.3, "cm"))+
     #  theme_bw()+
     xlab('RefSeq Genera')+ylab('log(Fold Change)\n(RS/CON)')+
     coord_flip()+
     scale_x_discrete(limits=rev(levels(USDA_result_table_refseq_string_sig_top40$genus)))+
     guides(fill = guide_legend(title = "RefSeq Order", ncol = 3))+
     theme(legend.position = c(0.25,-0.25),
           plot.margin = unit(c(0.1,3.5,3.5,0.1), "cm"))+
     annotation_custom(ggplotGrob(counts_plot_by_refseq), xmin = -5.5, xmax = 41, ymin = 6, ymax = 9.5)

tiff("CON_RS.edgeR.refseq.taxa_barchart_logFC.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 600)
refseq_barchart_function_logFC_top40
dev.off()

pdf("CON_RS.edgeR.refseq.taxa_barchart_logFC.pdf", width = 6.5, height = 6)
refseq_barchart_function_logFC_top40
dev.off()


#heatmap of taxonomic changes####

heatmap.2(prop.table(genus_refseq_TMM, 2),
          scale = "row",
          ColSideColors = c(rep("gray",4), rep("orange", 5), rep("blue", 6)),
          trace = "none")

#5. Reading in and formatting COG table####
#had trouble reading in COG file (b/c spaces in sample names).  Converted to CSV to import properly 
USDA_table_COG <- read.csv("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/COGs all samples.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE)
View(USDA_table_COG)
names(USDA_table_COG) <- gsub("Rat.Cecal.Metatranscriptome.", "", names(USDA_table_COG))
names(USDA_table_COG)<- gsub( "X68_out.extendedFrags.fastq", "68", names(USDA_table_COG))

USDA_table_COG <- USDA_table_COG[, c(names(USDA_table_COG)[1:3], 
                                                             CON_sample_names,
                                                             WG_sample_names,
                                                             RS_sample_names)]

USDA_table_COG_fxn <- USDA_table_COG[,-c(1:2)]
rownames(USDA_table_COG_fxn) <- USDA_table_COG_fxn[,"function."]
USDA_table_COG_fxn <- USDA_table_COG_fxn[,-1]
colSums(USDA_table_COG_fxn)
#View(USDA_table_refseq_COG)


#6. edgeR normalization and differential expression of COGs annotations at function level and volcano plots####
rnaseqMatrix_COG = round(USDA_table_COG_fxn) 
rnaseqMatrix_COG = rnaseqMatrix_COG[rowSums(cpm(rnaseqMatrix_COG) > 1) >= 2,]
nrow(rnaseqMatrix_COG) #2627 functions after filtering out low abundance features

#set up DGE 
exp_study_COG = DGEList(counts=rnaseqMatrix_COG, group=group)
#calculate TMM normalization factors for library size adjustment
exp_study_COG = calcNormFactors(exp_study_COG, method = "TMM")
#calculate the average and genewise dispersion values, then plot 
exp_study_COG = estimateCommonDisp(exp_study_COG)
exp_study_COG = estimateTagwiseDisp(exp_study_COG)
#plotBCV(exp_study_COG) 

#clear negative trend between CV and CPM. I will fit a linear model to the residuals before doing DE analysis
exp_study_COG <- estimateGLMCommonDisp(exp_study_COG, design)
exp_study_COG <- estimateGLMTrendedDisp(exp_study_COG, design, method="power")
exp_study_COG <- estimateGLMTagwiseDisp(exp_study_COG, design)
#plotBCV(exp_study_COG) 
#GLM looks much better than the common linear model used initially, so fitting data to that linear model
fit_cog <- glmQLFit(exp_study_COG, design)
lrt_cog <- glmLRT(fit_cog) 
COG_TMM <- cpm.DGEList(exp_study_COG) #use the lrt_object to obtain normalized values

#comparing the CON and WG groups with coef = 2
fit_cog_CON_WG <- glmQLFTest(fit_cog, coef = 2) 
fit_cog_CON_WG_toptags <- topTags(fit_cog_CON_WG, n = nrow(fit_cog_CON_WG))
nrow(fit_cog_CON_WG_toptags$table) #2627
sum(fit_cog_CON_WG_toptags$table$FDR < 0.05) #0 DE COGs total
sum(fit_cog_CON_WG_toptags$table$FDR < 0.05 & fit_cog_CON_WG_toptags$table$logFC > 0) #0 upregulated COGs
sum(fit_cog_CON_WG_toptags$table$FDR < 0.05 & fit_cog_CON_WG_toptags$table$logFC < 0) #0 downregulated COGs
sum(fit_cog_CON_WG_toptags$table$FDR < 0.1 & fit_cog_CON_WG_toptags$table$logFC > 0) #10 near upregulated COGs
sum(fit_cog_CON_WG_toptags$table$FDR < 0.1 & fit_cog_CON_WG_toptags$table$logFC < 0) #181 near downregulated COGs
sum(fit_cog_CON_WG_toptags$table$FDR > 0.05) #2627 ND COGs 
View(fit_cog_CON_WG_toptags$table)
#comparing the CON and WG groups with coef = 2
fit_cog_CON_RS <- glmQLFTest(fit_cog, coef = 3) 
fit_cog_CON_RS_toptags <- topTags(fit_cog_CON_RS, n = nrow(fit_cog_CON_RS))
View(fit_cog_CON_RS_toptags$table)
nrow(fit_cog_CON_RS_toptags$table) #2627
sum(fit_cog_CON_RS_toptags$table$FDR < 0.05) #1218 DE COGs total
sum(fit_cog_CON_RS_toptags$table$FDR < 0.05 & fit_cog_CON_RS_toptags$table$logFC > 0) #483 upregulated COGs
sum(fit_cog_CON_RS_toptags$table$FDR < 0.05 & fit_cog_CON_RS_toptags$table$logFC < 0) #735 downregulated COGs
sum(fit_cog_CON_RS_toptags$table$FDR > 0.05) #1409 ND COGs 

#save edgeR results to files 
setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
write.csv(fit_cog_CON_WG_toptags$table, file = "CON WG COG DE functions.csv", row.names = TRUE)
write.csv(fit_cog_CON_RS_toptags$table, file = "CON RS COG DE functions.csv", row.names = TRUE)

#7. Volcano plots of WG and RS comparisons (COGs)####
setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("CON WG COG Volcano.tiff", width = 2, height = 1.5, units = "in", res = 1200)
par(mgp=c(0,0.1,0), mar=c(2,1.5,0.1,0.1))
par(oma=c(0,0,0,0))

plot(fit_cog_CON_WG_toptags$table$logFC, 
     -1*log10(fit_cog_CON_WG_toptags$table$FDR), 
     xlab=NA, 
     ylab=NA,
     pch=16,
     xaxt = "none", yaxt = "none",
     col=ifelse(fit_cog_CON_WG_toptags$table$FDR<=0.1 & fit_cog_CON_WG_toptags$table$logFC < 0, "gray50",
                ifelse(fit_cog_CON_WG_toptags$table$FDR<=0.1 & fit_cog_CON_WG_toptags$table$logFC > 0, "darkorange", "black")),
     cex = 0.5,
     ylim = c(0,1.65))

abline(-1*log10(0.05),0, lty = 2, lwd = 1, col = "grey40") #add line at FDR = 0.05 
abline(-1*log10(0.1),0, lty = 2, lwd = 1, col = "grey40")
abline(v = 0, lty = 2, lwd = 1, col = "grey40") #add line at logFC = 0

legend(x = -3, 
       y =  1.85, 
       legend = "       CON \nn=0 (181)", 
       text.col = "gray50",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)
legend(x = 1.5, 
       y = 1.85, 
       legend = "WG \nn=0 (10)", 
       text.col = "darkorange",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)
legend(x = -7.25, 
       y = 0.55, 
       legend = "N.D.\nn=2627", 
       text.col = "black",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)
legend(x = 5.5, 
       y = -1*log10(0.05) + 0.2, 
       legend = "p = 0.05", 
       text.col = "red",
       bty = "n",
       cex = 0.4,
       xjust = 0.5,
       text.font = 3)
legend(x = 5.5, 
       y = -1*log10(0.1) + 0.05, 
       legend = "p = 0.1  ", 
       text.col = "red",
       bty = "n",
       cex = 0.4,
       xjust = 0.5,
       text.font = 3)
axis(side = 1, cex.axis = 0.6, tck=-.025)
axis(side = 2, cex.axis = 0.6, tck=-.025)

mtext(side = 1, line = 0.75, "log(Fold Change)", cex = 0.75)
mtext(side = 2, line = 0.75, "-log(FDR)", cex = 0.75)

dev.off()


tiff("CON RS COG Volcano.tiff", width = 2, height = 1.5, units = "in", res = 1200)
par(mgp=c(0,0.1,0), mar=c(2,1.5,0.1,0.1))
par(oma=c(0,0,0,0))

plot(fit_cog_CON_RS_toptags$table$logFC, 
     -1*log10(fit_cog_CON_RS_toptags$table$FDR), 
     xlab=NA, 
     ylab=NA,
     pch=16,
     xaxt = "none", yaxt = "none",
     col=ifelse(fit_cog_CON_RS_toptags$table$FDR<=0.05 & fit_cog_CON_RS_toptags$table$logFC < 0, "gray50",
                ifelse(fit_cog_CON_RS_toptags$table$FDR<=0.05 & fit_cog_CON_RS_toptags$table$logFC > 0, "blue", "black")),
     cex = 0.5,
     ylim = c(0,6.45))


abline(-1*log10(0.05),0, lty = 2, lwd = 1, col = "grey40") #add line at FDR = 0.05 
abline(v = 0, lty = 2, lwd = 1, col = "grey40") #add line at logFC = 0

legend(x = -2.5, 
       y = 7.25, 
       legend = "  CON \nn=735", 
       text.col = "gray50",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

legend(x = 1.5, 
       y = 7.25, 
       legend = "RS \nn=483", 
       text.col = "blue",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

legend(x = -7.25, 
       y = 1.85, 
       legend = "N.D. \nn=1409", 
       text.col = "black",
       bty = "n",
       cex = 0.65,
       xjust = 0.5)

axis(side = 1, cex.axis = 0.6, tck=-.025)
axis(side = 2, cex.axis = 0.6, tck=-.025)

mtext(side = 1, line = 0.75, "log(Fold Change)", cex = 0.75)
mtext(side = 2, line = 0.75, "-log(FDR)", cex = 0.75)

dev.off()

#8. Barplot of COG annotations####
length(unique(USDA_table_COG$function.)) #2932
length(unique(USDA_table_COG$level2)) #25
length(unique(USDA_table_COG$level1)) #4

USDA_table_COG_level2 <- USDA_table_COG[,-c(1,3)]
USDA_table_COG_level2 <- aggregate(. ~ USDA_table_COG_level2$level2, 
                                      data = USDA_table_COG_level2[2:ncol(USDA_table_COG_level2)], 
                                      FUN = "sum")
rownames(USDA_table_COG_level2) <- USDA_table_COG_level2$`USDA_table_COG_level2$level2`
USDA_table_COG_level2 <- USDA_table_COG_level2[,-1]
#nrow(USDA_table_COG_level2) #25 as expected

#edgeR normalization of level2 level counts 
rnaseqMatrix_COG_level2 = round(USDA_table_COG_level2) 
rnaseqMatrix_COG_level2 = rnaseqMatrix_COG_level2[rowSums(cpm(rnaseqMatrix_COG_level2) > 1) >= 2,]
nrow(rnaseqMatrix_COG_level2) #25

#set up DGE 
exp_study_COG_level2 = DGEList(counts=rnaseqMatrix_COG_level2, group=group)
#calculate TMM normalization factors for library size adjustment
exp_study_COG_level2 = calcNormFactors(exp_study_COG_level2, method = "TMM")
level2_COG_TMM <- cpm.DGEList(exp_study_COG_level2)

COG_level2_sums <- data.frame(rowSums(level2_COG_TMM))
COG_level2_sums_ordered <- order(COG_level2_sums,
                                    COG_level2_sums$rowSums.level2_COG_TMM.)
COG_level2_sums <- COG_level2_sums[COG_level2_sums_ordered, ,drop = FALSE]

#setting 300,000 TMM counts as cutoff, aggregating all other level2s within 'Other'. Will give 9 taxa + aggregated group 
high_abun_level2s <- row.names(
     COG_level2_sums[COG_level2_sums$rowSums.level2_COG_TMM. > 400000, , drop = FALSE]
)
low_abun_level2s <- row.names(
     COG_level2_sums[COG_level2_sums$rowSums.level2_COG_TMM. < 400000, , drop = FALSE]
)

COG_level2_condensed <- level2_COG_TMM[row.names(level2_COG_TMM) %in% high_abun_level2s, , drop = FALSE]
COG_level2_condensed <- COG_level2_condensed[order(rowSums(COG_level2_condensed), decreasing = TRUE), ] 

Other <- level2_COG_TMM[row.names(level2_COG_TMM) %in% low_abun_level2s, , drop = FALSE]
Other <- colSums(Other)

COG_level2_condensed <- rbind(COG_level2_condensed, Other)
COG_level2_condensed_melt <- melt(COG_level2_condensed)
names(COG_level2_condensed_melt) <- c("COG level2", "Sample", "TMM")

#set level2s as factors, specify level2
COG_level2_condensed_melt$`COG level2` <- factor(COG_level2_condensed_melt$`COG level2`,
                                                       levels = unique(COG_level2_condensed_melt$`COG level2`))
levels(COG_level2_condensed_melt$`COG level2`)[1] <- "Translation, ribosomal \nstructure and biogenesis"
levels(COG_level2_condensed_melt$`COG level2`)[5] <- "Posttranslational modification, \nprotein turnover, chaperones"

COG_level2_condensed_melt$Sample <- factor(COG_level2_condensed_melt$Sample,
                                              levels = unique(COG_level2_condensed_melt$Sample))

col_palette2 <- rainbow(13)
col_palette2[5] <- "chartreuse3"
#col_palette[11] <- "black"

COG_level2_all_samples_TMM_abundance <- ggplot(data=COG_level2_condensed_melt, 
                                                  aes(x=Sample, y=TMM, fill=`COG level2`)) +
     geom_rect(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 4.5, fill = "gray", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 4.5, xmax = 9.5, fill = "orange", alpha = 0.5)+
     geom_rect(ymin = -Inf, ymax = Inf, xmin = 9.5, xmax = Inf, fill = "lightblue", alpha = 0.5)+
     geom_bar(stat="identity", width = 1, col = "black") +
     scale_fill_manual(values = col_palette2)+
     theme_bw()+
     annotate("text", x = 2.5, y = 1200000, label = "bold(CON)", parse = TRUE, size = 3.5)+
     annotate("text", x = 7, y = 1200000, label = "bold(WG)", parse = TRUE, size = 3.5)+
     annotate("text", x = 12.5, y = 1200000, label = "bold(RS)", parse = TRUE, size = 3.5)+
     theme(axis.text.y = element_text(colour = "black"),
           axis.text.x = element_blank(),
           legend.text = element_text(size=6),
           legend.title = element_text(size=10),
           legend.key.height = unit(0.5, "line"),
           legend.key.width = unit(0.75, "line"))

setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("USDA_COG_level2_all_samples_TMM_abundance.tiff", width = 4.5, height = 3.0, units = "in", compression = "lzw", res = 1200)
COG_level2_all_samples_TMM_abundance
dev.off()


#9. Plot log FC of highly expressed significant COGs ####
#subset COG table to just include the functions included for DE analysis
USDA_table_COG_toptags <- USDA_table_COG[USDA_table_COG$function. %in% rownames(fit_cog_CON_RS_toptags$table),]
nrow(USDA_table_COG_toptags) #2627 functions as expected
#order the functional annotations alphabetically in both objects 
USDA_table_COG_toptags <- USDA_table_COG_toptags[order(USDA_table_COG_toptags$function.),]
USDA_result_table <- fit_cog_CON_RS_toptags$table[order(rownames(fit_cog_CON_RS_toptags$table)),]
nrow(USDA_result_table) #2627 functions
sum(rownames(USDA_result_table) == USDA_table_COG_toptags$function.)
#add the full strings to result table, write to CSV
USDA_result_table_COG_string <- cbind(USDA_table_COG_toptags[,1:3], USDA_result_table)
View(USDA_result_table_COG_string)
write.csv(USDA_result_table_COG_string, file = "CON RS COG DE functions.csv", row.names = FALSE)
#sort by logFC
USDA_result_table_COG_string <- USDA_result_table_COG_string[order(USDA_result_table_COG_string$logFC, decreasing = TRUE),]
#extract significant features
USDA_result_table_COG_string_sig <- USDA_result_table_COG_string[USDA_result_table_COG_string$FDR < 0.05,]
#sort from highest to least normalized counts
View(head(USDA_result_table_COG_string_sig[order(USDA_result_table_COG_string_sig$logCPM, decreasing = TRUE),], n = 40))
#extract the 40 most abundant COGs
USDA_result_table_COG_string_sig_top40 <- head(USDA_result_table_COG_string_sig[order(USDA_result_table_COG_string_sig$logCPM, decreasing = TRUE),], n = 40)

#make a column with taxa names as factors so ggplot plots taxa in the specified order (sorted by logFC)
USDA_result_table_COG_string_sig_top40 <- USDA_result_table_COG_string_sig_top40[order(USDA_result_table_COG_string_sig_top40$logFC, decreasing = TRUE),]
USDA_result_table_COG_string_sig_top40$function. <- factor(USDA_result_table_COG_string_sig_top40$function., levels = unique(USDA_result_table_COG_string_sig_top40$function.))
View(USDA_result_table_COG_string_sig_top40)
unique(USDA_result_table_COG_string_sig_top40$level2) #16 unique functions at level 2 

ggplot(data = USDA_result_table_COG_string_sig_top40, aes(x = function., y =logFC)) +
     geom_bar(stat = "identity", aes(fill = level2))+
     #  theme(axis.text.y=element_text(color = result_table_fam_sig_only$color))+
     coord_flip()

#want to shorten some of the COG function names so it doesn't take up so much space.
#extract first 50 characters from COG function names 
USDA_result_table_COG_string_sig_top40$function.short <- substr(USDA_result_table_COG_string_sig_top40$function., 1, 50)
#wherever the nchar == 50 after subsetting, add "..." in the next column then concatenate the strings 
USDA_result_table_COG_string_sig_top40$dots <- ""
USDA_result_table_COG_string_sig_top40$dots[nchar(USDA_result_table_COG_string_sig_top40$function.short) == 50] <- "..."
#use paste to concatenate strings in adjacent columns
USDA_result_table_COG_string_sig_top40$function.short <- paste(USDA_result_table_COG_string_sig_top40$function.short,
                                                               USDA_result_table_COG_string_sig_top40$dots,
                                                                       sep = "")
#set as factor so plotted in the same order 
USDA_result_table_COG_string_sig_top40$function.short <- factor(USDA_result_table_COG_string_sig_top40$function.short,
                                                                        levels = USDA_result_table_COG_string_sig_top40$function.short)

#create a plot of the normalized counts for each gene (object already ordered properly so just plot the values)
counts_plot_by_COG <- ggplot(data = USDA_result_table_COG_string_sig_top40,
                             aes(x = logCPM, y = function.short)) + 
     geom_point(size = 2.5, aes(color = level2), show.legend = FALSE)+
     geom_segment(size = 1.25, aes(x = 10,
                                  xend = logCPM,
                                  y = function.short, 
                                  yend = function.short,
                                  color = level2), show.legend = FALSE)+
     theme_classic() +
     scale_y_discrete(limits=rev(levels(USDA_result_table_COG_string_sig_top40$function.short))) +
     xlab("log\n(TMM)")+ylab("")+xlim(c(10,15))+
     theme(axis.text.y = element_blank(),
           axis.text.x = element_text(color = "black"),
           plot.margin = unit(c(0.1,0.1,0.175,0.0), "cm")) # , , bottom, 

COG_barchart_function_logFC_top40 <- ggplot(data = USDA_result_table_COG_string_sig_top40, aes(x = function.short, y =logFC)) +
     geom_bar(stat = "identity", aes(fill = level2), color = "black")+
     theme(axis.text.y=element_text(color = "black", size = 8), 
           axis.text.x=element_text(color = "black"),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_line(color = "gray80"),
           panel.background = element_rect(fill = "white", colour = "black"),
           legend.key.size = unit(0.3, "cm"))+
     #  theme_bw()+
     xlab('COG Function')+ylab('log(Fold Change)\n(RS/CON)')+
     coord_flip()+
     scale_x_discrete(limits=rev(levels(USDA_result_table_COG_string_sig_top40$function.short)))+
     guides(fill = guide_legend(title = "COG level 2", ncol = 2))+
     theme(legend.position = c(0.10,-0.3),
           plot.margin = unit(c(0.1,3.5,3.5,0.1), "cm"))+
     annotation_custom(ggplotGrob(counts_plot_by_COG), xmin = -5.5, xmax = 41, ymin = 7.15, ymax = 14)

tiff("CON_RS.edgeR.COG.function_barchart_logFC.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 600)
COG_barchart_function_logFC_top40
dev.off()

pdf("CON_RS.edgeR.COG.function_barchart_logFC.pdf", width = 6.5, height = 6)
COG_barchart_function_logFC_top40
dev.off()


#heatmap of functional changes####

heatmap.2(prop.table(COG_TMM, 2),
          scale = "row",
          ColSideColors = c(rep("gray",4), rep("orange", 5), rep("blue", 6)),
          trace = "none")


#hierarchical clustering of RefSeq and COG data####
set.seed(1)
View(t(prop.table(COG_TMM, 2)))
#convert TMM counts to within-sample abundance, transpose (dist compares rows), calculate distance matrix, then UPGMA clustering, and finally make a dendrogram
COG_TMM_dendro <- t(prop.table(COG_TMM, 2)) %>%
     dist %>%
     hclust(method = "average") %>%
     as.dendrogram
plot(COG_TMM_dendro)     

COG_TMM_dendro_colored <- COG_TMM_dendro %>% 
     set("by_labels_branches_col", value = RS_sample_names, TF_values = c("blue", Inf)) %>% 
     set("by_labels_branches_col", value = CON_sample_names, TF_values = c("gray", Inf)) %>%
     set("by_labels_branches_col", value = WG_sample_names, TF_values = c("orange", Inf)) %>%
     set("by_labels_branches_lwd", value = RS_sample_names, TF_values = c(2, 1)) %>%
     set("by_labels_branches_lwd", value = c(RS_sample_names, WG_sample_names, CON_sample_names), TF_values = c(4, 1)) %>%
     ggplot
     
RefSeq_TMM_dendro <- t(prop.table(genus_refseq_TMM, 2)) %>%
     dist %>%
     hclust(method = "average") %>%
     as.dendrogram
plot(RefSeq_TMM_dendro)     

RefSeq_TMM_dendro_colored <- RefSeq_TMM_dendro %>% 
     set("by_labels_branches_col", value = RS_sample_names, TF_values = c("blue", Inf)) %>% 
     set("by_labels_branches_col", value = CON_sample_names, TF_values = c("gray", Inf)) %>%
     set("by_labels_branches_col", value = WG_sample_names, TF_values = c("orange", Inf)) %>%
     set("by_labels_branches_lwd", value = RS_sample_names, TF_values = c(2, 1)) %>%
     set("by_labels_branches_lwd", value = c(RS_sample_names, WG_sample_names, CON_sample_names), TF_values = c(4, 1)) %>%
     set("labels", rep("", 15)) %>%
     ggplot

#making tanglegram comparing RefSeq taxonomy and COG Function clustering
COG_RefSeq_dendlist <- dendlist(RefSeq_TMM_dendro, COG_TMM_dendro)
#making color vector for connectling lines 
RefSeq_TMM_dendro_labels <- RefSeq_TMM_dendro %>%
     set("labels_to_char")%>%
     labels
#manually adjusting order
RefSeq_TMM_dendro_labels <- data.frame(c("84","76", "18", "90", "40", "32", "48", "68", "58", "26", "24", "16", "60", "12", "20"))
names(RefSeq_TMM_dendro_labels) <- "sample"
RefSeq_TMM_dendro_labels$colors <- NA     
RefSeq_TMM_dendro_labels$colors[RefSeq_TMM_dendro_labels$sample %in% RS_sample_names] <- "blue"
RefSeq_TMM_dendro_labels$colors[RefSeq_TMM_dendro_labels$sample %in% WG_sample_names] <- "orange"
RefSeq_TMM_dendro_labels$colors[RefSeq_TMM_dendro_labels$sample %in% CON_sample_names] <- "gray"

COG_RefSeq_dendlist %>% 
     set("by_labels_branches_col", value = RS_sample_names, type = "any", TF_values = c("blue", Inf)) %>% 
     set("by_labels_branches_col", value = CON_sample_names, TF_values = c("gray", Inf)) %>%
     set("by_labels_branches_col", value = WG_sample_names, TF_values = c("orange", Inf)) %>%
     set("by_labels_branches_lwd", value = c(RS_sample_names, WG_sample_names, CON_sample_names), TF_values = c(2, 1)) %>%
#     set("labels_col", "white") %>%
     tanglegram(sort = TRUE, 
                common_subtrees_color_lines = FALSE, 
                highlight_distinct_edges  = FALSE, 
                highlight_branches_lwd = FALSE, 
                color_lines = rev(RefSeq_TMM_dendro_labels$colors),
                main_right = "COG Functions",
                main_left = "RefSeq Taxonomy",
                columns_width = c(5.5,2,5.5))
     legend(x = -1, y = 0.075, 
            pch = 15, 
            legend = c("CON", "WG", "RS"), 
            bty = "n", 
            horiz = FALSE, ncol = 1, cex = 1, pt.cex = 2,
            col = c("gray", "orange", "blue"),
            text.width = 0.1)

setwd("G:/My Drive/TRANSCEND/USDA/metatranscriptomics/new MG RAST analysis/all samples/")
tiff("tanglegram upgma refseq cog.tiff", width = 5, height = 3, units = "in", compression = "lzw", res = 1200)
COG_RefSeq_dendlist %>% 
     set("by_labels_branches_col", value = RS_sample_names, TF_values = c("blue", Inf)) %>% 
     set("by_labels_branches_col", value = CON_sample_names, TF_values = c("gray", Inf)) %>%
     set("by_labels_branches_col", value = WG_sample_names, TF_values = c("orange", Inf)) %>%
     set("by_labels_branches_lwd", value = c(RS_sample_names, WG_sample_names, CON_sample_names), TF_values = c(3, 1)) %>%
     set("labels_col", "white") %>%
     tanglegram(sort = TRUE, 
                common_subtrees_color_lines = FALSE, 
                highlight_distinct_edges  = FALSE, 
                highlight_branches_lwd = FALSE, 
                color_lines = rev(RefSeq_TMM_dendro_labels$colors),
                main_right = "COG\nFunctions",
                cex_main_right = 1.5,
                main_left = "RefSeq\nTaxonomy",
                cex_main_left = 1.5,
                margin_inner = 1.5,
                columns_width = c(5.25,2,5.25))%>%
     legend(x = -0.035, y = 0.85, 
            pch = 15, 
            legend = c("CON", "WG", "RS"), 
            bty = "n", 
            horiz = FALSE, ncol = 1, cex = 0.75, pt.cex = 2,
            col = c("gray", "orange", "blue"),
            text.width = 0.1)
dev.off()

pdf("tanglegram upgma refseq cog.pdf", width = 5, height = 3)
COG_RefSeq_dendlist %>% 
     set("by_labels_branches_col", value = RS_sample_names, TF_values = c("blue", Inf)) %>% 
     set("by_labels_branches_col", value = CON_sample_names, TF_values = c("gray", Inf)) %>%
     set("by_labels_branches_col", value = WG_sample_names, TF_values = c("orange", Inf)) %>%
     set("by_labels_branches_lwd", value = c(RS_sample_names, WG_sample_names, CON_sample_names), TF_values = c(3, 1)) %>%
     set("labels_col", "white") %>%
     tanglegram(sort = TRUE, 
                common_subtrees_color_lines = FALSE, 
                highlight_distinct_edges  = FALSE, 
                highlight_branches_lwd = FALSE, 
                color_lines = rev(RefSeq_TMM_dendro_labels$colors),
                main_right = "COG\nFunctions",
                cex_main_right = 1.5,
                main_left = "RefSeq\nTaxonomy",
                cex_main_left = 1.5,
                margin_inner = 1.5,
                columns_width = c(5.25,2,5.25))%>%
     legend(x = -0.035, y = 0.85, 
            pch = 15, 
            legend = c("CON", "WG", "RS"), 
            bty = "n", 
            horiz = FALSE, ncol = 1, cex = 0.75, pt.cex = 2,
            col = c("gray", "orange", "blue"),
            text.width = 0.1)
dev.off()


