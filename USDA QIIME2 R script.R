#USDA QIIME2 Analysis 
#December 2019 
#updated February 2020 with newer R version 
#Zachary Bendiks

#procesed R1 file in QIIME2, generated results at 17,500 seqs
#location of local files: /G:/My Drive/TRANSCEND/USDA/Qiime2/
#load.Rdata2("G:/My Drive/TRANSCEND/USDA/Qiime2/USDA QIIME2 R environment v2.RData")


#update R version using R GUI
#install.packages("installr")
#library(installr)
#updateR()

#need to 'car' package to perform ANOVA with type III SS for unbalanced designs
#can't do on lab PC (R version too old) so updating tablet and performing ANOVAs there
#install.packages("car")
library(car) #only runs on newer R versions 

#also need multcomp package for posthoc tests for unbalanced designs
#install.packages("multcomp")
#library(multcomp)
#NOTE 2/4/2020: statistician noted that posthoc tests aren't all that useful when interaction effects are reported 


#TAXONOMY SUMMARY AT SPECIES LEVEL####
#1.1 Import and clean up table####
asv_table <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/rarefied_table_17500_L7_USDA_filtered.csv", 
                        header = F,
                        stringsAsFactors = F,
                        comment.char = "") #includes # symbol in table instead of ignoring it 
#remove extra header line plus the empty "taxonomy" column
asv_table <- asv_table[-1,-ncol(asv_table)]
asv_table[1,1] <- "taxonomy" #rename "#OTU ID" column name
#set first row as column names them remove the first now 
names(asv_table) <- asv_table[1,]
asv_table <- asv_table[-1,]
#nrow(asv_table) #131 features 

#1.2 cleaning up taxa strings to give deepest classification####
asv_table$taxonomy <- gsub("__", "_", asv_table$taxonomy)

#split taxa strings by semicolon
OTU_split <- strsplit(asv_table$taxonomy, ";")
maxLen <- max(sapply(OTU_split, length))

#convert to matrix, fill with NAs to the sixth column (seventh col if species info included)
OTU_split_filled <- t(sapply(OTU_split, function(x)
  c(x, rep(NA, maxLen - length(x)))))

#convert to data frame
OTU_split_filled_df <- data.frame(OTU_split_filled)

#Step-wise addition of deepest taxonomic assignments from species to kingdom level
#create a new "taxonomy" column, add species-level assignment to it when length of the species-level character string > 4)
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X7)) > 4] <- 
  as.matrix(OTU_split_filled_df$X7[nchar(matrix(OTU_split_filled_df$X7)) > 4])
#add genus-level designations to taxonomy when the ASV has genus but not species assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X7)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X6)) > 4] <- 
  matrix(OTU_split_filled_df$X6[nchar(matrix(OTU_split_filled_df$X7)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X6)) > 4])
#add family-level designations to taxonomy when the ASV has family but not genus assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X6)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X5)) > 4] <- 
  matrix(OTU_split_filled_df$X5[nchar(matrix(OTU_split_filled_df$X6)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X5)) > 4])
#add order-level designations to taxonomy when the ASV has order but not family assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X5)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X4)) > 4] <- 
  matrix(OTU_split_filled_df$X4[nchar(matrix(OTU_split_filled_df$X5)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X4)) > 4])
#add class-level designations to taxonomy when the ASV has class but not order assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X4)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X3)) > 4] <- 
  matrix(OTU_split_filled_df$X3[nchar(matrix(OTU_split_filled_df$X4)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X3)) > 4])
#add phylum-level designations to taxonomy when the ASV has phylum but not class assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X3)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X2)) > 4] <- 
  matrix(OTU_split_filled_df$X2[nchar(matrix(OTU_split_filled_df$X3)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X2)) > 4])
#add kingdom-level designations to taxonomy when the ASV has kingdom but not phylum assignment   
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X2)) < 4 &
                               nchar(matrix(OTU_split_filled_df$X1)) > 4] <- 
  matrix(OTU_split_filled_df$X1[nchar(matrix(OTU_split_filled_df$X2)) < 4 &
                                  nchar(matrix(OTU_split_filled_df$X1)) > 4])

#need to add genus level to the species-level designations
OTU_split_filled_df$taxonomy[nchar(matrix(OTU_split_filled_df$X7)) > 4] <- #in taxonomy column where there is a species designation
  paste(OTU_split_filled_df$X6[nchar(matrix(OTU_split_filled_df$X7)) > 4], #paste the values in the genus column
        OTU_split_filled_df$X7[nchar(matrix(OTU_split_filled_df$X7)) > 4]) #paste the values in the species column 
#removing the prefixes from the new taxa names 
OTU_split_filled_df$taxonomy_clean <- gsub("[a-z]_", "", OTU_split_filled_df$taxonomy) #replaces prefixes with nothing
#adding cleaned up prefixes back to the asv table, replacing the old strings 
asv_table$taxonomy <- OTU_split_filled_df$taxonomy_clean

#1.3. identify the most and least abundant taxa ####
#str(asv_table) need to convert to numeric to combine everything
asv_table[,2:ncol(asv_table)] <- as.data.frame(lapply(asv_table[,2:ncol(asv_table)], as.numeric)) #convert sample sequence counts to numeric, output as a list, then convert to a dataframe 
#str(asv_table) counts as numeric, taxopnomy strings as characters 

#combine features with same taxa annotation
asv_table_SumByTaxa <- aggregate(asv_table[,2:ncol(asv_table)], 
                                 by = list(asv_table$taxonomy),
                                 FUN = "sum") 

taxa_sums <- data.frame(rowSums(asv_table_SumByTaxa[2:ncol(asv_table_SumByTaxa)])) #sum of all taxa 
taxa_sums$taxa <- asv_table_SumByTaxa$Group.1 #associating sum with actual taxa name
#order table based on sum 
taxa_sums <- taxa_sums[order(taxa_sums$rowSums.asv_table_SumByTaxa.2.ncol.asv_table_SumByTaxa..., 
                             decreasing = TRUE),]
taxa_sums_top15 <- taxa_sums[1:15,]
sum(as.numeric(taxa_sums_top15$rowSums.asv_table_SumByTaxa.2.ncol.asv_table_SumByTaxa...))
taxa_sums_low_abun <- taxa_sums[16:nrow(taxa_sums),]
#sum low abundance taxa into "Other", add to the high taxa object 
sum_low_abundance <- sum(taxa_sums_low_abun$rowSums.asv_table_SumByTaxa.2.ncol.asv_table_SumByTaxa...)
taxa_sums_top15[nrow(taxa_sums_top15) + 1,] <- c(sum_low_abundance, "Other")

#Convert table to relative abundance (divide by the sum of counts in each column)
#colSums(asv_table_SumByTaxa[, 2:ncol(asv_table_SumByTaxa)]) #still 17,500 for each sample
asv_table_SumByTaxa_prop <- asv_table_SumByTaxa[, 2:ncol(asv_table_SumByTaxa)] / colSums(asv_table_SumByTaxa[, 2:ncol(asv_table_SumByTaxa)])
asv_table_SumByTaxa_prop$taxonomy <- asv_table_SumByTaxa$Group.1

#Subset table to just include taxa of interest 
asv_table_SumByTaxa_prop_subset <- asv_table_SumByTaxa_prop[asv_table_SumByTaxa_prop$taxonomy %in%
                                                              taxa_sums_top15$taxa,]
#reorder subset table to list taxa in descending order
asv_table_SumByTaxa_prop_subset <- asv_table_SumByTaxa_prop_subset[match(taxa_sums[1:15,]$taxa,
                                                                         asv_table_SumByTaxa_prop_subset$taxonomy),]

column_order <- order(asv_table_SumByTaxa_prop_subset[1,-ncol(asv_table_SumByTaxa_prop_subset)]) #order of columns based on Bacteroides abundance 
asv_table_SumByTaxa_prop_subset_ordered <- asv_table_SumByTaxa_prop_subset[,column_order]
asv_table_SumByTaxa_prop_subset_ordered$taxonomy <- asv_table_SumByTaxa_prop_subset$taxonomy
#make row names the taxonomy, remove column from table 
rownames(asv_table_SumByTaxa_prop_subset_ordered) <- asv_table_SumByTaxa_prop_subset_ordered$taxonomy
asv_table_SumByTaxa_prop_subset_ordered_taxa_rownames <- asv_table_SumByTaxa_prop_subset_ordered[,-ncol(asv_table_SumByTaxa_prop_subset_ordered)]

#1.4 Associate group information with sample IDs using mapping file ####
asv_mapfile <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/USDA_mapfile.csv", 
                          header = T, #note that I'm including the header this time 
                          stringsAsFactors = F,
                          comment.char = "") #includes # symbol in table instead of ignoring it 

#subset mapping file to only include samples in the table 
asv_mapfile_subset <- asv_mapfile[asv_mapfile$X.SampleID %in% names(asv_table_SumByTaxa_prop_subset_ordered_taxa_rownames),]
#reorder rows in mapfile to match order of samples in the table columns 
asv_mapfile_subset <- asv_mapfile_subset[match(names(asv_table_SumByTaxa_prop_subset_ordered_taxa_rownames),
                                               asv_mapfile_subset$X.SampleID),]
#transposing table (columns list all data of importance, with individual samples as rows )
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames <- data.frame(t(asv_table_SumByTaxa_prop_subset_ordered_taxa_rownames))
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames$Diet <- asv_mapfile_subset$Code2
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames$SampleID <- rownames(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames)

#1.5 Reformat table with reshape2 ####
#install.packages("reshape2")
library(reshape2)
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted <- melt(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames,
                                                                      id.vars = c("SampleID","Diet"))
#remove the "." itnroduced in the taxa names 
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable <- gsub("\\.", " ", tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable)
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable[tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "S24 7"] <- "S24-7"
#ggplot automatically reorders legend, so I need to adjust at factor level  
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable <- as.factor(
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable
)

#levels(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable)

tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable <- factor(
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable, 
  levels = taxa_sums_top15[-nrow(taxa_sums_top15),]$taxa
)

#order of samples is also wrong, so need to adjust at factor level  
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$SampleID <- factor(
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$SampleID, 
  levels = names(asv_table_SumByTaxa_prop_subset_ordered_taxa_rownames)
)

#1.6 generating plot with ggplot2 ####
#install.packages("ggplot2")
library(ggplot2)

#when ggplot makes figure legend, it automatically adds gaps between the entries.  I want to remove that behavior
#from: https://github.com/tidyverse/ggplot2/issues/2844
draw_key_polygon2 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
#register new key drawing function, effect is global and persistent throughout R session!
#automatically applied when geom_bar called 
GeomBar$draw_key = draw_key_polygon2 

#default rainbow color scheme makes it hard to resolve the different green colors. I manually adjust to increase the contrast 
col_palette <- rainbow(15)
col_palette[6] <- "green3"
col_palette[8] <- "cyan"
col_palette[14] <- "hotpink"

#rename diets and specify diet order as factor levels for facet grid 
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet[tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet == "AC"] <- "CON"
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet[tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet == "WWGC"] <- "WG"
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet[tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet == "HM260"] <- "RS"
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet[tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet == "HMWG"] <- "WG+RS"

tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet_f <- factor(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Diet,
       levels = c("CON", "WG", "RS", "WG+RS"))

#making plot 
usda_taxa_summary_by_group <- ggplot(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted,
                                       aes(x = SampleID, y = value)) +
  geom_bar(stat="identity",aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = col_palette)+
  facet_grid(~ Diet_f, space = "free", scales = "free")+
  theme_bw()+
  ylab("Relative Proportion")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5,"cm"))+
  ylim(c(0,1))+
  guides(fill=guide_legend(title="Taxonomy"))+
  ggtitle("Rarefied to 17,500 seqs")

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by group.pdf", width = 10, height = 5)
usda_taxa_summary_by_group
dev.off()  

tiff("taxa summary usda facet by group.tiff", width = 10, height = 5, units = "in", res = 1200)
usda_taxa_summary_by_group
dev.off()  

#want to incorporate fat information 
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Fat <- ""
HFD_table  <- tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted[grep("H[0-9]{1,2}$", tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$SampleID),] 
HFD_table$Fat <- "HFD"
MFD_table  <- tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted[grep("M[0-9]{1,2}$", tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$SampleID),] 
MFD_table$Fat <- "MFD"
melted_table_w_fat <- rbind(HFD_table, MFD_table)
nrow(melted_table_w_fat) == nrow(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted)

usda_taxa_summary_by_group_by_fat <- ggplot(melted_table_w_fat,
       aes(x = SampleID, y = value)) +
  geom_rect(data = melted_table_w_fat,
            fill = "white", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)+
    geom_bar(stat="identity",aes(fill = variable), width = 0.9)+
  scale_fill_manual(values = col_palette)+
  facet_wrap(Diet_f ~ Fat, scale = "free_x", ncol = 2)+
  theme_bw()+
  ylab("Relative Proportion")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5,"cm"))+
  ylim(c(0,1))+
  guides(fill=guide_legend(title="Taxonomy"))

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by group by fat.pdf", width = 6.0, height = 7.5)
usda_taxa_summary_by_group_by_fat
dev.off()  

tiff("taxa summary usda facet by group by fat.tiff", width = 6.0, height = 7.5, units = "in", res = 1200)
usda_taxa_summary_by_group_by_fat
dev.off()  


#1.7 remaking taxa plot but faceting plot by individual taxa ####
#addind a column I can use to color the facet backgrounds by enrichment 
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Enrichment <- ""

#identify RS-enriched taxa. Focusing on taxa levels that were explicitly enriched, not aggregated enrichment like how LEFSE does it 
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Enrichment[
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "S24-7" |
    tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Bacteroides uniformis"] <- "RS"
#identify RS-reduced taxa  
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Enrichment[
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Clostridiales" |
    tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Ruminococcaceae" |
    tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Lachnospiraceae" |
    tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Coprococcus" |
    tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Oscillospira"] <- "CON"
#identify WG-enriched taxa  
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$Enrichment[
  tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable == "Lactobacillus"] <- "WG"

#put newline expression in place of spaces so facet headers aren't so long for species, then set as factor so taxa order is preserved
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable2 <- gsub(" ", "\n", 
                                                                                tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable)
  
tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable2 <- factor(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable2,
     ordered = TRUE,
     levels = unique(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted$variable2))

#generate plot using geom_rect to color backgrounds of individual facets 
usda_taxa_summary_faceted_by_taxa <- ggplot(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted,
       aes(x = Diet_f, y = value)) +
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "RS"),
                          fill = "lightskyblue1", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "CON"),
            fill = "lightgoldenrod", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "WG"),
            fill = "lightpink", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_boxplot(col = NA, alpha = 0, outlier.color = "black", outlier.alpha = 1, outlier.size = 1) +
  geom_bar(stat = "summary", fun.y="mean", aes(fill = Diet_f))+ #plot the means of each taxa within each group 
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.5, col = rep(c("azure4", "red", "blue", "purple"), 15))+ #add standard error bars and repeat the color scheme 15x to fill in the plot
  facet_wrap(~ variable2, ncol = 5)+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x.bottom = element_blank(),
        panel.grid.major = element_line(color = "gray65"),
        panel.grid.minor = element_line(color = "white"))+
  labs(fill = "Diet")+
  ylab("Relative Proportion")+xlab(NULL)
  
setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by taxa.pdf", width = 7.5, height = 5)
usda_taxa_summary_faceted_by_taxa
dev.off()  

tiff("taxa summary usda facet by taxa.tiff", width = 7.5, height = 5, units = "in", res = 1200)
usda_taxa_summary_faceted_by_taxa
dev.off()  

#free scales to see all taxa more clearly 
usda_taxa_summary_faceted_by_taxa_free_scale <- ggplot(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted,
                                            aes(x = Diet_f, y = value)) + 
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "RS"),
            fill = "lightskyblue1", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "CON"),
            fill = "lightgoldenrod", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "WG"),
            fill = "lightpink", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_boxplot(col = NA, alpha = 0, outlier.color = "black", outlier.alpha = 1, outlier.size = 1) +
  geom_bar(stat = "summary", fun.y="mean", aes(fill = Diet_f))+ #plot the means of each taxa within each group 
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.5, col = rep(c("azure4", "red", "blue", "purple"), 15))+ #add standard error bars and repeat the color scheme 15x to fill in the plot
  facet_wrap(~ variable2, ncol = 5, scales = "free_y")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x.bottom = element_blank(),
        panel.grid.major = element_line(color = "gray65"),
        panel.grid.minor = element_line(color = "white"))+
  labs(fill = "Diet")+
  ylab("Relative Proportion")+xlab(NULL)

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by taxa free scale.pdf", width = 8.5, height = 5)
usda_taxa_summary_faceted_by_taxa_free_scale
dev.off()  

tiff("taxa summary usda facet by taxa free scale.tiff", width = 8.5, height = 5, units = "in", res = 1200)
usda_taxa_summary_faceted_by_taxa_free_scale
dev.off()  

usda_taxa_summary_faceted_by_taxa_free_scale_no_outliers <- ggplot(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted,
                                                       aes(x = Diet_f, y = value)) + 
#  geom_boxplot(col = NA, alpha = 0, outlier.color = "black", outlier.alpha = 1, outlier.size = 1) +
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "RS"),
            fill = "lightskyblue1", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "CON"),
            fill = "lightgoldenrod", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "WG"),
            fill = "lightpink", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.05)+
  geom_bar(stat = "summary", fun.y="mean", aes(fill = Diet_f))+ #plot the means of each taxa within each group 
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.5, col = rep(c("azure4", "red", "blue", "purple"), 15))+ #add standard error bars and repeat the color scheme 15x to fill in the plot
  facet_wrap(~ variable2, ncol = 5, scales = "free_y")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x.bottom = element_blank(),
        panel.grid.major = element_line(color = "gray65"),
        panel.grid.minor = element_line(color = "white"))+
  labs(fill = "Diet")+
  ylab("Relative Proportion")+xlab(NULL)  

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by taxa free scale no outliers.pdf", width = 8.5, height = 5)
usda_taxa_summary_faceted_by_taxa_free_scale_no_outliers
dev.off()  

tiff("taxa summary usda facet by taxa free scale no outliers.tiff", width = 8.5, height = 5, units = "in", res = 1200)
usda_taxa_summary_faceted_by_taxa_free_scale_no_outliers
dev.off()  


usda_taxa_summary_faceted_by_taxa_free_scale_boxplots <- ggplot(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted,
                                                       aes(x = Diet_f, y = value)) + 
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "RS"),
            fill = "lightskyblue1", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "CON"),
            fill = "lightpink", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == "WG"),
            fill = "lightgoldenrod2", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)+
  geom_rect(data = subset(tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames_melted, Enrichment == ""),
            fill = "white", xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)+
  geom_boxplot(aes(fill = Diet_f),
               col = "black",
               outlier.alpha = 1, 
               outlier.size = 0.8,
               width = 0.9,
               lwd = 0.20) +
  facet_wrap(~ variable2, ncol = 5, scales = "free_y")+
  scale_fill_manual(values = c("gray75", "darkorange", "blue", "green"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x.bottom = element_blank(),
        axis.text.y = element_text(size = 7.0),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 7.5),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white", linetype = "solid"),
        panel.grid = element_line(color = NULL))+
  labs(fill = "Diet")+
  ylab("Relative Proportion")+xlab(NULL)

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("taxa summary usda facet by taxa free scale boxplots.pdf", width = 6.5, height = 4)
usda_taxa_summary_faceted_by_taxa_free_scale_boxplots
dev.off()  

tiff("taxa summary usda facet by taxa free scale boxplots.tiff", width = 6.5, height = 4, units = "in", res = 1200)
usda_taxa_summary_faceted_by_taxa_free_scale_boxplots
dev.off()  


#ORDINATION PLOT OF BRAY CURTIS DM WITH CENTROIDS AND CONF. ELLIPSES ####
#install.packages("vegan")
#install.packages("MASS")
#install.packages("RColorBrewer")

library(vegan)
library(MASS)
library(RColorBrewer)
#2.1 Formatting input objects ####
#converting QIIME dm into 'dist' objects that can be used in vegan 
#mapping file header = FALSE so first data line isn't read as header 
asv_mapfile
QIIMEdm <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/bray_curtis_DM_USDA_filtered.csv",
                    stringsAsFactors = FALSE, 
                    header = TRUE,
                    row.names = 1) #this flag sets the first columnn into the row names

#convert distance matrix to 'dist' object type 
QIIMEdist <- as.dist(QIIMEdm)
MDS <- isoMDS(QIIMEdist)
ordiplot(MDS)
#generates orincipal coordinates 
PCoA <- cmdscale(QIIMEdist)
#two wyas to plot
ordiplot(PCoA)
plot(PCoA)
#create vector of sample names in the order they appear in the DM
sample_names <- rownames(QIIMEdm)
#subset mapping file by sample names in DM
nrow(asv_mapfile) #96
nrow(QIIMEdm) #93 
asv_mapfile_sub <- asv_mapfile[asv_mapfile$X.SampleID %in% sample_names,]
nrow(asv_mapfile_sub) #93 
#reorder map object based on sample names in PCoA using the match argument 
asv_mapfile_sub_ordered <- asv_mapfile_sub[match(sample_names, asv_mapfile_sub$X.SampleID),]
asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "AC"] <- "CON"
asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WWGC"] <- "WG"
asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "HM260"] <- "RS"
asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "HMWG"] <- "WG+RS"
asv_mapfile_sub_ordered$Code2 <- factor(asv_mapfile_sub_ordered$Code2,
                                        levels = c("CON","WG", "RS", "WG+RS"))
#adding % variance to column names so they appear in plot 
colnames(PCoA) <- c("PC1 - 27.4%", "PC2 - 9.5%")

#2.2 generating plot ####
tiff("vegan PCoA usda qiime2.tiff", width = 4.25, height = 4.25, units = "in", compression = "lzw", res = 1200)

par(mgp=c(1,0,0) ,mar=c(2,2,0.5,0.5)) #clear all plots in viewer to reset parameters 

ordiplot(PCoA, axes = FALSE, frame.plot = TRUE)
#want SD ellipses plotted first so pts and centroid lines are plotted over them 
ordiellipse(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "azure4", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "red", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "purple", lwd = 1.75, kind = "sd", draw = "polygon")
#adding centroid lines for each group
ordispider(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "azure4", lwd = 1.25)
ordispider(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "red", lwd = 1.25)
ordispider(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.25)
ordispider(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "purple", lwd = 1.25)
#coloring individual points 
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "azure4", cex = 0.8))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "red", cex = 0.8))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "blue", cex = 0.8))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "purple", cex = 0.8))
#add a generic legend 
legend(x = -0.45, y = 0.54, 
       pch = 21, 
       legend = c("CON", "WG", "RS", "WG+RS"), 
       bty = "n", 
       horiz = FALSE, ncol = 2, cex = 0.8, pt.cex = 1.0,
       pt.bg = c("azure4", "red", "blue", "purple"), 
       col = "black",
       text.width = 0.1)

dev.off()

#2.3 Use vegan package to calculate PERMDISP (dispersal of points about centroid, use to argue that RS samples more variable) ####
nrow(asv_mapfile_sub_ordered)
nrow(QIIMEdm)
asv_mapfile_sub_ordered_matched <- asv_mapfile_sub_ordered[match(asv_mapfile_sub_ordered$X.SampleID,
                                      rownames(QIIMEdm)),]

betadisper(QIIMEdist,
           asv_mapfile_sub_ordered_matched$Code,
           type = "centroid")

PERMDISP <- betadisper(QIIMEdist,
           asv_mapfile_sub_ordered_matched$Code2,
           type = "centroid")

anova(PERMDISP)
permutest(PERMDISP, pairwise = TRUE, permutations = 999)

#three way ANOVA
PERMDISP_distances <- melt(PERMDISP$distances)
PERMDISP_distances$Starch <- ""
PERMDISP_distances$Starch[grep("AC",rownames(PERMDISP_distances))] <- "No RS"
PERMDISP_distances$Starch[grep("WWGC",rownames(PERMDISP_distances))] <- "No RS"
PERMDISP_distances$Starch[grep("HM260",rownames(PERMDISP_distances))] <- "RS"
PERMDISP_distances$Starch[grep("HMWG",rownames(PERMDISP_distances))] <- "RS"
PERMDISP_distances$WholeGrain <- ""
PERMDISP_distances$WholeGrain[grep("AC",rownames(PERMDISP_distances))] <- "No WG"
PERMDISP_distances$WholeGrain[grep("WWGC",rownames(PERMDISP_distances))] <- "WG"
PERMDISP_distances$WholeGrain[grep("HM260",rownames(PERMDISP_distances))] <- "No WG"
PERMDISP_distances$WholeGrain[grep("HMWG",rownames(PERMDISP_distances))] <- "WG"
PERMDISP_distances$Fat <- ""
#this regex says look for H or M, followed by 1 or 2 numbers, then the end of the string
PERMDISP_distances[grep("H[0-9]{1,2}$",rownames(PERMDISP_distances)), ] #only HFD samples
PERMDISP_distances$Fat[grep("H[0-9]{1,2}$",rownames(PERMDISP_distances))] <- "HighFat"
PERMDISP_distances[grep("M[0-9]{1,2}$",rownames(PERMDISP_distances)),] #only MFD samples
PERMDISP_distances$Fat[grep("M[0-9]{1,2}$",rownames(PERMDISP_distances))] <- "MedFat"

PERMDISP_2way_AOV_model <- aov(value ~ Starch + WholeGrain + Fat + Starch:WholeGrain + Starch:Fat + WholeGrain:Fat + Starch:WholeGrain:Fat, data = PERMDISP_distances)
summary(PERMDISP_2way_AOV_model)
plot(PERMDISP_2way_AOV_model)
TukeyHSD(PERMDISP_2way_AOV_model, which = "Starch:WholeGrain")

#three way ANOVA type III SS for unbalanced designs using car package 
#cannot model three way interaction (aliased coefficients)
alias(lm(value ~ Starch + 
           WholeGrain + 
           Fat + 
           Starch:WholeGrain + 
           Starch:Fat + 
           WholeGrain:Fat + 
           Starch:WholeGrain:Fat,
         data = PERMDISP_distances))

#create model with contrasts for type III SS 
#really not sure what I'm doing, seems like an awful lot of work for something that isn't going to improve the analysis much 
options(contrasts = c("contr.sum", "contr.poly"))
PERMDISP_distances$Starch <- factor(PERMDISP_distances$Starch)
PERMDISP_distances$WholeGrain <- factor(PERMDISP_distances$WholeGrain)
PERMDISP_distances$Fat <- factor(PERMDISP_distances$Fat)

PERMDISP_2way_AOV_model_type3 <- lm(value ~ Starch * 
                                       WholeGrain * 
                                       Fat -
                                       Starch:WholeGrain:Fat,
                                     data = PERMDISP_distances)

lm(value ~ Starch + 
     WholeGrain + 
     Fat +
     Starch:WholeGrain + 
     WholeGrain:Fat +
     Starch:Fat,
   data = PERMDISP_distances)

#see: https://www.r-bloggers.com/anova-%E2%80%93-type-iiiiii-ss-explained/

#2.4 Same ordination plot but colored by pH (IN _PROGRESS) ####
color_continuous <- colorRampPalette(c("green", "red"))(nrow(asv_mapfile_sub_ordered))
asv_mapfile_sub_ordered_pH <- asv_mapfile_sub_ordered[order(asv_mapfile_sub_ordered$Cecal_pH),]

ordiplot(PCoA)
with(asv_mapfile_sub_ordered_pH,
     points(PCoA[,1:2, drop = FALSE], 
            col = color_continuous, 
            bg = color_continuous,
            pch = 21, cex = 1.4))

#2.5 adding fat as shape factor to PCoA plot ####
tiff("vegan PCoA usda qiime2 w fat shape 2.tiff", width = 3.0, height = 3.0, units = "in", compression = "lzw", res = 1200)

par(mgp=c(1,0,0) ,mar=c(2,2,0.5,0.5)) #clear all plots in viewer to reset parameters 

ordiplot(PCoA, axes = FALSE, frame.plot = TRUE, xlab = NA, ylab = NA, type = "none")
#want SD ellipses plotted first so pts and centroid lines are plotted over them 
ordiellipse(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "gray75", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "darkorange", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "green", lwd = 1.75, kind = "sd", draw = "polygon")
#adding centroid lines for each group
ordispider(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "azure4", lwd = 1.25)
ordispider(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "darkorange", lwd = 1.25)
ordispider(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.25)
ordispider(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "green", lwd = 1.25)
#coloring individual points but changing shape based on fat 
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEACM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "gray75", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEWWGCM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "darkorange", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHM260M", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "blue", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHMWGM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "green", cex = 1))

with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEACH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "gray75", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEWWGCH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "darkorange", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHM260H", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "blue", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHMWGH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "green", cex = 1))
#add a generic legend for diet
legend(x = -0.375, y = 0.565, 
       pch = 15, 
       legend = c("CON", "WG", "RS", "WG+RS"), 
       bty = "n", 
       horiz = FALSE, ncol = 2, cex = 0.75, pt.cex = 0.75,
       col = c("gray75", "darkorange", "blue", "green"), 
       text.width = 0.1)
#add a generic legend for fat
legend(x = -0.375, y = 0.425, 
       pch = c(17, 19), 
       legend = c("MFD", "HFD"), 
       bty = "n", 
       horiz = FALSE, ncol = 2, cex = 0.75, pt.cex = 0.75,
       col = "black",
       text.width = 0.1)
#adding x and y axis labels (do it this way so you cna adjust the size) 
mtext(side = 1, "PC1 - 27.4%", cex = 1)
mtext(side = 2, "PC2 - 9.5%", cex = 1)

dev.off()


#make a PDF version too 
pdf("vegan PCoA usda qiime2 w fat shape 2.pdf", width = 3.0, height = 3.0)

par(mgp=c(1,0,0) ,mar=c(2,2,0.5,0.5)) #clear all plots in viewer to reset parameters 

ordiplot(PCoA, axes = FALSE, frame.plot = TRUE, xlab = NA, ylab = NA, type = "none")
#want SD ellipses plotted first so pts and centroid lines are plotted over them 
ordiellipse(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "gray75", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "darkorange", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.75, kind = "sd", draw = "polygon")
ordiellipse(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "green", lwd = 1.75, kind = "sd", draw = "polygon")
#adding centroid lines for each group
ordispider(PCoA[grep("KEEAC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "CON"], col = "azure4", lwd = 1.25)
ordispider(PCoA[grep("KEEWWGC", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG"], col = "darkorange", lwd = 1.25)
ordispider(PCoA[grep("KEEHM260", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "RS"], col = "blue", lwd = 1.25)
ordispider(PCoA[grep("KEEHMWG", rownames(PCoA)), 1:2, drop = FALSE], group = asv_mapfile_sub_ordered$Code2[asv_mapfile_sub_ordered$Code2 == "WG+RS"], col = "green", lwd = 1.25)
#coloring individual points but changing shape based on fat 
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEACM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "gray75", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEWWGCM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "darkorange", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHM260M", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "blue", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHMWGM", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 24, bg = "green", cex = 1))

with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEACH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "gray75", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEWWGCH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "darkorange", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHM260H", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "blue", cex = 1))
with(asv_mapfile_sub_ordered, points(PCoA[grep("KEEHMWGH", rownames(PCoA)), 1:2, drop = FALSE], col = "black", pch = 21, bg = "green", cex = 1))
#add a generic legend for diet
legend(x = -0.375, y = 0.565, 
       pch = 15, 
       legend = c("CON", "WG", "RS", "WG+RS"), 
       bty = "n", 
       horiz = FALSE, ncol = 2, cex = 0.75, pt.cex = 0.75,
       col = c("gray75", "darkorange", "blue", "green"), 
       text.width = 0.1)
#add a generic legend for fat
legend(x = -0.375, y = 0.425, 
       pch = c(17, 19), 
       legend = c("MFD", "HFD"), 
       bty = "n", 
       horiz = FALSE, ncol = 2, cex = 0.75, pt.cex = 0.75,
       col = "black",
       text.width = 0.1)
#adding x and y axis labels (do it this way so you cna adjust the size) 
mtext(side = 1, "PC1 - 27.4%", cex = 1)
mtext(side = 2, "PC2 - 9.5%", cex = 1)

dev.off()



#2.6 three way ANOVA of BC distances to CON group ####
QIIMEdm_melted <- cbind(names(QIIMEdm), melt(QIIMEdm))
#remove distances between teh same sample
QIIMEdm_melted <- QIIMEdm_melted[QIIMEdm_melted$`names(QIIMEdm)` != 
                 QIIMEdm_melted$variable,]

#obtain just the intra-CON-M distances 
#identify indices where specific groups are found, then use those indices to separate the different comparisons 
CON_M_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                           grep("KEEACM", QIIMEdm_melted$variable))
QIIMEdm_melted_CON_M <- QIIMEdm_melted[CON_M_indices,]
QIIMEdm_melted_CON_M$Starch <- "No RS"
QIIMEdm_melted_CON_M$WG <- "No WG"
QIIMEdm_melted_CON_M$Fat <- "MFD"
QIIMEdm_melted_CON_M$Group <- "CON-MFD"
#every value is duplicated, so I organize by value and delete every even-numbered row 
QIIMEdm_melted_CON_M <- QIIMEdm_melted_CON_M[order(QIIMEdm_melted_CON_M$value, decreasing = TRUE),]
toDelete <- seq(from = 1, to = nrow(QIIMEdm_melted), by = 2)
QIIMEdm_melted_CON_M <- QIIMEdm_melted_CON_M[-toDelete, ]
nrow(QIIMEdm_melted_CON_M)
#obtain the CON-M:CON-H distances 
CON_H_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                           grep("KEEACH", QIIMEdm_melted$variable))
#subsetting this way prevents duplicates from appearing 
QIIMEdm_melted_CON_H <- QIIMEdm_melted[CON_H_indices,]
QIIMEdm_melted_CON_H <- QIIMEdm_melted_CON_H[order(QIIMEdm_melted_CON_H$value, decreasing = TRUE),]
QIIMEdm_melted_CON_H$Starch <- "No RS"
QIIMEdm_melted_CON_H$WG <- "No WG"
QIIMEdm_melted_CON_H$Fat <- "HFD"
QIIMEdm_melted_CON_H$Group <- "CON-HFD"
nrow(QIIMEdm_melted_CON_H)

#obtain the CON-M:WG-M distances 
WG_M_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                           grep("KEEWWGCM", QIIMEdm_melted$variable))
QIIMEdm_melted_WG_M <- QIIMEdm_melted[WG_M_indices,]
QIIMEdm_melted_WG_M <- QIIMEdm_melted_WG_M[order(QIIMEdm_melted_WG_M$value, decreasing = TRUE),]
QIIMEdm_melted_WG_M$Starch <- "No RS"
QIIMEdm_melted_WG_M$WG <- "WG"
QIIMEdm_melted_WG_M$Fat <- "MFD"
QIIMEdm_melted_WG_M$Group <- "WG-MFD"
nrow(QIIMEdm_melted_WG_M)

#obtain the CON-M:WG-H distances 
WG_H_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                          grep("KEEWWGCH", QIIMEdm_melted$variable))
QIIMEdm_melted_WG_H <- QIIMEdm_melted[WG_H_indices,]
QIIMEdm_melted_WG_H <- QIIMEdm_melted_WG_H[order(QIIMEdm_melted_WG_H$value, decreasing = TRUE),]
QIIMEdm_melted_WG_H$Starch <- "No RS"
QIIMEdm_melted_WG_H$WG <- "WG"
QIIMEdm_melted_WG_H$Fat <- "HFD"
QIIMEdm_melted_WG_H$Group <- "WG-HFD"
nrow(QIIMEdm_melted_WG_H)

#obtain the CON-M:RS-M distances 
RS_M_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                          grep("KEEHM260M", QIIMEdm_melted$variable))
QIIMEdm_melted_RS_M <- QIIMEdm_melted[RS_M_indices,]
QIIMEdm_melted_RS_M <- QIIMEdm_melted_RS_M[order(QIIMEdm_melted_RS_M$value, decreasing = TRUE),]
QIIMEdm_melted_RS_M$Starch <- "RS"
QIIMEdm_melted_RS_M$WG <- "No WG"
QIIMEdm_melted_RS_M$Fat <- "MFD"
QIIMEdm_melted_RS_M$Group <- "RS-MFD"
nrow(QIIMEdm_melted_RS_M)

#obtain the CON-M:RS-H distances 
RS_H_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                          grep("KEEHM260H", QIIMEdm_melted$variable))
QIIMEdm_melted_RS_H <- QIIMEdm_melted[RS_H_indices,]
QIIMEdm_melted_RS_H <- QIIMEdm_melted_RS_H[order(QIIMEdm_melted_RS_H$value, decreasing = TRUE),]
QIIMEdm_melted_RS_H$Starch <- "RS"
QIIMEdm_melted_RS_H$WG <- "No WG"
QIIMEdm_melted_RS_H$Fat <- "HFD"
QIIMEdm_melted_RS_H$Group <- "RS-HFD"
nrow(QIIMEdm_melted_RS_H)

#obtain the CON-M:WGRS-M distances 
WGRS_M_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                          grep("KEEHMWGM", QIIMEdm_melted$variable))
QIIMEdm_melted_WGRS_M <- QIIMEdm_melted[WGRS_M_indices,]
QIIMEdm_melted_WGRS_M <- QIIMEdm_melted_WGRS_M[order(QIIMEdm_melted_WGRS_M$value, decreasing = TRUE),]
QIIMEdm_melted_WGRS_M$Starch <- "RS"
QIIMEdm_melted_WGRS_M$WG <- "WG"
QIIMEdm_melted_WGRS_M$Fat <- "MFD"
QIIMEdm_melted_WGRS_M$Group <- "WG+RS-MFD"
nrow(QIIMEdm_melted_WGRS_M)

#obtain the CON-M:WGRS-H distances 
WGRS_H_indices <- intersect(grep("KEEACM", QIIMEdm_melted$`names(QIIMEdm)`),
                            grep("KEEHMWGH", QIIMEdm_melted$variable))
QIIMEdm_melted_WGRS_H <- QIIMEdm_melted[WGRS_H_indices,]
QIIMEdm_melted_WGRS_H <- QIIMEdm_melted_WGRS_H[order(QIIMEdm_melted_WGRS_H$value, decreasing = TRUE),]
QIIMEdm_melted_WGRS_H$Starch <- "RS"
QIIMEdm_melted_WGRS_H$WG <- "WG"
QIIMEdm_melted_WGRS_H$Fat <- "HFD"
QIIMEdm_melted_WGRS_H$Group <- "WG+RS-HFD"
nrow(QIIMEdm_melted_WGRS_H)

#combine all group dataframes together 
QIIMEdm_melted_3factors <- rbind(QIIMEdm_melted_CON_M,
      QIIMEdm_melted_CON_H,
      QIIMEdm_melted_WG_M,
      QIIMEdm_melted_WG_H,
      QIIMEdm_melted_RS_M,
      QIIMEdm_melted_RS_H,
      QIIMEdm_melted_WGRS_M,
      QIIMEdm_melted_WGRS_H)

#three way ANOVA with type III SS
options(contrasts = c("contr.sum", "contr.poly")) #must change these options to get legitimate type III SS results 
#first create model 
QIIMEdm_3wayANOVA_model <- aov(value ~ Starch + 
                           WG + 
                           Fat + 
                           Starch:WG + 
                           Starch:Fat + 
                           WG:Fat + 
                           Starch:WG:Fat,
                         data = QIIMEdm_melted_3factors)
summary(QIIMEdm_3wayANOVA_model)
#then run ANOVA with type III SS 
QIIMEdm_3wayANOVA_model_typeIII <- Anova(QIIMEdm_3wayANOVA_model, type = "III")
#Starch most descriptive, WG also descriptive, both interact with each other and with fat. May want to explore further 
ggplot(QIIMEdm_melted_3factors, aes(x = Group, y = value))+
  geom_boxplot(fill = c("gray75", "gray75",
                       "darkorange", "darkorange",
                       "blue", "blue",
                       "green", "green"),
               col = "black", lwd = 1)+
  ylab("Bray-Curtis Distance to CON-MFD")+
  theme_bw()

ggplot(QIIMEdm_melted_3factors, aes(x = Group, y = value))+
  geom_violin(aes(fill = Group), col = "black", lwd = 1)+
  scale_fill_manual(values = c("gray75", "gray75",
                      "darkorange", "darkorange",
                      "blue", "blue",
                      "green", "green"))+
  ylab("Bray-Curtis Distance to CON-MFD")+
  theme_bw()

#ALPHA DIVERSITY PLOTS ####
USDA_faithPD <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/faith-pd-USDA-alpha-div_filtered.csv",
         stringsAsFactors = FALSE, 
         header = TRUE,
         row.names = 1)
#already in the correct order so just bind the columns 
USDA_faithPD_w_diet <- cbind(USDA_faithPD, asv_mapfile_sub_ordered_matched$Code2)
names(USDA_faithPD_w_diet)[2] <- "Diet"

faithPD_plot <- ggplot(USDA_faithPD_w_diet, aes(x = Diet, y = faith_pd))+
  geom_boxplot(col = c("gray30", "red", "blue", "darkorchid4"), 
               fill = c("gray", "pink", "cyan", "darkorchid1"),
               lwd = 1,
               width = 0.9) + 
#  geom_point(col = "black") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 13),
        axis.title.y = element_text(size = 14)) + 
  ylab("Faith PD") + xlab(NULL)

faithPD_plot <- ggplot(USDA_faithPD_w_diet, aes(x = Diet, y = faith_pd))+
  geom_violin(lwd = 0.25, width = 1, col = "gray20", aes(fill = Diet)) + 
  geom_boxplot(col = "black", width = 0.20, outlier.size = 1.5, lwd = 0.75) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.y = element_text(size = 12),
        legend.position = "none") + 
  ylab("Faith PD") + xlab(NULL)+
  scale_fill_manual(values = c("gray75", "darkorange", "blue", "green"))

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("usda qiime2 alpha div faith pd.pdf", width = 3, height = 3)
faithPD_plot
dev.off()

tiff("usda qiime2 alpha div faith pd.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 1200)
faithPD_plot
dev.off()


USDA_observedASVs <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/observed-asvs-USDA-alpha-div_filtered.csv",
                         stringsAsFactors = FALSE, 
                         header = TRUE,
                         row.names = 1)
#already in the correct order so just bind the columns 
USDA_observedASVs_w_diet <- cbind(USDA_observedASVs, asv_mapfile_sub_ordered_matched$Code2)
names(USDA_observedASVs_w_diet)[2] <- "Diet"

observedASVs_plot <- ggplot(USDA_observedASVs_w_diet, aes(x = Diet, y = observed_otus))+
  geom_boxplot(col = c("gray30", "red", "blue", "darkorchid4"),
               fill = c("gray", "pink", "cyan", "darkorchid1"),
               lwd = 1,
               width = 0.9) + 
  #  geom_point(col = "black") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 13),
        axis.title.y = element_text(size = 14)) + 
  ylab("Observed ASVs") + xlab(NULL)

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
tiff("usda qiime2 alpha div observed asvs.tiff", width = 3, height = 4, units = "in", compression = "lzw", res = 1200)
observedASVs_plot
dev.off()

pdf("usda qiime2 alpha div observed asvs.pdf", width = 3, height = 4)
observedASVs_plot
dev.off()

#3.1 Two way ANOVA to see how RS and WG affect alpha diversity ####
USDA_faithPD_w_diet$Starch <- ""
USDA_faithPD_w_diet$Starch[USDA_faithPD_w_diet$Diet == "CON" | USDA_faithPD_w_diet$Diet == "WG"] <- "No"
USDA_faithPD_w_diet$Starch[USDA_faithPD_w_diet$Diet == "RS" | USDA_faithPD_w_diet$Diet == "WG+RS"] <- "Yes"
USDA_faithPD_w_diet$WholeGrain <- ""
USDA_faithPD_w_diet$WholeGrain[USDA_faithPD_w_diet$Diet == "CON" | USDA_faithPD_w_diet$Diet == "RS"] <- "No"
USDA_faithPD_w_diet$WholeGrain[USDA_faithPD_w_diet$Diet == "WG" | USDA_faithPD_w_diet$Diet == "WG+RS"] <- "Yes"
USDA_faithPD_w_diet$Fat <- ""
USDA_faithPD_w_diet$Fat[grep("H[0-9]{1,2}$",rownames(USDA_faithPD_w_diet))] <- "HighFat"
USDA_faithPD_w_diet$Fat[grep("M[0-9]{1,2}$",rownames(USDA_faithPD_w_diet))] <- "MedFat"


faithPD_3wayANOVA <- aov(faith_pd ~ Starch + 
                           WholeGrain + 
                           Fat + 
                           Starch:WholeGrain + 
                           Starch:Fat + 
                           WholeGrain:Fat + 
                           Starch:WholeGrain:Fat,
                         data = USDA_faithPD_w_diet)
summary(faithPD_3wayANOVA)
plot(faithPD_3wayANOVA)
TukeyHSD(faithPD_3wayANOVA, which = "Starch:WholeGrain")

USDA_observedASVs_w_diet$Starch <- ""
USDA_observedASVs_w_diet$Starch[USDA_observedASVs_w_diet$Diet == "CON" | USDA_observedASVs_w_diet$Diet == "WG"] <- "No"
USDA_observedASVs_w_diet$Starch[USDA_observedASVs_w_diet$Diet == "RS" | USDA_observedASVs_w_diet$Diet == "WG+RS"] <- "Yes"
USDA_observedASVs_w_diet$WholeGrain <- ""
USDA_observedASVs_w_diet$WholeGrain[USDA_observedASVs_w_diet$Diet == "CON" | USDA_observedASVs_w_diet$Diet == "RS"] <- "No"
USDA_observedASVs_w_diet$WholeGrain[USDA_observedASVs_w_diet$Diet == "WG" | USDA_observedASVs_w_diet$Diet == "WG+RS"] <- "Yes"
USDA_observedASVs_w_diet$Fat <- ""
USDA_observedASVs_w_diet$Fat[grep("H[0-9]{1,2}$",rownames(USDA_observedASVs_w_diet))] <- "HighFat"
USDA_observedASVs_w_diet$Fat[grep("M[0-9]{1,2}$",rownames(USDA_observedASVs_w_diet))] <- "MedFat"

observedASVs_2wayANOVA <- aov(observed_otus ~ Starch + 
                                WholeGrain + 
                                Fat + 
                                Starch:WholeGrain + 
                                Starch:Fat + 
                                WholeGrain:Fat + 
                                Starch:WholeGrain:Fat, data = USDA_observedASVs_w_diet)
summary(observedASVs_2wayANOVA)
plot(observedASVs_2wayANOVA)
TukeyHSD(observedASVs_2wayANOVA, which = "Starch:WholeGrain")

#CORRELATION PLOTS BETWEEN PH & BUTYRATE-ASSOCIATED TAXA FROM GNEISS####
#4.1 S24-7  ####
nrow(asv_mapfile_sub)
S247_abundance <- tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames["S24.7"]
nrow(S247_abundance)

#match rowname order in mapfile
S247_abundance <- S247_abundance[match(asv_mapfile_sub$X.SampleID,
                                       rownames(S247_abundance)), ,drop = FALSE]
rownames(S247_abundance) == asv_mapfile_sub$X.SampleID

asv_mapfile_sub$Code2[asv_mapfile_sub$Code2 == "AC"] <- "CON"
asv_mapfile_sub$Code2[asv_mapfile_sub$Code2 == "HM260"] <- "RS"
asv_mapfile_sub$Code2[asv_mapfile_sub$Code2 == "HMWG"] <- "WG+RS"
asv_mapfile_sub$Code2[asv_mapfile_sub$Code2 == "WWGC"] <- "WG"
asv_mapfile_sub$Code2 <- factor(asv_mapfile_sub$Code2, levels = c("CON", "WG", "RS", "WG+RS"))
#add butyrate and cecal pH data
S247_abundance_metadata <- cbind(S247_abundance, asv_mapfile_sub$Butyrate, asv_mapfile_sub$Cecal_pH, asv_mapfile_sub$Code2)
names(S247_abundance_metadata)[2:4] <- c("Butyrate", "Cecal pH", "Diet")

ggplot(S247_abundance_metadata, aes(x = S24.7,
       y = Butyrate)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black", formula = y ~ poly(x,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("[Butyrate]") + xlab("S24-7 Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))
  

ggplot(S247_abundance_metadata, aes(x = S24.7,
                                    y = `Cecal pH`)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black", formula = y ~ poly(x,2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("Cecal pH") + xlab("S24-7 Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))

#4.2 R. bromii  ####
Rbromii_abundance <- tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames["Ruminococcus.bromii"]
nrow(Rbromii_abundance)

#match rowname order in mapfile
Rbromii_abundance <- Rbromii_abundance[match(asv_mapfile_sub$X.SampleID,
                                       rownames(Rbromii_abundance)), ,drop = FALSE]
rownames(Rbromii_abundance) == asv_mapfile_sub$X.SampleID

#add butyrate and cecal pH data
Rbromii_abundance_metadata <- cbind(Rbromii_abundance, asv_mapfile_sub$Butyrate, asv_mapfile_sub$Cecal_pH, asv_mapfile_sub$Code2)
names(Rbromii_abundance_metadata)[2:4] <- c("Butyrate", "Cecal pH", "Diet")

ggplot(Rbromii_abundance_metadata, aes(x = Ruminococcus.bromii,
                                       y = Butyrate)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black", formula = y ~ poly(x,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("[Butyrate]") + xlab("R. bromii Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))

ggplot(Rbromii_abundance_metadata, aes(x = Ruminococcus.bromii,
                                    y = `Cecal pH`)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black", formula = y ~ poly(x,2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("Cecal pH") + xlab("R. bromii Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))

#4.3 Oscillospira  ####
oscillo_abundance <- tasv_table_SumByTaxa_prop_subset_ordered_taxa_rownames["Oscillospira"]
nrow(oscillo_abundance)

#match rowname order in mapfile
oscillo_abundance <- oscillo_abundance[match(asv_mapfile_sub$X.SampleID,
                                             rownames(oscillo_abundance)), ,drop = FALSE]
rownames(oscillo_abundance) == asv_mapfile_sub$X.SampleID

#add butyrate and cecal pH data
oscillo_abundance_metadata <- cbind(oscillo_abundance, asv_mapfile_sub$Butyrate, asv_mapfile_sub$Cecal_pH, asv_mapfile_sub$Code2)
names(oscillo_abundance_metadata)[2:4] <- c("Butyrate", "Cecal pH", "Diet")

ggplot(oscillo_abundance_metadata, aes(x = Oscillospira,
                                       y = Butyrate)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black", formula = y ~ poly(x,2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("[Butyrate]") + xlab("Oscillospira Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))

ggplot(oscillo_abundance_metadata, aes(x = Oscillospira,
                                       y = `Cecal pH`)) +
  geom_point(aes(fill = Diet), col = "black", alpha = 0.75, size = 3, pch = 21)+
  geom_smooth(method="lm", col = "black")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  ylab("Cecal pH") + xlab("Oscillospira Abundance")+
  scale_fill_manual(values = c("azure4", "red", "blue", "purple"))



#MAKING ASV-level FIGURES FOR GNEISS MODEL####
ASV_level_table <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/rarefied_ASV_table_17500_USDA.csv",
                              header = F,
                              stringsAsFactors = F,
                              comment.char = "")

#remove extra header line plus the "taxonomy" column
ASV_level_table <- ASV_level_table[-1,-ncol(ASV_level_table)]
ASV_level_table[1,1] <- "ASV" #rename "#OTU ID" column name
names(ASV_level_table) <- ASV_level_table[1,]
ASV_level_table <- ASV_level_table[-1,]
nrow(ASV_level_table) #1819 features 

#5.1 cleaning up y2 balance file ####
#read in balances of interest 
y2_asvs_pH_assoc <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/y2 ASVs cecal ph gneiss.csv",
         header = T,
         stringsAsFactors = F)

#Step-wise addition of deepest taxonomic assignments from species to kingdom level
#create a new "taxonomy" column, add species-level assignment to it when length of the species-level character string > 4)
y2_asvs_pH_assoc$taxonomy <- ""
#if row has no family level, put order level in taxonomy column
y2_asvs_pH_assoc$taxonomy[nchar(y2_asvs_pH_assoc$X4) <= 4] <- y2_asvs_pH_assoc$X3[nchar(y2_asvs_pH_assoc$X4) <= 4]
#if row has no genus level but does have family level, put family level in taxonomy column
y2_asvs_pH_assoc$taxonomy[nchar(y2_asvs_pH_assoc$X5) <= 4 &
                            nchar(y2_asvs_pH_assoc$X4) > 4] <- y2_asvs_pH_assoc$X4[nchar(y2_asvs_pH_assoc$X5) <= 4 &
                                                                                     nchar(y2_asvs_pH_assoc$X4) > 4]
#if row has no species level but does have genus level, put genus level in taxonomy column
y2_asvs_pH_assoc$taxonomy[nchar(y2_asvs_pH_assoc$X6) <= 4 &
                            nchar(y2_asvs_pH_assoc$X5) > 4] <- y2_asvs_pH_assoc$X5[nchar(y2_asvs_pH_assoc$X6) <= 4 &
                                                                                      nchar(y2_asvs_pH_assoc$X5) > 4]
#if row has species level, concatenate genus and species level strings and add a space in between.
#also specifying not to do this when taxonomic levels are repeated (weird behavior in QIIME2)
y2_asvs_pH_assoc$taxonomy[nchar(y2_asvs_pH_assoc$X6) > 4 &
                            y2_asvs_pH_assoc$X6 != y2_asvs_pH_assoc$X5] <- paste(y2_asvs_pH_assoc$X5[nchar(y2_asvs_pH_assoc$X6) > 4 &
                                                                                                             y2_asvs_pH_assoc$X6 != y2_asvs_pH_assoc$X5],
                                                                               y2_asvs_pH_assoc$X6[nchar(y2_asvs_pH_assoc$X6) > 4 &
                                                                                                     y2_asvs_pH_assoc$X6 != y2_asvs_pH_assoc$X5])
#adding in the repeated taxonomy levels 
y2_asvs_pH_assoc$taxonomy[y2_asvs_pH_assoc$X6 == y2_asvs_pH_assoc$X5] <- y2_asvs_pH_assoc$X6[y2_asvs_pH_assoc$X6 == y2_asvs_pH_assoc$X5]
#View(y2_asvs_pH_assoc)
#cleaning up prefixes
y2_asvs_pH_assoc$taxonomy <- gsub("[a-z]__", "", y2_asvs_pH_assoc$taxonomy)
#blank space in front of taxa names, and double spaces in between genus and species names.  Cleaning up 
y2_asvs_pH_assoc$taxonomy <- gsub("^ ", "", y2_asvs_pH_assoc$taxonomy)
y2_asvs_pH_assoc$taxonomy <- gsub("  ", " ", y2_asvs_pH_assoc$taxonomy)


#5.2 subsetting ASV table, generating taxa plot ####
ASV_level_table_y2_only <- ASV_level_table[ASV_level_table$ASV %in% y2_asvs_pH_assoc$Feature.ID,]
y2_asvs_pH_assoc_ordered <- y2_asvs_pH_assoc[match(ASV_level_table_y2_only$ASV, y2_asvs_pH_assoc$Feature.ID), ]
ASV_level_table_y2_only$taxonomy <- y2_asvs_pH_assoc_ordered$taxonomy
#View(ASV_level_table_y2_only)
#convert to abundance of the overall dataset (rarefaction depth 17,500)
ASV_level_table_y2_only[,-c(1,ncol(ASV_level_table_y2_only))] <- lapply(ASV_level_table_y2_only[,-c(1,ncol(ASV_level_table_y2_only))], function(x) as.numeric(as.character(x)))
ASV_level_table_y2_only_proportion <- ASV_level_table_y2_only
ASV_level_table_y2_only_proportion[,-c(1,ncol(ASV_level_table_y2_only))] <- ASV_level_table_y2_only[,-c(1,ncol(ASV_level_table_y2_only))]/17500

ASV_level_table_y2_only_proportion_melted <- melt(ASV_level_table_y2_only_proportion)

ggplot(ASV_level_table_y2_only_proportion_melted,
       aes(x = variable, y = value)) +
  geom_bar(stat="identity",aes(fill = taxonomy), width = 0.9, col = "black")+
  scale_fill_manual(values = rainbow(21))+
  theme_bw()+
  ylab("Relative Proportion")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5,"cm"))+
  ylim(c(0,1))+
  guides(fill=guide_legend(title="Taxonomy"))+
  ggtitle("Balance 'y2'")


#I want a pie chart/bar chart summarizing the composition of this balance
View(ASV_level_table_y2_only)
#create an object with the ASV#, taxonomy, and frequency across the dataset
y2_balance_taxonomy_sums <- data.frame(cbind(ASV_level_table_y2_only$ASV,
      ASV_level_table_y2_only$taxonomy,
      rowSums(ASV_level_table_y2_only[,3:ncol(ASV_level_table_y2_only)-1])))
y2_balance_taxonomy_sums[,3] <- as.numeric(as.character(y2_balance_taxonomy_sums[,3]))
str(y2_balance_taxonomy_sums)
#calculate the sum for each unique taxonomic group
y2_balance_taxonomy_sums_combined <- aggregate(x = y2_balance_taxonomy_sums$X3,
          by = list(y2_balance_taxonomy_sums$X2),
          FUN = "sum")
#reorder based on taxa level sum
y2_balance_taxonomy_sums_combined <- y2_balance_taxonomy_sums_combined[order(y2_balance_taxonomy_sums_combined$x, decreasing = TRUE),]
View(y2_balance_taxonomy_sums_combined)
y2_balance_taxonomy_sums_combined$proportion <- y2_balance_taxonomy_sums_combined$x / sum(y2_balance_taxonomy_sums_combined$x)
#combine low-abundance taxa into "Others".  Looks like 5% is a good cutoff 
y2_balance_taxonomy_sums_combined_low_abundance <- y2_balance_taxonomy_sums_combined[-c(1:5),]

y2_balance_taxonomy_sums_combined_low_abundance <- c("Other", sum(y2_balance_taxonomy_sums_combined_low_abundance$x),
      sum(y2_balance_taxonomy_sums_combined_low_abundance$proportion))

piechart_object_y2 <- rbind(y2_balance_taxonomy_sums_combined[c(1:5),], y2_balance_taxonomy_sums_combined_low_abundance)
piechart_object_y2$proportion <- as.numeric(piechart_object_y2$proportion)
piechart_object_y2$percent <- round(piechart_object_y2$proportion * 100, 1)
#set piechart slice sizes as equal to the counts 
slices <- as.numeric(as.character(piechart_object_y2$x))
#build labels for each slice 
labels <- c(as.character(piechart_object_y2$Group.1))
labels[6] <- "Others"
labels <- paste(labels, piechart_object_y2$percent)
labels <- paste(labels, "%", sep = "")
labels <- gsub(" ", "\n", labels) #add newline breaks so labels don't run so long 
#establish color palette.  For "Others" set to a neutral color 
y2_color_palette <- rainbow(length(labels))
y2_color_palette[6] <- "white"

#base R function to make piecharts 
pie(slices,
    labels = labels,
    col = y2_color_palette,
    main = "y2 denominator")

#ggplot method - I may need to use this method as I cannot save pie plot as object 
#ggplot(piechart_object_y2, aes(x = "", y = percent, fill = Group.1))+
#  geom_bar(stat = "identity") + 
#  coord_polar("y", start = 0)

#B. uniformis label is too close and there's no simple way to adjust label distances.  Decided to add a newline character in front of it to pull the whole label down
labels[4] <- "\n\nBacteroides\nuniformis\n8.6%"
#S24-7 label needs to be a tad higher 
labels[1] <- "S24-7\n40.6%\n"

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("usda qiime2 y2 piechart.pdf", width = 2.5, height = 2.5)
par(mar=c(0,2,0,0))

pie(slices,
    labels = labels,
    col = y2_color_palette,
#    main = "y2 denominator",
    cex = 0.5)

dev.off()

tiff("usda qiime2 y2 piechart.tiff", width = 2.5, height = 2.5, units = "in", compression = "lzw", res = 1200)
par(mar=c(0,2,0,0))
pie(slices,
    labels = labels,
    col = y2_color_palette,
    #    main = "y2 denominator",
    cex = 0.5)
dev.off()


#5.3 cleaning up y4 balance file ####
#butyrate 
y4_asvs_butyrate_assoc <- read.csv("G:/My Drive/TRANSCEND/USDA/Qiime2/y4 ASVs butyrate gneiss.csv",
                                   header = T,
                                   stringsAsFactors = F)

#Step-wise addition of deepest taxonomic assignments from species to kingdom level
#create a new "taxonomy" column, add species-level assignment to it when length of the species-level character string > 4)
y4_asvs_butyrate_assoc$taxonomy <- ""
#if row has no family level, put order level in taxonomy column
y4_asvs_butyrate_assoc$taxonomy[nchar(y4_asvs_butyrate_assoc$X4) <= 5] <- y4_asvs_butyrate_assoc$X3[nchar(y4_asvs_butyrate_assoc$X4) <= 5]
#if row has no genus level but does have family level, put family level in taxonomy column
y4_asvs_butyrate_assoc$taxonomy[nchar(y4_asvs_butyrate_assoc$X5) <= 5 &
                            nchar(y4_asvs_butyrate_assoc$X4) > 5] <- y4_asvs_butyrate_assoc$X4[nchar(y4_asvs_butyrate_assoc$X5) <= 5 &
                                                                                                 nchar(y4_asvs_butyrate_assoc$X4) > 5]
#if row has no species level but does have genus level, put genus level in taxonomy column
y4_asvs_butyrate_assoc$taxonomy[nchar(y4_asvs_butyrate_assoc$X6) <= 5 &
                            nchar(y4_asvs_butyrate_assoc$X5) > 5] <- y4_asvs_butyrate_assoc$X5[nchar(y4_asvs_butyrate_assoc$X6) <= 5 &
                                                                                     nchar(y4_asvs_butyrate_assoc$X5) > 5]
#if row has species level, concatenate genus and species level strings and add a space in between.
#also specifying not to do this when taxonomic levels are repeated (weird behavior in QIIME2)
y4_asvs_butyrate_assoc$taxonomy[nchar(y4_asvs_butyrate_assoc$X6) > 5 &
                            y4_asvs_butyrate_assoc$X6 != y4_asvs_butyrate_assoc$X5] <- paste(y4_asvs_butyrate_assoc$X5[nchar(y4_asvs_butyrate_assoc$X6) > 5 &
                                                                                                       y4_asvs_butyrate_assoc$X6 != y4_asvs_butyrate_assoc$X5],
                                                                                 y4_asvs_butyrate_assoc$X6[nchar(y4_asvs_butyrate_assoc$X6) > 5 &
                                                                                                       y4_asvs_butyrate_assoc$X6 != y4_asvs_butyrate_assoc$X5])
#adding in the repeated taxonomy levels 
y4_asvs_butyrate_assoc$taxonomy[y4_asvs_butyrate_assoc$X6 == y4_asvs_butyrate_assoc$X5] <- y4_asvs_butyrate_assoc$X6[y4_asvs_butyrate_assoc$X6 == y4_asvs_butyrate_assoc$X5]

#View(y4_asvs_butyrate_assoc)
#cleaning up prefixes
y4_asvs_butyrate_assoc$taxonomy <- gsub("[a-z]__", "", y4_asvs_butyrate_assoc$taxonomy)
#blank space in front of taxa names, and double spaces in between genus and species names.  Cleaning up 
y4_asvs_butyrate_assoc$taxonomy <- gsub("^ ", "", y4_asvs_butyrate_assoc$taxonomy)
y4_asvs_butyrate_assoc$taxonomy <- gsub("  ", " ", y4_asvs_butyrate_assoc$taxonomy)
#View(y4_asvs_butyrate_assoc)

#5.4 subsetting ASV table, generating taxa plot ####
ASV_level_table_y4_only <- ASV_level_table[ASV_level_table$ASV %in% y4_asvs_butyrate_assoc$Feature.ID,]
y4_asvs_butyrate_assoc_ordered <- y4_asvs_butyrate_assoc[match(ASV_level_table_y4_only$ASV, y4_asvs_butyrate_assoc$Feature.ID), ]
ASV_level_table_y4_only$taxonomy <- y4_asvs_butyrate_assoc_ordered$taxonomy
#View(ASV_level_table_y4_only)
#convert to abundance of the overall dataset (rarefaction depth 17,500)
ASV_level_table_y4_only[,-c(1,ncol(ASV_level_table_y4_only))] <- lapply(ASV_level_table_y4_only[,-c(1,ncol(ASV_level_table_y4_only))], function(x) as.numeric(as.character(x)))
#I want a pie chart/bar chart summarizing the composition of this balance
View(ASV_level_table_y4_only)
#create an object with the ASV#, taxonomy, and frequency across the dataset
y4_balance_taxonomy_sums <- data.frame(cbind(ASV_level_table_y4_only$ASV,
                                             ASV_level_table_y4_only$taxonomy,
                                             rowSums(ASV_level_table_y4_only[,3:ncol(ASV_level_table_y4_only)-1])))
y4_balance_taxonomy_sums[,3] <- as.numeric(as.character(y4_balance_taxonomy_sums[,3]))
str(y4_balance_taxonomy_sums)
#calculate the sum for each unique taxonomic group
y4_balance_taxonomy_sums_combined <- aggregate(x = y4_balance_taxonomy_sums$X3,
                                               by = list(y4_balance_taxonomy_sums$X2),
                                               FUN = "sum")
#reorder based on taxa level sum
y4_balance_taxonomy_sums_combined <- y4_balance_taxonomy_sums_combined[order(y4_balance_taxonomy_sums_combined$x, decreasing = TRUE),]
View(y4_balance_taxonomy_sums_combined)
y4_balance_taxonomy_sums_combined$proportion <- y4_balance_taxonomy_sums_combined$x / sum(y4_balance_taxonomy_sums_combined$x)
#combine low-abundance taxa into "Others".  Looks like 5% is a good cutoff 
y4_balance_taxonomy_sums_combined_low_abundance <- y4_balance_taxonomy_sums_combined[-c(1:5),]

y4_balance_taxonomy_sums_combined_low_abundance <- c("Other", sum(y4_balance_taxonomy_sums_combined_low_abundance$x),
                                                     sum(y4_balance_taxonomy_sums_combined_low_abundance$proportion))

piechart_object_y4 <- rbind(y4_balance_taxonomy_sums_combined[c(1:5),], y4_balance_taxonomy_sums_combined_low_abundance)
piechart_object_y4$proportion <- as.numeric(piechart_object_y4$proportion)
piechart_object_y4$percent <- round(piechart_object_y4$proportion * 100, 1)
#set piechart slice sizes as equal to the counts 
slices <- as.numeric(as.character(piechart_object_y4$x))
#build labels for each slice 
labels <- c(as.character(piechart_object_y4$Group.1))
labels[6] <- "Others"
labels <- paste(labels, piechart_object_y4$percent)
labels <- paste(labels, "%", sep = "")
labels <- gsub(" ", "\n", labels) #add newline breaks so labels don't run so long 
#establish color palette.  For "Others" set to a neutral color. Make sure any overlapping taxa from the other pie chart have the same color 
y4_color_palette <- rainbow(length(labels))
y4_color_palette[1] <- "purple"
y4_color_palette[2] <- "orange"
y4_color_palette[3] <- "pink"
y4_color_palette[4] <- "#FF0000FF"
y4_color_palette[5] <- "gray70"
y4_color_palette[6] <- "white"

#base R function to make piecharts 
pie(slices,
    labels = labels,
    col = y4_color_palette,
    main = "y4 denominator")

#ggplot method - I may need to use this method as I cannot save pie plot as object 
#ggplot(piechart_object_y4, aes(x = "", y = percent, fill = Group.1))+
#  geom_bar(stat = "identity") + 
#  coord_polar("y", start = 0)

#Clostiridiales label should be shifted up
labels[1] <- "Clostridiales\n32.9%\n"
#S24-7 label needs to be a tad lower 
labels[4] <- "\n\nS24-7\n11.1%\n"

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("usda qiime2 y4 piechart.pdf", width = 2.5, height = 2.5)
par(mar=c(0,2,0,0))

piechart <- pie(slices,
    labels = labels,
    col = y4_color_palette,
    #    main = "y4 denominator",
    cex = 0.5)

dev.off()


tiff("usda qiime2 y4 piechart.tiff", width = 2.5, height = 2.5, units = "in", compression = "lzw", res = 1200)
par(mar=c(0,2,0,0))
pie(slices,
    labels = labels,
    col = y4_color_palette,
    #    main = "y4 denominator",
    cex = 0.5)
dev.off()




#5.3 making ASV-cecal pH plot ####
#add cecal pH data to the melted table by matching the sample names across the table and map file 
ASV_level_table_y2_only_proportion[,-c(1,ncol(ASV_level_table_y2_only))] <- ASV_level_table_y2_only[,-c(1,ncol(ASV_level_table_y2_only))]/17500

#sum of all y2 denominator taxa 
y2_balance_props_by_sample <- colSums(ASV_level_table_y2_only_proportion[,3:ncol(ASV_level_table_y2_only_proportion)-1])
y2_balance_props_by_sample <- melt(y2_balance_props_by_sample)

#explicitly adding a column with sample names to make merging easier 
y2_balance_props_by_sample$sample <- rownames(y2_balance_props_by_sample)
y2_balance_props_by_sample_w_metadata <- merge(y2_balance_props_by_sample,
                                                        asv_mapfile_sub[,c("X.SampleID","Cecal_pH"),drop=F],
                                                        by.x = "sample",
                                                        by.y = "X.SampleID")
#WWGCH9 missing pH values so removing from table 
y2_balance_props_by_sample_w_metadata <- y2_balance_props_by_sample_w_metadata[
  y2_balance_props_by_sample_w_metadata$sample != "KEEWWGCH9",
  ]

#View(y2_balance_props_by_sample_w_metadata)

#simple linear regression of dnominator abundance and cecla pH 
y2_pH_lm <- lm(value ~ Cecal_pH, data = y2_balance_props_by_sample_w_metadata)
summary(y2_pH_lm)
names(summary(y2_pH_lm))
y2_Fstat <- round(summary(y2_pH_lm)$fstatistic["value"], 2)
y2_adjR <- round(summary(y2_pH_lm)$adj.r.squared, 2)
y2_pval <- formatC(summary(y2_pH_lm)$coefficients[2,4], format = "e", digits = 2) #formats p-val output to only have two decimal places 
#create annotation text layer for stats results 
y2_stats_annotation <- c(paste0("F ==", y2_Fstat),
  paste0("Adj.~italic(R)^2 ==", y2_adjR),
  paste0("italic(p) ==", y2_pval))

#set parse = TRUE to honor the formatting annotations in the annotation object 
y2_pH_ggplot <- ggplot(data = y2_balance_props_by_sample_w_metadata,
       aes(x = Cecal_pH, y = value)) +
  geom_smooth(method = "lm", se = TRUE)+
  geom_point(col = "black")+
    theme_bw()+
  ylab("y2 Denominator\nRelative Proportion")+xlab("Cecal pH")+
  ylim(0,1)+
  theme(axis.text = element_text(colour = "black"),
        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),
        panel.grid = element_blank())+
  annotate("text",
           x = 8.5, y = c(0.95, 0.80, 0.65),
           label = y2_stats_annotation, parse = TRUE,
           col = "blue", size = 4,
           hjust = 1)
  
setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("usda qiime2 y2 pH chart.pdf", width = 4.0, height = 2.5)
y2_pH_ggplot
dev.off()

tiff("usda qiime2 y2 pH chart.tiff", width = 4.0, height = 2.5, units = "in", compression = "lzw", res = 1200)
y2_pH_ggplot
dev.off()

#5.4 making ASV-cecal butyrate plot ####
#add cecal pH data to the melted table by matching the sample names across the table and map file 
ASV_level_table_y4_only_proportion <- ASV_level_table_y4_only
ASV_level_table_y4_only_proportion[,-c(1,ncol(ASV_level_table_y4_only))] <- ASV_level_table_y4_only[,-c(1,ncol(ASV_level_table_y4_only))]/17500

#sum of all y4 denominator taxa 
y4_balance_props_by_sample <- colSums(ASV_level_table_y4_only_proportion[,3:ncol(ASV_level_table_y4_only_proportion)-1])
y4_balance_props_by_sample <- melt(y4_balance_props_by_sample)

#explicitly adding a column with sample names to make merging easier 
y4_balance_props_by_sample$sample <- rownames(y4_balance_props_by_sample)
y4_balance_props_by_sample_w_metadata <- merge(y4_balance_props_by_sample,
                                               asv_mapfile_sub[,c("X.SampleID","Butyrate"),drop=F],
                                               by.x = "sample",
                                               by.y = "X.SampleID")
#WWGCH9 missing pH values so removing from table 
y4_balance_props_by_sample_w_metadata <- y4_balance_props_by_sample_w_metadata[
  y4_balance_props_by_sample_w_metadata$sample != "KEEWWGCH9",
  ]

#View(y4_balance_props_by_sample_w_metadata)

y4_pH_ggplot <- ggplot(data = y4_balance_props_by_sample_w_metadata,
                       aes(x = Butyrate, y = value)) +
  geom_smooth(method = "lm", se = TRUE)+
  geom_point(col = "black")+
  theme_bw()+
  ylab("y4 Numerator\nRelative Proportion")+xlab("Cecal butyrate conc.")+
  ylim(0,0.5)+
  theme(axis.text = element_text(colour = "black"),
        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))

setwd("G:/My Drive/TRANSCEND/USDA/Qiime2/")
pdf("usda qiime2 y4 pH chart.pdf", width = 4.0, height = 2.5)
y4_pH_ggplot
dev.off()

tiff("usda qiime2 y4 pH chart.tiff", width = 4.0, height = 2.5, units = "in", compression = "lzw", res = 1200)
y4_pH_ggplot
dev.off()


#Save Workspace  ####
#I'm using R 3.4 so use version 2. Version 3.5 is when version 3 was used 
save.image(file = "G:/My Drive/TRANSCEND/USDA/Qiime2/USDA QIIME2 environment v2.RData",
           version = 2)
