library(DESeq2)
library(openxlsx)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library("EnhancedVolcano")
library('pheatmap')
library('RColorBrewer')

files <- list.files(pattern = '.txt')    
for(i in files) {
  x <-read_delim(i, delim = "\t", escape_double = FALSE, 
                 col_types = cols(Chr = col_skip(), Start = col_skip(), 
                                  End = col_skip(), Strand = col_skip()),
                 trim_ws = TRUE, skip = 1)
  assign(i,x)  
}
Pub_v<- featureCounts_GRCm39.107_PD.txt %>% 
  rename_at(vars(contains('Final_GRCm39.107_')), list( ~ gsub('Final_GRCm39.107_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))

Pub_v_ann= getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                            filters = c('ensembl_gene_id') ,
                            values = Pub_v$Geneid,
                            mart = ensembl107)
							
							
Pub_v_wgn <- Rat_surge %>%
  dplyr::select(!Length) %>% 
  data.frame(row.names = 1)

Pub_v_cpm = as.data.frame(cpm(Pub_v_wgn))
colnames(Pub_v_cpm) = paste0('cpm_', colnames(Pub_v_cpm))
Pub_v_cpm <- rownames_to_column(Pub_v_cpm,var = "Geneid")

Pub_v_2 <- Pub_v%>% 
  inner_join(Pub_v_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Pub_v_cpm, by = "Geneid") %>%
  dplyr::select(Geneid, external_gene_name, description, everything()) %>% 
  mutate(round(.[, 4:30], digit = 2))


rawCounts <- data.frame(Pub_v2)
datapub <- data.frame(rawCounts, row.names = 1)  
desC <- rawCounts %>% select(1,2,3)
row.names(rawCounts) <- rawCounts$Geneid
#select raw reads of P22, Metestrous & Proestrous samples                       
rawCounts <- rawCounts %>% select(10:12,17:19,21:23)
colnames(rawCounts) = paste0('raw_',colnames(rawCounts))
rawCounts$Geneid <- row.names(rawCounts)
countS <- right_join(desC, rawCounts)
datapub <- datapub[-c(1,2)]
#select raw reads of P22, Metestrous & Proestrous samples                       
datapub <- datapub %>% select(7:9,14:16,18:20)
colnames(datapub)  <- c("P22_1","P22_2","P22_3","M_1","M_2","M_3","P_1","P_2","P_3")
infopub = data.frame(smpl=colnames(datapub))
rownames(infopub) = infopub$smpl
infopub$type = gsub('s', '', matrix(unlist(strsplit(infopub$smpl, '_')), nc=2, byrow=T)[,1])
infopub$type = factor(infopub$type)
infopub$type = relevel(infopub$type, ref = "P22")
  datap1 <- datapub[,1:6]
  datap2 <- datapub[,c(1:3,7:9)]
  infop1 <- infopub[1:6,]
  infop2 <- infopub[c(1:3,7:9),]
ddspub<- DESeqDataSetFromMatrix(countData =datapub, colData = infopub, design = ~type)
keep <- rowSums(counts(ddspub)> 5) >= 2
ddspub <- ddspub[keep,]
rld <- rlog(ddspub, blind = FALSE)
ddspub<- DESeqDataSetFromMatrix(countData =datap1, colData = infop1, design = ~type)
keep <- rowSums(counts(ddspub)> 5) >= 2
ddspub <- ddspub[keep,]
ddsDEpub <- DESeq(ddspub)    
vresTC = tibble(Geneid=rownames(datapub)) 
normCounts <- counts(ddsDEpub, normalized = T)
normCounts <- data.frame(normCounts)
colnames(normCounts) = paste0('norm_',colnames(normCounts))
normCounts$Geneid <- row.names(normCounts)
countS <- left_join(countS, normCounts)
resTC1 <- as.data.frame(results(ddsDEpub, name = "type_M_vs_P22", alpha=0.05))

ddspub<- DESeqDataSetFromMatrix(countData =datap2, colData = infop2, design = ~type)
keep <- rowSums(counts(ddspub)> 5) >= 2
ddspub <- ddspub[keep,]
ddsDEpub <- DESeq(ddspub)    
normCounts <- counts(ddsDEpub, normalized = T)
normCounts <- data.frame(normCounts[,4:6])
colnames(normCounts) = paste0('norm_',colnames(normCounts))
normCounts$Geneid <- row.names(normCounts)

countS <- left_join(countS, normCounts)
resTC2 <- as.data.frame(results(ddsDEpub, name = "type_P_vs_P22", alpha=0.05))

colnames(resTC1) = paste0("M_vs_P22",'_',colnames(resTC1))
colnames(resTC2) = paste0("P_vs_P22",'_',colnames(resTC2))
resTC1 = resTC1 %>% mutate(Geneid=rownames(.)) %>% as_tibble() %>% select(7,2,5,6)
resTC2 = resTC2 %>% mutate(Geneid=rownames(.)) %>% as_tibble() %>% select(7,2,5,6)
vresTC = left_join(vresTC, full_join(resTC1,resTC2))
all = left_join(countS, vresTC)

#KeggBrite ion channels [(ensmus)Geneid, external_gene_name]
Category_lists <- read_excel("~/R workplace/deseq2/Pubertás/Alapadatok/Category lists.xlsx", sheet = "Ionchannels")##########################################################
Ionchannels <- inner_join(all, Category_lists)
keep <- rowSums(Ionchannels[4:12] > 5) >= 2
Ionchannels <- Ionchannels[keep,]
write.xlsx(Ionchannels, "Supplementary Ion Channels KEGG Brite.xlsx")

Cacna1 <- Ionchannels[grepl("Cacna1", Ionchannels$external_gene_name),]
Scn <- Ionchannels[(grepl("Scn", Ionchannels$external_gene_name)&grepl("a", Ionchannels$external_gene_name)),]
Scn <- Scn %>% 
  filter(Scn$external_gene_name!= "Scn7a",  Scn$external_gene_name!= "Scn5a", Scn$external_gene_name!= "Scnn1a")
Hcn <- Ionchannels[grepl("Hcn", Ionchannels$external_gene_name),]
Kcn <- full_join(Ionchannels[grepl("Kcna", Ionchannels$external_gene_name),],
                 full_join(Ionchannels[grepl("Kcnb", Ionchannels$external_gene_name),],
                           full_join(Ionchannels[grepl("Kcnc", Ionchannels$external_gene_name),],
                                     full_join(Ionchannels[grepl("Kcnd", Ionchannels$external_gene_name),],
                                               full_join(Ionchannels[grepl("Kcnm", Ionchannels$external_gene_name),],
                                                         full_join(Ionchannels[grepl("Kcnn", Ionchannels$external_gene_name),],
                                                                   full_join(Ionchannels[grepl("Kcnq", Ionchannels$external_gene_name),],
                                                                             Ionchannels[grepl("Kcnt", Ionchannels$external_gene_name),])))))))
Kcn <- Kcn %>% 
  filter(Kcn$external_gene_name!= "Kcna5",
         Kcn$external_gene_name!= "Kcnab1",
         Kcn$external_gene_name!= "Kcnab2",
         Kcn$external_gene_name!= "Kcnab3",
         Kcn$external_gene_name!= "Kcnq4",
         Kcn$external_gene_name!= "Kcnt2")
SelectedI <- rbind(Cacna1, Scn, Hcn, Kcn)
#Heatmap
rlog <- as.data.frame(assay(rld)) 
rlog <- rlog %>% mutate(Geneid = row.names(rlog)) %>% select(10,1:9)
rlog <- inner_join(rlog,SelectedI[1:2])
row.names(rlog) <- rlog$external_gene_name
rm <- rlog %>% select(2:10)
test <- data.frame(rm)
pheatmap(test,					
         scale = "row", 
         color = colorRampPalette(rev(brewer.pal(n =11, name = "RdBu")))(200),
         cellwidth = 8,
         cellheight = 8,
         fontsize = 6,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 90,
         gaps_col = c(3,6),
         filename = "heatmap_SelectedIonChannels.pdf")
#Volcano
asd <- SelectedI %>% select(2,22,24,25,27)
asd[is.na(asd)] <- 1
row.names(asd) <- asd$external_gene_name
asd1 <- asd[grepl("Kcn", asd$external_gene_name) == FALSE,] %>% select(,2:5)
asd3 <- asd[grepl("Kcn", asd$external_gene_name),] %>% select(,2:5)
colnames(asd1)  <- c("Log2FC_P22vsM","padj_P22vsM","Log2FC_P22vsP","padj_P22vsP")
colnames(asd3)  <- colnames(asd1)
Gene1 <- row.names(asd1)
Gene3 <- row.names(asd3)
plot1 <- EnhancedVolcano(asd1,
                lab = Gene1,
                x = 'Log2FC_P22vsM',
                y = 'padj_P22vsM',
                xlim = c(-3, 3),
                ylim = c(0, 10),
                selectLab = Gene1,
                labSize = 3.0,
                shape= c(16),
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                title='Calcium, sodium and HCN channels (3w vs. M)', subtitle='',
                pCutoff = 0.05,
                FCcutoff = 1, 
                col = c('grey30', 'grey30', 'royalblue', 'red2'),
                pointSize = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                lengthConnectors = unit(0.01, "npc"),
                arrowheads = TRUE,
                boxedLabels = FALSE,
                max.overlaps = 20,
                legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-adj.~and~log[2]~FC)),
                legendPosition = 'top',
                directionConnectors = "both",
                legendLabSize = 10,
                legendIconSize = 2.0,
                caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05")
)
plot2 <- EnhancedVolcano(asd1,
                lab = Gene1,
                x = 'Log2FC_P22vsP',
                y = 'padj_P22vsP',
                xlim = c(-3, 3),
                selectLab = Gene1,
                labSize = 3.0,
                shape= c(16),
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 1,
                title='Calcium, sodium and HCN channels (3w vs. P)', subtitle='',
                pCutoff = 0.05,
                FCcutoff = 1, 
                col = c('grey30', 'grey30', 'royalblue', 'red2'),
                pointSize = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                lengthConnectors = unit(0.01, "npc"),
                arrowheads = TRUE,
                boxedLabels = FALSE,
                max.overlaps = 20,
                legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-adj.~and~log[2]~FC)),
                legendPosition = 'top',
                directionConnectors = "both",
                legendLabSize = 10,
                legendIconSize = 2.0,
                caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05")
)
plot3 <- EnhancedVolcano(asd3,
                         lab = Gene3,
                         x = 'Log2FC_P22vsM',
                         y = 'padj_P22vsM',
                         xlim = c(-3, 3),
                         ylim = c(0, 10),
                         selectLab = Gene3,
                         labSize = 3.0,
                         shape= c(16),
                         labCol = 'black',
                         labFace = 'bold',
                         colAlpha = 1,
                         title='Potassium channels (3w vs. M)', subtitle='',
                         pCutoff = 0.05,
                         FCcutoff = 1, 
                         col = c('grey30', 'grey30', 'royalblue', 'red2'),
                         pointSize = 1,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         lengthConnectors = unit(0.01, "npc"),
                         arrowheads = TRUE,
                         boxedLabels = FALSE,
                         max.overlaps = 20,
                         legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-adj.~and~log[2]~FC)),
                         legendPosition = 'top',
                         directionConnectors = "both",
                         legendLabSize = 10,
                         legendIconSize = 2.0,
                         caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05")
)
plot4 <- EnhancedVolcano(asd3,
                         lab = Gene3,
                         x = 'Log2FC_P22vsP',
                         y = 'padj_P22vsP',
                         xlim = c(-3, 3),
                         selectLab = Gene3,
                         labSize = 3.0,
                         shape= c(16),
                         labCol = 'black',
                         labFace = 'bold',
                         colAlpha = 1,
                         title='Potassium channels (3w vs. P)', subtitle='',
                         pCutoff = 0.05,
                         FCcutoff = 1, 
                         col = c('grey30', 'grey30', 'royalblue', 'red2'),
                         pointSize = 1,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         lengthConnectors = unit(0.01, "npc"),
                         arrowheads = TRUE,
                         boxedLabels = FALSE,
                         max.overlaps = 20,
                         legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-adj.~and~log[2]~FC)),
                         legendPosition = 'top',
                         directionConnectors = "both",
                         legendLabSize = 10,
                         legendIconSize = 2.0,
                         caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05")
)
plotlist <- list(plot1,plot2,plot3,plot4)
pdf("volcano selected ion channels.pdf")
for (x in 1:length(plotlist)) {
  print(plotlist[[x]])
  }
dev.off()
