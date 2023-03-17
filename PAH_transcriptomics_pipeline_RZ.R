if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2",force = TRUE)
BiocManager::install("biomaRt",force = TRUE)
library(DESeq2)
library(biomaRt)
library(readxl)
library(tidyverse)

##########################################
###       read in clinical data        ###
##########################################
subject_summary <- read_excel("E:/capstone/data/UNWA-01-20VW Data Table - 091021.xlsx", 
                              sheet="Sample Meta Data")

subject_summary <- as.data.frame(subject_summary)
colnames(subject_summary)

colnames(subject_summary)[1] <- "parent_sample_name"
subject_summary$idno <- as.numeric(substr(subject_summary$CLIENT_SAMPLE_ID, 2, 5))
subject_summary$ETIOLOGY3 <- subject_summary$ETIOLOGY
subject_summary$ETIOLOGY3[!(subject_summary$ETIOLOGY %in% c("CTEPH", "Group235"))] <- "PAH" # adding a new ETIOLOGY variable with 3 levels

subject_survey <- read_excel("E:/capstone/data/UNWA-01-20VW Data Table - 091021.xlsx", 
                             sheet="Master Clinical Data")
subject_survey <- as.data.frame(subject_survey)
subject_survey <- subset(subject_survey, select=c(idno, echo_rvdia, echo_tapse, labs_ntprobnp_RESEARCH, 
                                                  dth_tt_use_2020, visit_6mwdist, visit_nyha, race, visit_height, visit_weight,
                                                  reveal_score, reveal_severity, echo_rvdil, echo_rvdilsev, 
                                                  echo_rvfx, echo_rvfxsev, rhc_pvr_use))
subject_survey$visit_height[!is.na(subject_survey$visit_height) & (subject_survey$visit_height < 100)] <- 165 # 1 person; correct data input error: 1.65 from original
subject_survey$bmi <- subject_survey$visit_weight/((subject_survey$visit_height/100)**2) # calculate BMI
subject_survey$bmi[is.na(subject_survey$bmi)] <- median(subject_survey$bmi, na.rm=T) # impute a few missing BMIs with median values
subject_survey$visit_6mwdist_meter <- subject_survey$visit_6mwdist*0.3048 # change visit_6mwdist (feet) to meters
subject_survey$race3 <- NA
subject_survey$race3[subject_survey$race=="White"] <- "white"
subject_survey$race3[(subject_survey$race!="White") & !(is.na(subject_survey$race))] <- "non-white"
subject_survey$labs_ntprobnp_RESEARCH[which.max(subject_survey$labs_ntprobnp_RESEARCH)] <- 4736 # correct one outlier (confirmed with raw data)

# echo rv diameter (two criteria)
subject_survey$echo_rvdia_ge5 <- (subject_survey$echo_rvdia >= 5.0)

subject_survey$echo_rvdia_ge4.2 <- (subject_survey$echo_rvdia >= 4.2)

# echo rv dilation + severity as a measure for RV
subject_survey$echo_rvdilsev_new <-
  ifelse(
    is.na(subject_survey$echo_rvdilsev) |
      is.na(subject_survey$echo_rvdil), 0,
    ifelse(
      (subject_survey$echo_rvdilsev=="Severe") &
        (subject_survey$echo_rvdil=="Yes"), 1, 0
    )
  )
table(subject_survey$echo_rvdilsev_new, exclude=NULL)
# 0   1 
# 105  35 

# echo rv function + severity as a measure for RV
subject_survey$echo_rvfxsev_new <-
  ifelse(
    is.na(subject_survey$echo_rvfxsev) |
      is.na(subject_survey$echo_rvfx), 0,
    ifelse(
      (subject_survey$echo_rvfxsev=="Severe") &
        (subject_survey$echo_rvfx=="Yes"), 1, 0
    )
  )
table(subject_survey$echo_rvfxsev_new, exclude=NULL)
# 0   1 
# 127  13

# among PVR>=7, rvdia>=5 vs rvdia<4.2
subject_survey$echo_rvdia_pvr7 <- NA
subject_survey$echo_rvdia_pvr7[(subject_survey$rhc_pvr_use >= 7) & (subject_survey$echo_rvdia >= 5)] <- 1
subject_survey$echo_rvdia_pvr7[(subject_survey$rhc_pvr_use >= 7) & (subject_survey$echo_rvdia < 4.2)] <- 0
sum(is.na(subject_survey$rhc_pvr_use) | is.na(subject_survey$echo_rvdia)) # 61
table(subject_survey$echo_rvdia_pvr7, exclude=NULL)
# 0    1 <NA> 
# 13   17  110 

# tapse
subject_survey$echo_tapse_lt2 <- (subject_survey$echo_tapse < 2.0)

# NYHA (1&2, 3&4)
subject_survey$visit_nyha_34 <- (subject_survey$visit_nyha %in% c("III: Symptoms with mild exertion", "IV: Symptoms at rest"))
subject_survey$visit_nyha_34[is.na(subject_survey$visit_nyha)] <- NA

# NT_proBNP
subject_survey$ntprobnp_cat <- NA
subject_survey$ntprobnp_cat[subject_survey$labs_ntprobnp_RESEARCH < 300] <- 1
subject_survey$ntprobnp_cat[(subject_survey$labs_ntprobnp_RESEARCH >= 300) & (subject_survey$labs_ntprobnp_RESEARCH < 1400)] <- 2
subject_survey$ntprobnp_cat[subject_survey$labs_ntprobnp_RESEARCH >= 1400] <- 3

# 6MWT
subject_survey$visit_6mwdist_meter_cat <- NA
subject_survey$visit_6mwdist_meter_cat[subject_survey$visit_6mwdist_meter < 165] <- 1
subject_survey$visit_6mwdist_meter_cat[(subject_survey$visit_6mwdist_meter >= 165) & (subject_survey$visit_6mwdist_meter < 320)] <- 2
subject_survey$visit_6mwdist_meter_cat[(subject_survey$visit_6mwdist_meter >= 320) & (subject_survey$visit_6mwdist_meter < 440)] <- 3
subject_survey$visit_6mwdist_meter_cat[subject_survey$visit_6mwdist_meter >= 440] <- 4

dim(subject_survey) # 140 29

# merge subject_survey with subject_summary (n_non_metabo is determined by # of columns in subject_summary)

subject_summary <- merge(subject_summary, subject_survey, by="idno")
dim(subject_summary) # 140 58
colnames(subject_summary)

## Only keep those with PAH???

subject_summary <- subject_summary %>% filter(ETIOLOGY3=="PAH")

##########################################
###       read in RNAseq data        ###
##########################################
dds<-read.delim('E:/capstone/data/All data raw counts for Jo.txt') 

colnames(dds)[41:67] <- substr(colnames(dds)[41:67], 2,50) 

dds<-dds[,c(-65,-66,-67)] 



key<-read_excel("E:/capstone/data/All Conditions key for Jo.xlsx") 
key<-data.frame(key[c(-64,-65,-66),c(2,5)]) #63 subjects 

#echo_rvdilsev_new with adjust variables: age gender bmi

subject_survey_sub<-subject_summary[,c("idno","ALIVE_AT_3_YEAR","AGE","GENDER","bmi","ETIOLOGY")] 
subject_survey_sub$ETIOLOGY <- as.factor(subject_survey_sub$ETIOLOGY)
key_id<-as.numeric(key[,2]) 

subject_survey_sub<-filter(subject_survey_sub,idno %in% key_id) 






data1<-merge(subject_survey_sub,key,by.x='idno',by.y='record_id...5') 
data1$RNAseq_ID <- gsub('\\.', '-', data1$RNAseq_ID) 

index <- duplicated(dds$Gene_Symbol) #13131 and 13134重复:1-mar, 2-mar 

which(index=="TRUE") 

dds <- dds[c(-13129,-13130),] 
colnames(dds)<-gsub('\\.', '-', colnames(dds)) 
rownames(dds) <- dds$Gene_Symbol


########################
###### biomart read ####
########################

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                host = "https://www.ensembl.org")

gb <-getBM(attributes=c("external_gene_name","gene_biotype"),filters="external_gene_name",
          values=rownames(dds), mart=mart)


gb <- gb[!gb$gene_biotype %in% c("unprocessed_pseudogene","transcribed_unprocessed_pseudogene",
                                  "processed_pseudogene"),]
## 22871
length(unique(gb$external_gene_name))
## 22837

## filter out pseudogene
dds <- dds %>% filter(Gene_Symbol %in% gb$external_gene_name)

#countData: 以样品为列名,基因为行名的表达矩阵
dds <- dds[,-1]
dds <- dds[,data1$RNAseq_ID]


#colData的形式: 以样品作为行名,列是样品对应的分组类型
coldata <- data.frame(row.names = data1$RNAseq_ID, condition=as.factor(data1$ALIVE_AT_3_YEAR),age=data1$AGE,gender=as.factor(data1$GENDER),bmi=data1$bmi,ETIOLOGY=data1$ETIOLOGY) 
## DESeq2 object
data <- DESeqDataSetFromMatrix(countData=dds, colData=coldata, design=~age+gender+bmi+ETIOLOGY+condition) 

#prefiltering 

keep <- rowSums(counts(data)) >= 10 

data <- data[keep,] 



#To perform the median of ratios method of normalization,  

#DESeq2 has a single estimateSizeFactors() function that will generate size factors. 

data <- estimateSizeFactors(data) 

sizeFactors(data) 

normalized_counts <- counts(data, normalized=TRUE) 
normalized_dat <- as.data.frame(normalized_counts)
Gene_Symbol <- rownames(normalized_dat)
normalized_dat <- cbind(Gene_Symbol,normalized_dat)
rownames(normalized_dat) <- NULL
write.table(normalized_dat, file="E:/capstone/data/normalized_dat.txt", sep="\t", quote=F,row.names = F) 

#NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that perform differential expression analysis that use the negative binomial model. 

######################################
###### Differential Analysis #########
######################################

ddm <- DESeq(data)
res <- results(ddm)
res <- res[!is.na(res$padj),] 
ALIVE_AT_3_YEAR.0.01 <- res[res$padj<0.01,] ## 174 
## res[res$padj<0.05,] 318
summary(ALIVE_AT_3_YEAR.0.01)

######################################
############# Heat map ###############
######################################

library("pheatmap")
rld <- rlogTransformation(ddm, blind=FALSE) 
exprMatrix_rlog <- assay(ddm) 

diff_gene.ALIVE_AT_3_YEAR <- subset(ALIVE_AT_3_YEAR.0.01, padj < 0.01 & (log2FoldChange > 1 | log2FoldChange < -1))

diff_data.ALIVE_AT_3_YEAR <- as.matrix(exprMatrix_rlog[rownames(diff_gene.ALIVE_AT_3_YEAR[order(diff_gene.ALIVE_AT_3_YEAR$padj)[1:40],]), ]) 
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
df <- as.data.frame(colData(ddm)[,c("condition","gender")])
df <- df[order(df$condition),]
diff_data.ALIVE_AT_3_YEAR <- diff_data.ALIVE_AT_3_YEAR[,rownames(df)]
pheatmap::pheatmap(diff_data.ALIVE_AT_3_YEAR,show_rownames = F,show_colnames = F,annotation_col = df, scale="row",cluster_cols = F)
