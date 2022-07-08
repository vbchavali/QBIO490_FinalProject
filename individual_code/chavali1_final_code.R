if(!require(BiocManager)) install.packages("BiocManager")
if(!require(TCGAbiolinks)) BiocManager::install(version='3.14')
if(!require(BioinformaticsFMRP/TCGAbiolinksGUI.data)) BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
if(!require(BioinformaticsFMRP/TCGAbiolinks)) BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)

# downloading the query
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query, method="api")  
sum_exp = GDCprepare(query)

# downloading clinical data
clin_query <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Clinical",
                       file.type = "xml")

GDCdownload(clin_query)
# get the clinic dataframe
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"
# Changing the gender column name to Biological Sex
colnames(clinic)[colnames(clinic) == "gender"] <- "Biological_Sex"

# downloading MAF
query_maf <- GDCquery(project = "TCGA-LUAD", 
                      data.category = "Simple Nucleotide Variation", 
                      data.type = "Masked Somatic Mutation", 
                      legacy = F)

GDCdownload(query_maf)
maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep, clinicalData = clinic, isTCGA = TRUE)

# getting rid of NA values in gender column
bool_gender_na = is.na(colData(sum_exp)$gender)
sex_no_na = sum_exp$gender[!bool_gender_na]

# Getting rid of outliers and creating boxplot of CYP1A1 counts
"CYP1A1" %in% rowData(sum_exp)$gene_name
cyp1a1_ensembl_mask <- rowData(sum_exp)$gene_name == "CYP1A1"
ensembl_cyp1a1 <- rowData(sum_exp)$gene_id[cyp1a1_ensembl_mask]
cyp1a1_counts <- assays(sum_exp)$unstranded[ensembl_cyp1a1,!bool_gender_na]
cyp1a1_outliers <- boxplot.stats(cyp1a1_counts)$out
cyp1a1_out_indx <- which(cyp1a1_counts %in% cyp1a1_outliers)
cyp1a1_no_outliers <- cyp1a1_counts[-cyp1a1_out_indx]
sex_no_na1 <- sex_no_na[-cyp1a1_out_indx]

boxplot(cyp1a1_no_outliers~sex_no_na1, 
        main = "CYP1A1 Counts Without Outliers Between Male and Females",
        xlab = "Sex", 
        ylab = "CYP1A1 Counts")
cyp1a1_summary <- boxplot(cyp1a1_no_outliers~sex_no_na1)$stats

# Getting rid of outliers and creating boxplot of GSTM1 counts
"GSTM1" %in% rowData(sum_exp)$gene_name
gstm1_ensembl_mask <- rowData(sum_exp)$gene_name == "GSTM1"
ensembl_gstm1 <- rowData(sum_exp)$gene_id[gstm1_ensembl_mask]
gstm1_counts <- assays(sum_exp)$unstranded[ensembl_gstm1,]
gstm1_outliers <- boxplot.stats(gstm1_counts)$out
gstm1_out_indx <- which(gstm1_counts %in% gstm1_outliers)
gstm1_no_outliers <- gstm1_counts[-gstm1_out_indx]
sex_no_na2 <- sex_no_na[-gstm1_out_indx]

boxplot(gstm1_no_outliers~sex_no_na2, 
        main = "GSTM1 Counts Without Outliers Between Male and Females",
        xlab = "Sex", 
        ylab = "GSTM1 Counts")
gstm1_summary <- boxplot(gstm1_no_outliers~sex_no_na2)$stats

# Oncoplot & coOncoplot code
clinic = maf_object@clinical.data
male_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$Biological_Sex == "MALE"])
female_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$Biological_Sex == "FEMALE"])

male_maf = subsetMaf(maf = maf_object,
                      tsb = male_patients_ids)

female_maf = subsetMaf(maf = maf_object,
                     tsb = female_patients_ids)

library("ggplot2") 

oncoplot(maf = maf_object,
         top = 20) 

top_genes = c("TP53", "TTN", "MUC16", "CSMD3", "RYR2", "LRP1B", "ZFHX4", "USH2A",
              "KRAS", "XIRP2", "FLG", "SPTA1", "NAV3", "COL11A1", "ZNF536", "ANK2",
              "FAT3", "PCLO", "CSMD1", "PCDH15")
smoking_genes = c("GATA3", "MSH5", "TERT", "KLF6")

coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = "Male Patients", 
           m2Name = "Female Patients", 
           genes = top_genes)

coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = "Male Patients", 
           m2Name = "Female Patients", 
           genes = smoking_genes)

