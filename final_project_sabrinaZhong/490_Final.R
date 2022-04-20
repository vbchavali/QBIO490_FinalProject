BiocManager::install(version='3.14')
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query, method="api")  
sum_exp2 = GDCprepare(query)

bool_gender_na = is.na(colData(sum_exp2)$gender)
sex_no_na = sum_exp2$gender[!bool_gender_na]

head (colData(sum_exp2)$exposure_id)

dev.off()
"GSTM1" %in% rowData(sum_exp2)$gene_name
cyp1a1_ensembl_mask <- rowData(sum_exp2)$gene_name == "CYP1A1"
ensembl_cyp1a1 <- rowData(sum_exp2)$gene_id[cyp1a1_ensembl_mask]
cyp1a1_counts <- assays(sum_exp2)$unstranded[ensembl_cyp1a1,!bool_gender_na]
boxplot(cyp1a1_counts~sex_no_na, 
        cex.axis = 1,
        xlab = "Sex")


"CYP1A1" %in% rowData(sum_exp2)$gene_name
gstm1_ensembl_mask <- rowData(sum_exp2)$gene_name == "GSTM1"
ensembl_gstm1 <- rowData(sum_exp2)$gene_id[gstm1_ensembl_mask]
gstm1_counts <- assays(sum_exp2)$unstranded[ensembl_gstm1,]

boxplot(gstm1_counts~sex_no_na, 
        cex.axis = 1,
        xlab = "Sex")

"KLF6" %in% rowData(sum_exp2)$gene_name
klf6_ensembl_mask <- rowData(sum_exp2)$gene_name == "KLF6"
ensembl_klf6 <- rowData(sum_exp2)$gene_id[klf6_ensembl_mask]
klf6_counts <- assays(sum_exp2)$unstranded[ensembl_klf6,]

boxplot(klf6_counts~sex_no_na, 
        cex.axis = 1,
        main = "KLF6 Counts in Male vs. Female Patients",
        ylab = "KLF6 Gene Counts",
        xlab = "Sex")

"TERT" %in% rowData(sum_exp2)$gene_name
tert_ensembl_mask <- rowData(sum_exp2)$gene_name == "TERT"
ensembl_tert <- rowData(sum_exp2)$gene_id[tert_ensembl_mask]
tert_counts <- assays(sum_exp2)$unstranded[ensembl_tert,]

boxplot(tert_counts~sex_no_na, 
        cex.axis = 1,
        main = "TERT Counts in Male vs. Female Patients",
        ylab = "TERT Gene Counts",
        xlab = "Sex")

"GATA3" %in% rowData(sum_exp2)$gene_name
gata3_ensembl_mask <- rowData(sum_exp2)$gene_name == "GATA3"
ensembl_gata3 <- rowData(sum_exp2)$gene_id[gata3_ensembl_mask]
gata3_counts <- assays(sum_exp2)$unstranded[ensembl_gata3,]

boxplot(gata3_counts~sex_no_na, 
        cex.axis = 1,
        main = "GATA3 Counts in Male vs. Female Patients",
        ylab = "GATA3 Gene Counts",
        xlab = "Sex")

"MSH5" %in% rowData(sum_exp2)$gene_name
msh5_ensembl_mask <- rowData(sum_exp2)$gene_name == "MSH5"
ensembl_msh5 <- rowData(sum_exp2)$gene_id[msh5_ensembl_mask]
msh5_counts <- assays(sum_exp2)$unstranded[ensembl_msh5,]

boxplot(msh5_counts~sex_no_na, 
        cex.axis = 1,
        main = "MSH5 Counts in Male vs. Female Patients",
        ylab = "MSH5 Gene Counts",
        xlab = "Sex")

query_clin <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Clinical", 
                       file.type = "xml")
GDCdownload(query_clin)
clinic <- GDCprepare_clinic(query_clin, clinical.info = "patient")
colnames(clinic)[colnames(clinic)=="bcr_patient_barcode"] <-"Tumor_Sample_Barcode"

query_maf <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf)
maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep, clinicalData = clinic, isTCGA = TRUE)

head(assay(sum_exp2))

clinic = maf_object@clinical.data
male_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$gender == "MALE"])
female_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$gender == "FEMALE"])

male_maf = subsetMaf(maf = maf_object,
                      tsb = male_patients_ids)

female_maf = subsetMaf(maf = maf_object,
                     tsb = female_patients_ids)

library("ggplot2") 

oncoplot(maf = maf_object,
         top = 20) 

coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = "Male Patients", 
           m2Name = "Female Patients")

head(clinic)

male_patients_ids = c(clinic$gender == "male")

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients", 
              m2_name = "Female Patients",
              gene = "CYP1A1")

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients", 
              m2_name = "Female Patients",
              gene = "GSTM1")

