BiocManager::install(version='3.14')
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library(maftools)

# access to clinical data
clin_query <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Clinical",
                       file.type = "xml")

GDCdownload(clin_query)
# get the clinic dataframe
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"
colnames(clinic)[colnames(clinic) == "gender"] <- "Biological_Sex"

# access to MAF
query_maf <- GDCquery(project = "TCGA-LUAD", 
                      data.category = "Simple Nucleotide Variation", 
                      data.type = "Masked Somatic Mutation", 
                      legacy = F)

GDCdownload(query_maf)
maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep, clinicalData = clinic, isTCGA = TRUE)


# access to the transcriptomic dataframe
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(query, method="api")
# sum_exp = GDCprepare(query)



# Clinic Data Analysis
# clean the NA values
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_followup, clinic$days_to_death)
# make the new column called the survival time and change it in years
clinic$survival_time = (clinic$days_to_death - clinic$days_to_birth) / 365
clinic$survival_time
# draw the boxplot for the survival distribution between male and female 
boxplot(clinic$survival_time ~ clinic$Biological_Sex,
        main = "Survival Data",
        xlab="sex",
        ylab = "Survival Time",)


# Find the survival time between female and male who get the lung cancer
library(survival)
library(survminer)

# clean the NA values
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_followup, clinic$days_to_death)
# create the death_event column
clinic$death_event = ifelse(clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
gender_fit <- surv_fit( surv_object ~ clinic$Biological_Sex, data = clinic )

# create the survplot
survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                      xlab = "Days",
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=15), # increase font sizes
        axis.text = element_text(size=15),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

p


# Analyze the survival time of different tobacco_smoking_history between male and female
# find the dataframe with tobacco_smoking_history == 1 (not smoking)
clinic_smoking = clinic[clinic$tobacco_smoking_history == 1,]

clinic_smoking$days_to_death = ifelse(is.na(clinic_smoking$days_to_death), clinic_smoking$days_to_last_followup, clinic_smoking$days_to_death)

clinic_smoking$days_to_death = as.numeric(clinic_smoking$days_to_death)
clinic$death_event = ifelse(clinic$vital_status == "Alive", 0, 1)
  
surv_object <- Surv(time = clinic_smoking$days_to_death, 
                    event = clinic_smoking$death_event) 

gender_fit <- surv_fit( surv_object ~ clinic_smoking$Biological_Sex, data = clinic_smoking)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      xlab = "Days",
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=15), # increase font sizes
        axis.text = element_text(size=15),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

p

# find the dataframe with tobacco_smoking_history == 2 (current smoker)
clinic_smoking = clinic[clinic$tobacco_smoking_history == 2,]

clinic_smoking$days_to_death = ifelse(is.na(clinic_smoking$days_to_death), clinic_smoking$days_to_last_followup, clinic_smoking$days_to_death)
clinic_smoking$days_to_death = as.numeric(clinic_smoking$days_to_death)
clinic_smoking$death_event = ifelse(clinic_smoking$vital_status == "Alive", 0, 1)

surv_object <- Surv(time = clinic_smoking$days_to_death, 
                    event = clinic_smoking$death_event) 

gender_fit <- surv_fit( surv_object ~ clinic_smoking$Biological_Sex, data = clinic_smoking)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      xlab = "Days",
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=15), # increase font sizes
        axis.text = element_text(size=15),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

p

# find the dataframe with tobacco_smoking_history == 3 (Current reformed smoker for > 15 years)
clinic_smoking = clinic[clinic$tobacco_smoking_history == 3,]

clinic_smoking$days_to_death = ifelse(is.na(clinic_smoking$days_to_death), clinic_smoking$days_to_last_followup, clinic_smoking$days_to_death)
clinic_smoking$days_to_death = as.numeric(clinic_smoking$days_to_death)
clinic_smoking$death_event = ifelse(clinic_smoking$vital_status == "Alive", 0, 1)

surv_object <- Surv(time = clinic_smoking$days_to_death, 
                    event = clinic_smoking$death_event) 

gender_fit <- surv_fit( surv_object ~ clinic_smoking$Biological_Sex, data = clinic_smoking)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      xlab = "Days",
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=15), # increase font sizes
        axis.text = element_text(size=15),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

p

# find the dataframe with tobacco_smoking_history == 4 (Current reformed smoker for ≤15 years)
clinic_smoking = clinic[clinic$tobacco_smoking_history == 4,]

clinic_smoking$days_to_death = ifelse(is.na(clinic_smoking$days_to_death), clinic_smoking$days_to_last_followup, clinic_smoking$days_to_death)
clinic_smoking$days_to_death = as.numeric(clinic_smoking$days_to_death)
clinic_smoking$death_event = ifelse(clinic_smoking$vital_status == "Alive", 0, 1)

surv_object <- Surv(time = clinic_smoking$days_to_death, 
                    event = clinic_smoking$death_event) 

gender_fit <- surv_fit( surv_object ~ clinic_smoking$gender, data = clinic_smoking)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      title = "Surv time with less than 15 years smokers between sex",
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

p


# MAF Analysis

# play with the number of genes with the "top" argument
oncoplot(maf = maf_object,
         top = 10) 

# Write to clinic
clinic = maf_object@clinical.data

# create the male vector
male_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "MALE"]
# create the female vector
female_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "FEMALE"]
# create the male submaf
male_maf = subsetMaf(maf = maf_object,
                     tsb = male_patients_ids)
# create the female submaf
female_maf = subsetMaf(maf = maf_object,
                       tsb = female_patients_ids)

# create the new column called smokingTime
clinic$smokingTime = ifelse(clinic$tobacco_smoking_history == 1, "non-smoker", "smoker")
# create the non-smoker vector
nonSmoker_patients_ids = clinic$Tumor_Sample_Barcode[clinic$smokingTime == "non-smoker"]
# create the smoker vector
smoker_patients_ids = clinic$Tumor_Sample_Barcode[clinic$smokingTime == "smoker"]
# create the submaf for nonsmoker in male 
nonSmoker_male_maf = subsetMaf(maf = male_maf,
                          tsb = nonSmoker_patients_ids)
# create the submaf for smoker in male 
smoker_male_maf = subsetMaf(maf = male_maf,
                          tsb = smoker_patients_ids)

# create the submaf for nonsmoker in female 
nonSmoker_female_maf = subsetMaf(maf = female_maf,
                               tsb = nonSmoker_patients_ids)
# create the submaf for smoker in female 
smoker_female_maf = subsetMaf(maf = female_maf,
                            tsb = smoker_patients_ids)

# The link between lung cancer-related genes and tobacco smoking
lollipopPlot(maf_object, gene = "KLF6")
lollipopPlot(maf_object, gene = "TERT")
lollipopPlot(maf_object, gene = "MSH5")
lollipopPlot(maf_object, gene = "GATA3")

# create the lollipopPlot2 for gene KLF6
lollipopPlot2(m1 = nonSmoker_male_maf, 
              m2 = smoker_male_maf, 
              m1_name= "nonsmoker in male",
              m2_name = "smoker in male",
              gene = "KLF6")

lollipopPlot2(m1 = nonSmoker_female_maf, 
              m2 = smoker_female_maf, 
              m1_name= "nonsmoker in female",
              m2_name = "smoker in female",
              gene = "KLF6")

# create the lollipopPlot2 for gene TERT
lollipopPlot2(m1 = nonSmoker_male_maf, 
              m2 = smoker_male_maf, 
              m1_name= "nonsmoker in male",
              m2_name = "smoker in male",
              gene = "TERT")

lollipopPlot2(m1 = nonSmoker_female_maf, 
              m2 = smoker_female_maf, 
              m1_name= "nonsmoker in female",
              m2_name = "smoker in female",
              gene = "TERT")

# create the lollipopPlot2 for gene MSH5
lollipopPlot2(m1 = nonSmoker_male_maf, 
              m2 = smoker_male_maf, 
              m1_name= "nonsmoker in male",
              m2_name = "smoker in male",
              gene = "MSH5")

lollipopPlot2(m1 = nonSmoker_female_maf, 
              m2 = smoker_female_maf, 
              m1_name= "nonsmoker in female",
              m2_name = "smoker in female",
              gene = "MSH5")

# create the lollipopPlot2 for gene GATA3
lollipopPlot2(m1 = nonSmoker_male_maf, 
              m2 = smoker_male_maf, 
              m1_name= "nonsmoker in male",
              m2_name = "smoker in male",
              gene = "GATA3")

lollipopPlot2(m1 = nonSmoker_female_maf, 
              m2 = smoker_female_maf, 
              m1_name= "nonsmoker in female",
              m2_name = "smoker in female",
              gene = "GATA3")


