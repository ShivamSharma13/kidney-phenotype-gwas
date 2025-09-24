## This repository has scripts used for analyses in the following publication: 

[African ancestry-enriched variants in the GATM gene are associated with elevated serum creatinine levels](https://www.medrxiv.org/content/10.1101/2025.03.07.25323581v1) on _MedRxiv (2025)._

## Overview of scripts and results

### **Step 1:** Phenotype QC and modeling
+ _01_PhenotypeQC.r_ - QC of serum creatinine and stage setup for GWAS.
+ _02_PhenotypeModeling.r_ - Simple statistical modeling of serum creatinine and covariates.
<img src="https://github.com/user-attachments/assets/01544945-b4ed-4ff5-af2d-222bf74067c3"  width=50% height=50%/>

### **Step 2:** Standard and ancestry-specific GWAS using SAIGE
+ _03_GWASWorkflow.Saige.wdl_ - WDL workflow (ideal for cromwell) for GWAS using SAIGE.
+ _04_AnalyzeGWASResults.r_ - Analyzing and visualizing (manhattan plots) GWAS results. 
+ _05_PerformAdmixtureMapping.r_ - Admixture mapping using bioconductor packages.
+ _06_AnalyzeTractorGWAS.r_ - Ancestry-specific GWAS using tractor.
<img src="https://github.com/user-attachments/assets/11798ca9-b9db-4d27-9438-42c5c77e293c"  width=50% height=50%/>

### **Step 3:** Clumping and fine mapping
+ _08_LD&FineMapping.r_ - LD clumping and variant fine-mapping using Susie.
<img src="https://github.com/user-attachments/assets/88acd462-921e-4720-91ce-efaaa0e8c422"  width=50% height=50%/>

### **Step 4:** Polygenic scores estimation and modeling
+ _07_ClumpingAndPRS.ipynb_ - Stepwise PRS estimation and modeling.
+ _10_AnnotateVariants.wdl_ - Variant effect annotation using SnpEff.
<img src="https://github.com/user-attachments/assets/7d420c38-8bf8-4031-9973-9e9aba5abdb9" width=50% height=50%/>


