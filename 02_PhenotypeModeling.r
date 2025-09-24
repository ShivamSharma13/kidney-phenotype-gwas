rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')
library('data.table')
library('parallel')
library('car')
library('broom')
library('gridExtra')
library('tidyr')
library('scales')

options(repr.plot.width=16, repr.plot.height=10)

phenotype = as.data.frame(fread("Data/Saige.V2/CreatinineCohort.PhenotypeFile.tsv"))
rye = as.data.frame(fread("Data/Metadata/RyeEstimates.16From30PCS.UnrelatedParticipantsV7.Continental.7.Q"))

phenotype = merge(phenotype, rye, by = "SampleID")
phenotype = phenotype %>% mutate(AfricanFraction = African / 100)
mean(phenotype$African)
mean(phenotype$European)


dim(phenotype)
head(phenotype)

# Define the bins
breaks <- seq(0, 100, by = 2)  # Create breaks from 0 to 100 with a step of 2
labels <- paste(head(breaks, -1), tail(breaks, -1) , sep = "-")  # Create labels like "0-1", "2-3", etc.

# Create summary data frame with mean of C2 for each bin
admixreg_scr_binned = phenotype %>%
    mutate(AncestryBin = cut(African, breaks = breaks, labels = labels, include.lowest = TRUE)) %>%
    group_by(AncestryBin, Sex) %>%
    summarize(MeanCreatinineBin = mean(MedianCreatinine, na.rm = TRUE), .groups = 'drop') %>%
    mutate(AncestryBinLabel = as.numeric(sub("-.*", "", AncestryBin)) + 2)

head(admixreg_scr_binned)
tail(admixreg_scr_binned)



correlation_results <- admixreg_scr_binned %>%
  group_by(Sex) %>%
  summarize(Correlation = cor(MeanCreatinineBin, AncestryBinLabel, use = "complete.obs"),
           CorrelationSq = Correlation^2)

correlation_results

text_data <- data.frame(
  Sex = c("Female", "Male"),
  Label = c("R = 0.78", "R = 0.74")
)

g = ggplot(admixreg_scr_binned, aes(x = AncestryBinLabel, y = MeanCreatinineBin)) +
    geom_smooth(method='lm', formula= y~x, linewidth = 2, color = "darkblue")+
    geom_point(pch = 21, fill = "blue4", color = "black", size = 4) +
    scale_x_continuous(limits = c(0,100)) +
    labs(y = "Creatinine Levels (mg/dL)", x = "Percentage African Ancestry") +
    facet_wrap(~Sex) +
    geom_text(data = text_data, aes(x = -Inf, y = Inf, label = Label), 
               vjust = 1.5, hjust = -0.1, size = 10, color = "black") +
    theme_bw(base_size = 32)

ggsave(filename = "Data/Figures/AdmixtureRegression.Creatinine.png", g, device = "png", width = 16, height = 10)

g

model_scr = lm(formula = MedianCreatinine ~ AfricanFraction + Age + Sex, data = phenotype)
summary(model_scr)

#model_scr = lm(formula = MedianCreatinine ~ African + Age + Sex, data = phenotype)
#summary(model_scr)

dosages_15 = as.data.frame(fread("Data/Saige.V2/Dosage/extractedChr15.ACAF.GWASThreshold.Dosage.raw"))
head(dosages_15)
dosages = dosages_15 %>% select(IID, `chr15:45379909:C:T_C`, `chr15:45466280:A:G_A`)
names(dosages) = c("SampleID", "LeadSNP1", "LeadSNP3")

dosages = dosages %>% mutate(LeadSNP1 = 2 - LeadSNP1, LeadSNP3 = 2 - LeadSNP3)
phenotype = merge(phenotype, dosages %>% select(SampleID, LeadSNP1, LeadSNP3), 
                  by = "SampleID")

dim(phenotype)
head(phenotype)

phenotype %>% group_by(Sex) %>% summarize(mean(MedianCreatinine))

model_scr = lm(formula = MedianCreatinine ~ AfricanFraction + LeadSNP1 + Age + Sex, data = phenotype)
summary(model_scr)

model_scr = lm(formula = MedianCreatinine ~ AfricanFraction + LeadSNP3 + Age + Sex, data = phenotype)
summary(model_scr)

model_scr = lm(formula = MedianCreatinine ~ AfricanFraction + LeadSNP1 + LeadSNP3 + Age + Sex, data = phenotype)
summary(model_scr)

# Get variants for bunch of p-values.
cutoffs = c(0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.0000075, 0.000005, 0.000001, 0.0000005, 0.0000001, "RightQuadrant1", "RightQuadrant2", "RightQuadrant3")
cutoffs = c("RightQuadrant1", "RightQuadrant2", "RightQuadrant3")

results_ancestry = list()
results_prs = list()

for (cutoff_index in 1:length(cutoffs)) {
    prs = as.data.frame(fread(paste0("Data/Saige.V2/PGS/AfricanRegression/",cutoffs[cutoff_index],"/PRSice.best")))
    this_phenotype_prs = merge(phenotype, prs %>% select(IID, PRS), by.x = "SampleID", by.y = "IID")
    this_phenotype_prs = this_phenotype_prs %>% mutate(ScaledPRS = (PRS - min(PRS)) / (max(PRS) - min(PRS)))
    
    #Run the models.
    model_scr = lm(formula = MedianCreatinine ~ ScaledPRS + AfricanFraction + Age + Sex, data = this_phenotype_prs)
    summary_model_scr = summary(model_scr)
    #print(summary(model_scr))
    
    # Get the summary and tidy output and store effect size and beta in the results list
    summary_model_scr = tidy(summary_model_scr) %>% as.data.frame()
    results_ancestry[[cutoff_index]] = summary_model_scr %>%
                        filter(term %in% "AfricanFraction") %>%
                        select("estimate","std.error","statistic","p.value") %>%
                        mutate(Cutoff = paste0(cutoffs[cutoff_index]), Predictor = "AfricanAncestry")
    
    results_prs[[cutoff_index]] = summary_model_scr %>%
                        filter(term %in% "ScaledPRS") %>%
                        select("estimate","std.error","statistic","p.value") %>%
                        mutate(Cutoff = paste0(cutoffs[cutoff_index]), Predictor = "PRS")
    
    rm(this_phenotype_prs)
}

results_df_ancestry = do.call(rbind, results_ancestry)
names(results_df_ancestry) = c("Estimate","StdErr", "Statistic", "Pvalue", "Cutoff", "Predictor")

results_df_prs = do.call(rbind, results_prs)
names(results_df_prs) = c("Estimate","StdErr","Statistic","Pvalue", "Cutoff", "Predictor")

final_results = rbind(results_df_ancestry, results_df_prs)
final_results = final_results %>% mutate(Cutoff = factor(Cutoff, levels = c(0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.0000075, 0.000005, 0.000001, 0.0000005, 0.0000001, "RightQuadrant1", "RightQuadrant2", "RightQuadrant3")))



final_results 

g_beta = ggplot(final_results, aes(x = Cutoff, y = Estimate)) +
    geom_errorbar(aes(ymin = Estimate - StdErr, ymax = Estimate + StdErr),
                width = 0.2, color = "black") + 
    geom_point(size = 4, color = "darkmagenta") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "darkgreen") +
    facet_wrap(~Predictor, scales = "free_y") +
    labs(title = "Comparison of estimates",
       x = "-log10(Cutoffs) on GWAS summary stats (BetaAfr > 0)",
       y = "Estimate") +
    theme_bw(24) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

g_pv = ggplot(final_results, aes(x = Cutoff, y = log10(Pvalue))) +
    geom_point(size = 6, color = "darkmagenta") +  
    facet_wrap(~Predictor, scales = "free_y") +
    #scale_y_continuous(labels = scientific_format()) +
    labs(title = "Comparison of p-values",
         x = "-log10(Cutoffs) on GWAS summary stats (BetaAfr > 0)",
       y = "P-value") +
    theme_bw(24) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))


g3 = grid.arrange(g_beta, g_pv, nrow=2)
ggsave(filename = "Data/Figures/AdfricanAncestryXPGS.png", plot = g3, device = "png", width = 16, height = 20, dpi = 300)


cor(phenotype$ScaledPRS, phenotype$AfricanFraction)

vif(model_scr)

prs = as.data.frame(fread(paste0("Data/Saige.V2/PGS/WGSTractor/PRSice.best")))
this_phenotype_prs = merge(phenotype, prs %>% select(IID, PRS), by.x = "SampleID", by.y = "IID")
this_phenotype_prs = this_phenotype_prs %>% mutate(ScaledPRS = (PRS - min(PRS)) / (max(PRS) - min(PRS)))

#Run the models.
model_scr = lm(formula = MedianCreatinine ~ ScaledPRS + AfricanFraction + Age + Sex, data = this_phenotype_prs)
summary_model_scr = summary(model_scr)
print(summary(model_scr))

# Get the summary and tidy output and store effect size and beta in the results list
summary_model_scr = tidy(summary_model_scr) %>% as.data.frame()
summary_model_scr

#Implement a read function.
fread_df = function(x) {
    return(as.data.frame(fread(x))) 
}


prepare_clumping_files = function(clumping_files) { 
    # Read and store the files in a list
    clumping_files_read = lapply(clumping_files, fread_df)

    # Combine the list of data frames into one data frame
    clumping_results = do.call(rbind, clumping_files_read)

    #Reset index.
    rownames(clumping_results) = NULL
    
    
    return(clumping_results)
    
}


clumping_results = prepare_clumping_files(paste0('Data/Saige.V2/Clumping/Results/extractedChr', 
                                                    c(1:22), 
                                                    '.SaigeGWAS.Clumping.clumps'))

clumping_results %>% arrange(P) %>% head(20) %>% select(-SP2) %>% mutate(Log10P = -log10(P))

rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')
library('data.table')
library('parallel')
library('car')
library('broom')
library('gridExtra')
library('tidyr')
library('scales')

phenotype = as.data.frame(fread("Data/Saige.V2/CreatinineCohort.PhenotypeFile.tsv"))
rye = as.data.frame(fread("Data/Metadata/RyeEstimates.16From30PCS.UnrelatedParticipantsV7.Continental.7.Q"))

phenotype = merge(phenotype, rye, by = "SampleID")
phenotype = phenotype %>% mutate(AfricanFraction = African / 100)
mean(phenotype$African)
mean(phenotype$European)


dim(phenotype)
head(phenotype)
mean(phenotype$Age)

#Base model
summary(lm(formula = MedianCreatinine ~ AfricanFraction + Age + Sex, data = phenotype))

summary(lm(formula = MedianCreatinine ~ AfricanFraction + Age + Sex, data = phenotype))$coefficients[,4]  
confint(lm(formula = MedianCreatinine ~ AfricanFraction + Age + Sex, data = phenotype), level = 0.95)


# Get variants for bunch of p-values.
cutoffs = c(1:14)

results_ancestry = list()
results_prs = list()

for (cutoff_index in 1:length(cutoffs)) {
    prs = as.data.frame(fread(paste0("Data/Saige.V2/PGS/Liability/Walk/Product/PRSciseResults.Walked",
                                     cutoffs[cutoff_index],
                                     "Loci.txt")))
    this_phenotype_prs = merge(phenotype, prs %>% select(IID, PRS), by.x = "SampleID", by.y = "IID")
    this_phenotype_prs = this_phenotype_prs %>% mutate(ScaledPRS = (PRS - min(PRS)) / (max(PRS) - min(PRS)))
    
    #Run the models.
    model_scr = lm(formula = MedianCreatinine ~ ScaledPRS + AfricanFraction + Age + Sex, data = this_phenotype_prs)
    summary_model_scr = summary(model_scr)
    #print(summary(model_scr))

    # Compute 95% Confidence Intervals
    ci <- confint(model_scr, level = 0.95)

    # Get the summary and tidy output and store effect size and beta in the results list
    summary_model_scr = tidy(summary_model_scr) %>% as.data.frame()
    summary_model_scr = cbind(summary_model_scr,ci)
    row.names(summary_model_scr) = NULL
    names(summary_model_scr)[6:7] = c("Lower", "Upper")
    
    results_ancestry[[cutoff_index]] = summary_model_scr %>%
                        filter(term %in% "AfricanFraction") %>%
                        select("estimate","std.error","statistic","p.value", "Lower", "Upper") %>%
                        mutate(Cutoff = paste0(cutoffs[cutoff_index]), Predictor = "AfricanAncestry")
    
    results_prs[[cutoff_index]] = summary_model_scr %>%
                        filter(term %in% "ScaledPRS") %>%
                        select("estimate","std.error","statistic","p.value", "Lower", "Upper") %>%
                        mutate(Cutoff = paste0(cutoffs[cutoff_index]), Predictor = "PRS")
    
    rm(this_phenotype_prs)
}

results_df_ancestry = do.call(rbind, results_ancestry)
names(results_df_ancestry) = c("Estimate","StdErr", "Statistic", "Pvalue", "Lower", "Upper", "Cutoff", "Predictor")

results_df_prs = do.call(rbind, results_prs)
names(results_df_prs) = c("Estimate","StdErr","Statistic","Pvalue", "Lower", "Upper", "Cutoff", "Predictor")

base_result = data.frame(Estimate = c(0.0947168),
                                               StdErr = c(0.0077157),
                                               Statistic = c(12.3),
                                               Pvalue = c(0.000000000000000000000000000000000164927269842552),
                                               Lower = c(0.07959),
                                               Upper = c(0.10984),
                                               Cutoff = c(0),
                                               Predictor = c("AfricanAncestry"))

final_results = rbind(base_result, results_df_ancestry, results_df_prs)


final_results = final_results %>% mutate(Cutoff = factor(Cutoff, 
                                                         levels = c(0:14)))



dim(final_results)
final_results

options(repr.plot.width=20, repr.plot.height=16)

g_beta = ggplot(final_results %>% filter(Predictor == "AfricanAncestry"), 
                aes(x = Cutoff, y = Estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "darkgreen") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, color = "black") + 
    geom_point(size = 8, pch = 21, fill = "darkmagenta", color = "black", stroke = 1) +
    scale_x_discrete(breaks = seq(0,14,by=2)) +
    scale_y_continuous(limits = c(-0.0225, 0.11), 
                       breaks = c(-0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1), 
                       labels = c("-0.02", "0", "0.02", "0.04", "0.06", "0.08", "0.1")) +
    labs(x = "Number of variants (loci) in PRS", y = "African ancestry effect size (β)") +
    theme_bw(40)

g_beta

#ggsave(filename = "Data/Figures/Supplment/AfricanAncestryXPGS.LiabilityWalks.DeltaProduct.Beta.pdf", 
#       plot = g_beta, device = "pdf", 
#       width = 15, height = 12, dpi = 300)

g_pv = ggplot(final_results %>% filter(Predictor == "AfricanAncestry"), 
              aes(x = Cutoff, y = -log10(Pvalue), group = T)) +
    geom_line(color = "#3f003f", linewidth = 1.5) +
    geom_point(size = 8, pch = 21, fill = "darkmagenta", color = "black", stroke = 1) +
    scale_x_discrete(breaks = seq(0,14,by=2)) +
    labs(x = "Number of variants (loci) in PS", y = "African ancestry") +
    theme_bw(base_size = 40) 

g_pv

#ggsave(filename = "Data/Figures/Supplment/AfricanAncestryXPGS.LiabilityWalks.DeltaProduct.Significance.pdf", 
#       plot = g_pv, device = "pdf", 
#       width = 15, height = 12, dpi = 300)







g_beta = ggplot(final_results, aes(x = Cutoff, y = Estimate)) +
    geom_errorbar(aes(ymin = Estimate - StdErr, ymax = Estimate + StdErr),
                width = 0.2, color = "black") + 
    geom_point(size = 4, color = "darkmagenta") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "darkgreen") +
    scale_x_discrete(breaks = seq(0,70,by=5)) +
    scale_y_continuous(limits = c(-0.1, 0.15), 
                       breaks = c(-0.1, -0.05, 0, 0.05, 0.1, 0.15), 
                       labels = c("-0.1", "-0.05", "0", "0.05", "0.1", "0.15")) +
    facet_wrap(~Predictor, scales = "free_y") +
    labs(title = "Comparison of estimates",
       x = "Clumps considered based on 2xAFxBETA",
       y = "Estimate") +
    theme_bw(40) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

g_pv = ggplot(final_results, aes(x = Cutoff, y = -log10(Pvalue))) +
    geom_point(size = 6, color = "darkmagenta") +  
    facet_wrap(~Predictor, scales = "free_y") +
    scale_x_discrete(breaks = seq(0,70,by=5)) +
    labs(title = "Comparison of p-values",
         x = "Clumps considered based on 2xAFxBETA",
       y = "-log10(P-value)") +
    theme_bw(40) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))


#g3 = grid.arrange(g_beta, g_pv, nrow=2)
#g3
#ggsave(filename = "Data/Figures/AdfricanAncestryXPGS.LiabilityWalks.png", plot = g3, device = "png", 
#       width = 36, height = 20, dpi = 300)

g_beta
ggsave(filename = "Data/Figures/Supplment/AfricanAncestryXPGS.LiabilityWalks.Delta.png", plot = g_beta, device = "png", 
       width = 18, height = 10, dpi = 300)

















fam = as.data.frame(fread("Data/Sarah/ACAF_Subset/extractedChr6.ACAF.Plink.fam"))
status = as.data.frame(fread("Data/Phecodes/Cohorts/cohort.250.2.tsv"))
status = status %>% filter(PersonID %in% fam$V1 & Status %in% c("Case", "Control"))
status = status %>% select(PersonID, Status)
status = status %>% mutate(Status = ifelse(Status == "Case", 1, 0))
status %>% count(Status)

names(status) = c("IID", "Status")

write.table(x = status, file = "Data/Sarah/PhenotypeFile.tsv", sep ="\t", row.names = F, quote = F)
write.table(x = status %>% select(IID) %>% mutate(FID = IID), 
            file = "Data/Sarah/CohortSelected.txt", sep ="\t", row.names = F, quote = F)


head(fam)
head(status)

stats = as.data.frame(fread("Data/Sarah/t2dSNP_fromAoU.csv"))
stats = stats %>% mutate(MarkerID = paste("chr",id, sep = ""), P = 0)
names(stats) = c("SNP", "EA", "BETA", "CHR", "BP", "A1", "A2", "MarkerID", "P")
write.table(x = stats %>% select(-SNP), file = "Data/Sarah/KPtransformed.tsv", sep ="\t", row.names = F, quote = F)
head(stats)


bim = as.data.frame(fread("Data/Sarah/ACAF_Subset/extractedChr10.ACAF.Plink.bim"))
head(bim)

bim %>% filter(V4 %in% stats$BP)
stats %>% filter(BP %in% bim$V4)

#Run PRScise for this walk.
prscise_command = paste("Rscript ~/workspaces/admixturemapping/bin/PRSice/PRSice.R",
                    "--prsice ~/workspaces/admixturemapping/bin/PRSice/PRSice_linux",
                    "--base Data/Sarah/KPtransformed.tsv",
                    "--A1 A2",
                    "--A2 A1",
                    "--chr CHR",
                    "--bp BP",
                    "--stat BETA",
                    "--pvalue P",
                    "--snp MarkerID",
                    "--thread 12",
                    "--memory 97000",
                    "--binary-target T",
                    "--target Data/Sarah/ACAF_Subset/extractedChr6.ACAF.Plink",
                    "--pheno Data/Sarah/PhenotypeFile.tsv",
                    "--pheno-col Status",
                    "--ignore-fid",
                    "--fastscore",
                    "--keep-ambig",
                    sep = " ")

prs_command_output = system(prscise_command, intern = TRUE)


