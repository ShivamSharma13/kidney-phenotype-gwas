rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')
library('data.table')
library('parallel')
library('broom') 
library('gridExtra') 
library('ggVennDiagram')

options(repr.plot.width=16, repr.plot.height=10)

creatinine_covars = as.data.frame(fread("Data/Saige.European/CreatinineCohort.European.tsv"))
dim(creatinine_covars)
head(creatinine_covars)

##Add in the PCAs.
pca = as.data.frame(fread("Data/Saige.European/PCA/extractedChrAllPruned.PCA.eigenvec"))
pca = pca %>% select(IID, paste0("PC", c(1:10)))
creatinine_covars = merge(creatinine_covars, pca, by.x = "SampleID", by.y = "IID")
creatinine_covars = creatinine_covars %>% mutate(AgeSq = Age*Age)

creatinine_covars$SexCoded = ifelse(creatinine_covars$Sex == "Female", 0, 
                   ifelse(creatinine_covars$Sex == "Male", 1, -9))

dim(creatinine_covars)
head(creatinine_covars)

write.table(x = creatinine_covars, 'Data/Saige.European/CreatinineCohort.PhenotypeFile.tsv', 
            sep = '\t', quote = F, row.names = F)

creatinine_covars = as.data.frame(fread("Data/Saige.V2/CreatinineCohort.tsv"))
dim(creatinine_covars)
head(creatinine_covars)

#Look at the distribution of the original data.
g = ggplot(creatinine_covars, aes(x=MedianCreatinine))+
    facet_wrap(~Sex) +
    geom_histogram(position="identity", color = "black", fill = "darkmagenta", bins = 25)+
    labs(x = "Serum Creatinine (mg/dL)", y = "Number of observations")   +
    theme_bw(base_size = 36)

ggsave("Data/Figures/Figure1/CreatinineLevels.png", g, 
       width = 14, height = 8, dpi = 300, units = "in", device='png')
g

###--------------INT------------------###
#Perform INT on Median Creatinine levels.
creatinine_covars$MedianCreatinineINT = qnorm((rank(creatinine_covars$MedianCreatinine,na.last="keep")-0.5)/
                                              sum(!is.na(creatinine_covars$MedianCreatinine)))
head(creatinine_covars)

#Look at the distribution of INT.
g = ggplot(creatinine_covars, aes(x=MedianCreatinineINT))+
    geom_histogram(position="identity", color = "black", fill = "darkmagenta", bins = 25)+
    labs(x = "Inverse Normal Serum Creatinine (mg/dL)", y = "Number of observations")   +
    facet_wrap(~Sex) +
    theme_bw(base_size = 24)

ggsave("Data/Figures/Figure1/CreatinineLevels.INT.png", g, 
       width = 14, height = 8, dpi = 300, units = "in", device='png')

##Add in the PCAs.
pca = as.data.frame(fread("Data/Saige.V2/PCA/extractedChrAllPruned.CreatinineCohort.PCA.eigenvec"))
pca = pca %>% select(IID, paste0("PC", c(1:10)))
creatinine_covars = merge(creatinine_covars, pca, by.x = "SampleID", by.y = "IID")
creatinine_covars = creatinine_covars %>% mutate(AgeSq = Age*Age)

creatinine_covars$SexCoded = ifelse(creatinine_covars$Sex == "Female", 0, 
                   ifelse(creatinine_covars$Sex == "Male", 1, -9))

dim(creatinine_covars)
head(creatinine_covars)

write.table(x = creatinine_covars, 'Data/Saige.V2/CreatinineCohort.PhenotypeFile.tsv', 
            sep = '\t', quote = F, row.names = F)

#write.table(x = creatinine_covars %>% mutate(SampleID2 = SampleID) %>% select(SampleID,SampleID2), 
#            'Data/Saige/Cohort.Creatinine.ForPlink.tsv', 
#            sep = '\t', quote = F, row.names = F, col.names = F)


creatnine_pca = as.data.frame(fread("Data/PCA/extractedChrAllPruned.PCA.eigenvec"))
names(creatnine_pca)[1] = "Pop"

dim(creatnine_pca)
head(creatnine_pca)
unique(creatnine_pca$Pop)

creatnine_pca = within(creatnine_pca, Population <- 
            ifelse(Pop %in% c("forReferenceGBR","forReferenceIBS", "forReferenceTSI"), "European",
            ifelse(Pop %in% c("forReferenceESN","forReferenceGWD", "forReferenceYRI"), "African",
            ifelse(Pop %in% c("AOU"), "AllofUs", "None"))))

#Check the DF.
dim(creatnine_pca)
head(creatnine_pca)

colors = c("African" = "#283593",
          "European" = "#FFA000", 
          "AllofUs" = "#e8e8e8")

g = ggplot(creatnine_pca %>% filter(Population == "AllofUs"),
       aes(x = PC1*-1, y = PC2, fill = as.factor(Population))) +
    geom_point(stroke = 0.25, pch = 21, size = 2, color = "black") +  
    geom_point(data = creatnine_pca %>% filter(Population != "AllofUs"), 
               aes(x = PC1*-1, y = PC2), stroke = 0.25, pch = 21, size = 2) +  
    labs(x = "PC1", y = "PC2") +
    scale_fill_manual(values = colors, guide = "none") +
    theme_bw(base_size=40) +
    theme(legend.position = "none")

g
ggsave("Data/Figures/Figure1/PCAPlot.png", g, 
       width = 10, height = 10, dpi = 300, units = "in", device='png')


rye = as.data.frame(fread("Data/Metadata/RyeEstimates.16From30PCS.UnrelatedParticipantsV7.Continental.7.Q"))
rye = rye %>% filter(SampleID %in% (creatnine_pca %>% filter(Pop == "AOU") %>% pull(IID)))

dim(rye)
head(rye)

count(rye, SelfReportedRaceEthnicity)
16950/18979*100

#Define colors for admixture plots.
continental_ancestry_colors = c("SouthAsian" = "#B71C1C", "EastAsian" = "seagreen4", "WestAsian" = "#99887a", 
                               "European" = "#FFA000", "African" = "#283593",
                               "NativeAmerican" = "#41C9F8", "Oceania" = "#E040FB")

#Sort the dataset by max ancestry column.
dataset_sorted = rye[order(rye[,"African"], decreasing = TRUE),]

#this_df = dataset %>% select(African, European, SouthAsian, EastAsian,
#                           NativeAmerican, Oceania, WestAsian)
#this_order = hclust(dist(this_df))$order
#dataset_sorted = dataset[this_order,]

row.names(dataset_sorted) <- NULL

#Melt the dataset and preserve the index for the admixture plot.
dataset_sorted$index= as.numeric(rownames(dataset_sorted))

dataset_sorted_melt = reshape2::melt(data = dataset_sorted, 
                                    id.vars = c('SampleID', 'index'), 
                                    measure.vars = c('European', 'African', 'NativeAmerican',
                                               'SouthAsian', 'EastAsian', 'Oceania', 'WestAsian'))
colnames(dataset_sorted_melt)[3] <- 'Ancestry'
dataset_sorted_melt = dataset_sorted_melt %>% arrange(index)

#Make the ggplot.
rye_continental_admixture_plot = ggplot(data = dataset_sorted_melt,
                        aes(x=index, y=value, fill=Ancestry)) +
                geom_bar(stat="summary", width=1.001) + 
                scale_fill_manual(values=continental_ancestry_colors) + 
                theme_classic()+
                theme(line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
                      axis.text=element_blank(), axis.ticks=element_blank(),  
                     legend.position = "none") 

#Save the plot.
ggsave(filename = "Data/Figures/Figure1/AdmixturePlot.png", width = 16, height = 10, dpi = 600, 
       units = "in", device='png', limitsize = FALSE)

rye_continental_admixture_plot

creatinine_phenotype = as.data.frame(fread("Data/Saige.V2/CreatinineCohort.PhenotypeFile.tsv")) %>% select(SampleID, MedianCreatinine, MedianCreatinineINT, Age, Sex)
creatinine_phenotype = merge(creatinine_phenotype, rye, by = "SampleID")

mean(creatinine_phenotype$Age)
count(creatinine_phenotype, Sex)

dim(creatinine_phenotype)
head(creatinine_phenotype)

# Define the bins
breaks <- seq(0, 100, by = 2) 
labels <- paste((head(breaks, -1) + tail(breaks, -1)) / 2 , sep = "-")  

# Create summary data frame with mean of C2 for each bin
admixreg_scr_binned = creatinine_phenotype %>%
    mutate(AncestryBin = as.numeric(as.character(cut(African, breaks = breaks, labels = labels, include.lowest = TRUE)))) %>%
    group_by(AncestryBin, Sex) %>%
    summarize(MeanCreatinineBin = mean(MedianCreatinine, na.rm = TRUE),
              SECreatinineBin = sd(MedianCreatinine, na.rm = TRUE) / sqrt(n()),
              .groups = 'drop')

head(admixreg_scr_binned)
tail(admixreg_scr_binned)

correlation_results <- admixreg_scr_binned %>%
  group_by(Sex) %>%
  summarize(Correlation = cor(MeanCreatinineBin, as.numeric(AncestryBin), 
                              use = "everything", 
                              method = "spearman"),
           CorrelationSq = Correlation^2)

correlation_results

text_data <- data.frame(
  Sex = c("Female", "Male"),
  Label = c("Spearman \u03C1 = 0.79", "Spearman \u03C1 = 0.84")
)

g = ggplot(admixreg_scr_binned, aes(x = AncestryBin, y = MeanCreatinineBin)) +
    geom_smooth(method='lm', formula= y~x, linewidth = 2, color = "dodgerblue4")+
    geom_errorbar(aes(ymin = MeanCreatinineBin - SECreatinineBin, ymax = MeanCreatinineBin + SECreatinineBin), 
                  width = 0.2, color = "dodgerblue4", linewidth = 1) +  # Error bars for SE
    geom_point(pch = 21, fill = "dodgerblue", color = "black", size = 8) +
    scale_x_continuous(limits = c(0,100)) +
    labs(y = "Creatinine Levels (mg/dL)", x = "Percentage African Ancestry") +
    facet_wrap(~Sex) +
    geom_text(data = text_data, aes(x = -Inf, y = Inf, label = Label), 
               vjust = 1.5, hjust = -0.1, size = 10, color = "black") +
    theme_bw(base_size = 40)

ggsave(filename = "Data/Figures/Figure1/AdmixtureRegression.pdf", g, device = "pdf", 
       width = 16, height = 10, dpi = 300, units = "in")

g

#Get dosage for lead variant.
dosage = as.data.frame(fread("Data/Saige.V2/Dosage/extractedChr15.ACAF.GWASThreshold.Dosage.raw")) %>% select(IID, `chr15:45379909:C:T_C`)
names(dosage) = c("SampleID", "LeadVariant1")
dosage = dosage %>% mutate(LeadVariant1 = 2 - LeadVariant1)

#Merge with phenotype data.
creatinine_phenotype = merge(creatinine_phenotype, dosage, by = "SampleID")
head(creatinine_phenotype)

# Define the bins
breaks <- seq(0, 100, by = 2) 
labels <- paste((head(breaks, -1) + tail(breaks, -1)) / 2 , sep = "-")  

# Create summary data frame with mean of C2 for each bin
admixreg_scr_binned_condition = creatinine_phenotype %>%
    filter(!is.na(LeadVariant1)) %>%
    mutate(AncestryBin = as.numeric(as.character(cut(African, breaks = breaks, labels = labels, include.lowest = TRUE)))) %>%
    group_by(LeadVariant1, AncestryBin, Sex) %>%
    summarize(MeanCreatinineBin = mean(MedianCreatinine, na.rm = TRUE),
              SECreatinineBin = sd(MedianCreatinine, na.rm = TRUE) / sqrt(n()),
              .groups = 'drop')

head(admixreg_scr_binned_condition)
tail(admixreg_scr_binned_condition)

correlation_results <- admixreg_scr_binned_condition %>%
  group_by(LeadVariant1, Sex) %>%
  summarize(Correlation = cor(MeanCreatinineBin, as.numeric(AncestryBin), 
                              use = "everything", 
                              method = "spearman"),
           CorrelationSq = Correlation^2)

correlation_results

options(repr.plot.width=16, repr.plot.height=24)

text_data <- data.frame(
    LeadVariant1 = c("0", "0", "1", "1", "2", "2"),
  Sex = c("Female", "Male", "Female", "Male", "Female", "Male"),
  Label = c("Spearman \u03C1 = 0.66", "Spearman \u03C1 = 0.63",
           "Spearman \u03C1 = 0.66", "Spearman \u03C1 = 0.15",
           "Spearman \u03C1 = 0.35", "Spearman \u03C1 = 0.29")
)

g = ggplot(admixreg_scr_binned_condition, aes(x = AncestryBin, y = MeanCreatinineBin)) +
    geom_smooth(method='lm', formula= y~x, linewidth = 2, color = "dodgerblue4")+
    geom_errorbar(aes(ymin = MeanCreatinineBin - SECreatinineBin, ymax = MeanCreatinineBin + SECreatinineBin), 
                  width = 0.2, color = "dodgerblue4", linewidth = 1) +  # Error bars for SE
    geom_point(pch = 21, fill = "dodgerblue", color = "black", size = 8) +
    scale_x_continuous(limits = c(0,100)) +
    labs(y = "Creatinine Levels (mg/dL)", x = "Percentage African Ancestry") +
    facet_grid(LeadVariant1~Sex) +
    geom_text(data = text_data, aes(x = -Inf, y = Inf, label = Label), 
               vjust = 1.5, hjust = -0.1, size = 10, color = "black") +
    theme_bw(base_size = 40)

g

options(repr.plot.width=16, repr.plot.height=10)





creatinine_phenotype$AncestryBin <- cut(creatinine_phenotype$African, breaks = seq(0, 100, by = 10), 
                 #labels = c("0-20", "20-40", "40-60", "60-80", "80-100"), 
                 #labels = c("0-25", "25-50", "50-75", "75-100"), 
                 #labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"), 
                labels = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),  
                include.lowest = TRUE, right = FALSE)

dim(creatinine_phenotype)
head(creatinine_phenotype)

g = ggplot(creatinine_phenotype, aes(x = AncestryBin, y = MedianCreatinine)) +
    geom_boxplot(fill = "dodgerblue", outliers = F) +
    geom_smooth(aes(group = 1), method = "lm", color = "dodgerblue4", size = 2) +  # Smooth line
    labs(y = "Creatinine Levels (mg/dL)", x = "Percentage African Ancestry") +
    scale_y_continuous(limits = c(0.38,1.45), breaks = c(0.4,0.6,0.8,1.0,1.2,1.4))+
    facet_wrap(~Sex) +
    theme_bw(base_size = 36)
    #theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

ggsave(filename = "Data/Figures/Figure1/AdmixtureRegression.Boxplot.png", g, device = "png", width = 18, height = 10)
g







head(creatinine_phenotype)

# Fit a linear model
creatinine_phenotype = creatinine_phenotype %>% mutate(AfricanFraction = African / 100)
scr_model <- lm(MedianCreatinine ~ AfricanFraction + Age + Sex, data = creatinine_phenotype)

# Tidy the model output
scr_model_summary <- tidy(scr_model) %>%
                mutate(
                    term = gsub("Predictor", "Effect Size", term),  # Rename predictors
                    lower = estimate - 1.96 * std.error,
                    upper = estimate + 1.96 * std.error) %>%
                    filter(term != "(Intercept)") %>%
                    as.data.frame()

scr_model_summary$term = c("African\nancestry", "Age", "Male")
scr_model_summary$ScientificPalue = c("1.65 x 10\u207B\u00B3\u2074", "7.38 x 10\u207B\u00B9\u2079\u2078", "~0")
scr_model_summary

# Create the forest plot
f1 = ggplot(scr_model_summary, aes(y = term, x = estimate)) +
                geom_vline(xintercept = 0, linetype = "longdash", color = "darkred", linewidth = 1) +
                geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, linewidth = 1, color = "dodgerblue4") +  # Confidence intervals
                geom_point(size = 8, color = "dodgerblue") +  # Effect sizes
                labs(x = "Effect Size (\u3b2)") +
                theme_bw(base_size = 40) +
                theme(axis.title.y = element_blank(), plot.margin=unit(c(1,0,1,1), "cm"))


f2 = ggplot(scr_model_summary, aes(x = factor(1), y = term)) +  # Use a single x value for alignment
  geom_text(aes(label = paste("p =", ScientificPalue)), size = 7, vjust = 0.5) +  # Text with p-values
  labs(x = "P-value") +
  theme_bw(base_size = 40) +
  theme(axis.title.y = element_blank(),  # Remove y-axis title
        axis.ticks = element_blank(),     # Remove ticks
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
       plot.margin=unit(c(1,1,2.375,-0.35), "cm")) 

g = grid.arrange(f1, f2, nrow = 1, widths = c(0.75, 0.25), padding = unit(0, "cm"))

ggsave(filename = "Data/Figures/Figure1/ForestPlot.pdf", g, device = "pdf", width = 10, height = 10)




t2d = as.data.frame(fread("Data/Phecodes/Cohorts/cohort.250.2.tsv"))
head(t2d)
count(t2d, Status)

rye = as.data.frame(fread("Data/Metadata/RyeEstimates.16From30PCS.UnrelatedParticipantsV7.Continental.7.Q"))
dataset = merge(t2d, rye, by.x = "PersonID", by.y = "SampleID")
dataset = dataset %>% filter(African > 5 & (SouthAsian+EastAsian+NativeAmerican+Oceania+WestAsian<5))
dim(dataset)
head(dataset)

# Define the bins
breaks <- seq(0, 100, by = 2) 
labels <- paste((head(breaks, -1) + tail(breaks, -1)) / 2 , sep = "-")  

# Create summary data frame with mean of C2 for each bin
admixreg_t2d_binned = dataset %>%
    mutate(AncestryBin = as.numeric(as.character(cut(African, breaks = breaks, labels = labels, include.lowest = TRUE)))) %>%
    group_by(AncestryBin) %>%
    summarise(
        Cases = sum(Status == "Case"),
        Controls = sum(Status == "Control"),
        Prevalence = Cases / (Cases + Controls),
        SE = sqrt((Prevalence * (1 - Prevalence)) / (Cases + Controls))
    )

head(admixreg_t2d_binned)
tail(admixreg_t2d_binned)

g = ggplot(admixreg_t2d_binned, aes(x = AncestryBin, y = Prevalence)) +
    geom_errorbar(aes(ymin = Prevalence - SE, ymax = Prevalence + SE), 
                  width = 0.2, color = "dodgerblue4", linewidth = 1) +  # Error bars for SE
    geom_point(pch = 21, fill = "dodgerblue", color = "black", size = 8) +
    scale_x_continuous(limits = c(40,100)) +
    ylim(c(0.1,0.4))+
    labs(y = "Type 2 diabetes prevalence", x = "Percentage African Ancestry") +
    theme_bw(base_size = 36)

ggsave(filename = "Data/Figures/T2DRisk.V2.png", g, device = "png", width = 18, height = 12)


g

models = as.data.frame(fread("Data/tmp/AllSummary.tsv"))
dim(models)
head(models)

models %>% group_by(Ancestry) %>% summarise(Counts = sum(AncestryPvalue < 0.00002986857))

filtered_models = models %>% filter(AncestryPvalue < 0.00002986857) 
filtered_model_sets <- list(African = filtered_models %>% filter(Ancestry == "African") %>% pull(Phecode),
                        American = filtered_models %>% filter(Ancestry == "NativeAmerican") %>% pull(Phecode),
                        European = filtered_models %>% filter(Ancestry == "European") %>% pull(Phecode))


# Create a Venn diagram using ggVennDiagram
g = ggVennDiagram(filtered_model_sets, 
              category.names = c("African", "American", "European"),
              show.percent = TRUE,label_size = 12,label_percent_digit=2,set_size = 12) +
ggplot2::scale_fill_gradient(low="white",high = "darkred") +
theme_void(base_size = 24)+
theme(legend.position = "bottom")
ggsave(filename = "Data/Figures/Overlap.pdf", g, device = "pdf", width = 16, height = 16, dpi = 300)


?ggVennDiagram


