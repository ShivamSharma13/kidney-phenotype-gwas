rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')
library('ggthemes')
library('data.table')
library("GWASTools")
library('GENESIS')
library('stringr')
library('gdsfmt')
library('qqman')

options(repr.plot.width = 24, repr.plot.height = 12, repr.plot.res = 100)
options(repr.matrix.max.cols=50, repr.matrix.max.rows=500)

creatinine_covars = as.data.frame(fread("Data/Metadata/CreatinineCohort.V2.tsv"))
dim(creatinine_covars)
head(creatinine_covars)

#Look at the distribution of the original data.
g = ggplot(creatinine_covars, aes(x=MedianCreatinine))+
    geom_density(color="black", fill="lightblue") +
    theme_bw(base_size = 24)

ggsave(filename = "Data/Saige/Figures/MedianCreatinine.png", g, device = "png", width = 12, height = 8, dpi = 300)

###--------------INT------------------###
#Perform INT on Median Creatinine levels.
creatinine_covars$MedianCreatinineINT = qnorm((rank(creatinine_covars$MedianCreatinine,na.last="keep")-0.5)/
                                              sum(!is.na(creatinine_covars$MedianCreatinine)))
head(creatinine_covars)

#Look at the distribution of INT.
ggplot(creatinine_covars, aes(x=MedianCreatinineINT))+
    geom_density(color="black", fill="lightblue") +
    theme_bw(base_size = 24)

#Compare INT vs raw values for creatinine.
ggplot(creatinine_covars, aes(x = MedianCreatinineINT, y = MedianCreatinine)) +
    geom_point(pch = 21, fill = "blue4", color = "black", size = 3) +
    theme_classic(base_size = 24)

##Add in the PCAs.
pca = as.data.frame(fread("Data/Genetic/AllofUs/PCA/TwoWay.PCA.eigenvec"))
pca = pca %>% select(IID, PC1, PC2, PC3, PC4, PC5)
creatinine_covars = merge(creatinine_covars, pca, by.x = "SampleID", by.y = "IID")

dim(creatinine_covars)
head(creatinine_covars)

#Get all the msp files together and vertically stack them.
msp_files = paste0('Data/Genetic//LocalAncestry//TwoWay/Gnomix/chr', c(1:22), '/query_results.msp')
msp_files

## Stack all the MSP files.

#Implement a read function.
fread_df = function(x) {
    return(as.data.frame(fread(x))) 
}

# Read and store the files in a list
msp_files_read <- lapply(msp_files, fread_df)

# Combine the list of data frames into one data frame
msp <- do.call(rbind, msp_files_read)

#Reset index.
rownames(msp) <- NULL

dim(msp)
head(msp)

#Get the indices on which dosages exist for all the samples. It starts from 7 and ends at end of df.
dosage_index_in_msp = c(7:length(names(msp)))

#Get a LociID and ancestral block length.
msp$LociID = floor((msp$spos + msp$epos) / 2)
msp$Chr = msp$`#chm`

head(msp)

#Get the columns representing individuals. Values a re hardcoded but can be automated later.
individuals <- colnames(msp)[dosage_index_in_msp]
individual_first_copies = paste0(unique(gsub("\\..*", "", individuals)), ".0")
individual_dosages = paste0(individual_first_copies, "_Dosage")

head(individuals)
head(individual_first_copies)
tail(individual_first_copies)

length(individuals)
length(individual_first_copies)

#Add the haploid columns together to get the diploid dosages.
msp_dosages <- msp %>% 
    mutate(across(all_of(individual_first_copies), ~ .+get(str_c(gsub("\\..*", "",cur_column()),".1" )), 
                  .names = '{.col}_Dosage' )) %>%
    select(all_of(c("LociID", "Chr", individual_dosages)))
    

head(msp_dosages)
dim(msp_dosages)

#Fix colnames.
colnames(msp_dosages) = gsub("\\..*", "",colnames(msp_dosages))

dim(msp_dosages)
head(msp_dosages)


aou_individuals = colnames(msp_dosages)[grepl("^AOU", colnames(msp_dosages))]

#Filter if acceptable creatinine levels are present for these people.
aou_individuals = aou_individuals[aou_individuals %in% paste0("AOU_", creatinine_covars$SampleID)]
head(aou_individuals)
length(aou_individuals)

#Extract AOU particiapnts.
msp_summarized = msp_dosages %>% select(all_of(aou_individuals))
colnames(msp_summarized) = gsub("AOU_", "",colnames(msp_summarized))
msp_summarized = as.matrix(msp_summarized)

dim(msp_summarized)
head(msp_summarized)

#Create attributes which can be used for GDS later.
snpID = as.integer(gsub("^\\d+:", "", msp_dosages$LociID))
chromosome = as.integer(gsub(":\\d+", "", msp_dosages$Chr))
positions = as.integer(gsub("^\\d+:", "", msp_dosages$LociID))
scanID = as.integer(colnames(msp_summarized %>% as.data.frame()))

length(snpID)
length(chromosome)
length(positions)
length(scanID)

head(snpID)
head(chromosome)
head(positions)
head(scanID)

##Order covar df on scanID. Genesis needs dosage and "attributes" separately, so the order must match.
order_indices = match(scanID, creatinine_covars$SampleID)
creatinine_covars = creatinine_covars[order_indices, ]
rownames(creatinine_covars) <- NULL

dim(creatinine_covars)
head(creatinine_covars)

#Show the files.
showfile.gds(closeall=FALSE, verbose=TRUE)

#Show and close the files.
showfile.gds(closeall=TRUE, verbose=TRUE)

#Delete the GDS files.
unlink("Data/Genetic/LocalAncestry/TwoWay/Gnomix/GDS/Dosage.AFR.Creatinine.chrAll.gds", force=TRUE)

#Create a GDS file
afr_gds <- createfn.gds("Data/Genetic/LocalAncestry/TwoWay/Gnomix/GDS/Dosage.AFR.Creatinine.chrAll.gds")

#Add stuff for AFR GDS.
add.gdsn(afr_gds, "genotype", msp_summarized)
add.gdsn(afr_gds, "snp.id", snpID)
add.gdsn(afr_gds, "snp.chromosome", chromosome)
add.gdsn(afr_gds, "snp.position", positions)
add.gdsn(afr_gds, "sample.id", scanID)

#Show and close the files.
showfile.gds(closeall=TRUE, verbose=TRUE)

gdsList = lapply(list("Data/Genetic/LocalAncestry/TwoWay/Gnomix/GDS/Dosage.AFR.Creatinine.chrAll.gds"), 
                 GdsGenotypeReader)

scanAnnot <- ScanAnnotationDataFrame(data.frame(
    scanID=getScanID(gdsList[[1]]), stringsAsFactors=FALSE))

#Set phenotypes.
scanAnnot$MedianCreatinineINT = creatinine_covars$MedianCreatinineINT
scanAnnot$Age = creatinine_covars$Age
scanAnnot$Sex = creatinine_covars$Sex

#Add PCs.
scanAnnot$PC1 = creatinine_covars$PC1
scanAnnot$PC2 = creatinine_covars$PC2
scanAnnot$PC3 = creatinine_covars$PC3
scanAnnot$PC4 = creatinine_covars$PC4
scanAnnot$PC5 = creatinine_covars$PC5


#Add Genetic and Phenotypic information together.
genoDataList = lapply(gdsList, GenotypeData, scanAnnot=scanAnnot)

#Read the GRM header.
grm_matrix_header = as.data.frame(fread("Data/Genetic/AllofUs/PCA/TwoWay.PCA.rel.id")) %>% pull(IID)
head(grm_matrix_header)

#Read the GRM matrix.
grm_matrix = as.matrix(fread("Data/Genetic/AllofUs/PCA/TwoWay.PCA.rel"))

head(grm_matrix)

colnames(grm_matrix) <- grm_matrix_header
rownames(grm_matrix) <- grm_matrix_header

head(grm_matrix)

#Get individuals in the right order.
grm_matrix_aou <- grm_matrix[as.character(scanID), as.character(scanID)]
head(grm_matrix_aou)

null.model <- fitNullModel(scanAnnot, outcome="MedianCreatinineINT", 
                           covars=c("Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5"),
                          cov.mat = grm_matrix_aou)

genoIterators <- lapply(genoDataList, GenotypeBlockIterator)

admix_map_associations <- admixMap(genoIterators, null.model)

#Write the associations.
write.table(admix_map_associations, file = "Data/Genetic/AdmixtureMapping/TwoWay/Results/Creatinine.Results.tsv", row.names = F, quote = F, sep = "\t")


### In case we already ran the association tests.
#admix_map_associations = as.data.frame(fread("Data/Genetic/AdmixtureMapping/TwoWay/Results/Creatinine.Results.tsv"))
admix_map_associations = as.data.frame(fread("Data/AdmixtureMapping/Creatinine.Results.V1.tsv"))
head(admix_map_associations)

admix_map_associations %>% arrange(pval) %>% head()

ggplot_x = admix_map_associations %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=as.numeric(max(pos))) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(admix_map_associations, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate(XPosition=pos+tot)

  
ggplot_x_labels = ggplot_x %>% group_by(chr) %>% summarize(center=( max(XPosition) + min(XPosition) ) / 2 )

head(ggplot_x)
head(ggplot_x_labels)

g = ggplot(ggplot_x, aes(x=XPosition, y=-log10(pval))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("blue3", "darkorange1"), 22)) +
    
    labs(x = "Position", y = "-log10(p-value)") +

    # custom X axis:
    scale_x_continuous( label = ggplot_x_labels$chr, breaks= ggplot_x_labels$center ) +
    scale_y_continuous(limits = c(0, 5.5) ) +     # remove space between plot area and x axis
    
    #Horizontal line.
    geom_hline(yintercept=4.42)+

    # Custom the theme:
    theme_few(base_size = 32) +
    theme(legend.position="none")

ggsave(filename = "Data/Figures/AdmixtureMapping.Creatinine.png", g, device = "png", 
       width = 18, height = 14, dpi = 300)

g




