rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('ggplot2')
library('data.table')
library('parallel')
library('stringr')
library("ggpubr")
library("ggthemes")
library('ggrepel')

options(repr.plot.width=16, repr.plot.height=10)

mane = as.data.frame(fread("Data/Metadata/MANE.GRCh38.v1.2.refseq_genomic.NoComments.gff", header = F))
names(mane) = c("CHR", "Source", "Element", "Start", "Stop", "V6", "V7", "V8", "Description")
mane = mane %>% filter(Element == "gene") %>% select(CHR, Start, Stop, Description)

#Get the gene name please.
mane = mane %>% 
        mutate(Gene = str_extract(Description, "(?<=Name=)[^;]+")) %>% 
        mutate(CHR = as.integer(sub("chr", "", CHR))) %>%
        select(-Description)

dim(mane)
head(mane)

find_gene <- function(row, mane) {
    chr <- row['CHR']$CHR
    pos <- row['POS']$POS
    gene_info = mane %>% dplyr::filter(CHR == chr & Start < pos & Stop > pos)
      if (nrow(gene_info) == 0) {
            return(NA)
      } else {
            return(gene_info[1,]$Gene)
      }
}

#Implement a read function.
fread_df = function(x) {
    return(as.data.frame(fread(x))) 
}


prepare_exwas_files = function(exwas_result_files) { 
    # Read and store the files in a list
    exwas_files_read = lapply(exwas_result_files, fread_df)

    # Combine the list of data frames into one data frame
    exwas_results = do.call(rbind, exwas_files_read)

    #Reset index.
    rownames(exwas_results) = NULL
    
    ## Define the rows which need to be annotated.
    exwas_anno = exwas_results %>% filter(p.value <= 0.0000000501187234)
    exwas_noanno = exwas_results %>% filter(p.value > 0.0000000501187234) %>% 
                mutate(Gene = "None") 
    
    variant_lists <- split(exwas_anno, seq(nrow(exwas_anno)))
    GeneNames <- mclapply(variant_lists, find_gene, mane = mane, mc.cores = 1)

    # Combine the results back into DF1
    exwas_anno$Gene <- unlist(GeneNames)
    
    
    #Lable only the leading variant.
    exwas_anno = exwas_anno %>%    
                            group_by(Gene) %>%
                            mutate(Gene = ifelse(p.value == min(p.value), Gene, "None")) %>%
                            ungroup()
    
    exwas_results = rbind(exwas_anno, exwas_noanno)
    
    exwas_results = exwas_results %>% arrange(p.value)
    exwas_results_x = exwas_results %>% 
                    # Compute chromosome size
                    group_by(CHR) %>% 
                    summarise(Chr_len=as.numeric(max(POS))) %>%
                    mutate(Tot=cumsum(Chr_len)-Chr_len) %>%
                    select(-Chr_len) %>%

                    # Add this info to the initial dataset
                    left_join(exwas_results, ., by=c("CHR"="CHR")) %>%

                    # Add a cumulative position of each SNP
                    arrange(CHR, POS) %>%
                    mutate(XPosition=POS+Tot)

    
    exwas_results_x_ggplot = exwas_results_x %>% 
                        group_by(CHR) %>% 
                        summarize(center=( max(XPosition) + min(XPosition) ) / 2 ) %>% as.data.frame()
    
    return(list(exwas_results = exwas_results, 
                exwas_results_x = exwas_results_x, 
                exwas_results_x_ggplot = exwas_results_x_ggplot))
    
        
}

scr_gwas_results = prepare_exwas_files(paste0('Data/Saige.V2/Step2/extractedChr', 
                                                    c(1:22), 
                                                    '.SaigeGWAS'))

scr_gwas = scr_gwas_results$exwas_results
scr_gwas_x = scr_gwas_results$exwas_results_x
scr_gwas_ggplot = scr_gwas_results$exwas_results_x_ggplot

head(scr_gwas_x)

scr_gwas_x = scr_gwas_x %>% mutate(dbSNP = ifelse(POS == 45379909, "rs2467850", 
                                                 ifelse(POS == 45466280, "rs55673230",
                                                       ifelse(POS == 30755467, "rs34830202", "None"))))

head(scr_gwas_x %>% filter(POS %in% c(45379909, 45466280, 30755467)))

head(scr_gwas_x)

g_scr = ggplot(scr_gwas_x, aes(x=XPosition, y=-log10(p.value))) +
            geom_hline(yintercept=7.3, linetype = "longdash", color = "darkred", linewidth = 1.5) +
            geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
            scale_color_manual(values = rep(c("blue3", "darkorange1"), 22)) +
            geom_label_repel(data=subset(scr_gwas_x, dbSNP == "rs2467850"), 
                       aes(XPosition,-log10(p.value),label=dbSNP), 
                         size=10, max.overlaps = 20, nudge_x = -350000000, nudge_y = 0.25,
                             color = "black", min.segment.length = 0)  + 
            geom_label_repel(data=subset(scr_gwas_x, dbSNP == "rs55673230"), 
                       aes(XPosition,-log10(p.value),label=dbSNP), 
                         size=10, max.overlaps = 20, nudge_x = 300000000, nudge_y = 0.5,
                             color = "black", min.segment.length = 0)  + 
            geom_label_repel(data=subset(scr_gwas_x, dbSNP == "rs34830202"), 
                       aes(XPosition,-log10(p.value),label=dbSNP), 
                         size=10, max.overlaps = 20, nudge_x = -300000000, nudge_y = 0.5,
                             color = "black", min.segment.length = 0)  + 
            labs(x = "Chromosome", y = bquote(-log[10](p-value))) +
            scale_x_continuous(label = scr_gwas_ggplot$CHR, breaks=scr_gwas_ggplot$center) +  
            scale_y_continuous(limits = c(0,18), breaks = c(0,2,4,6,8,10,12,14,16,18), 
                               labels = c(0,2,4,6,8,10,12,14,16,18)) +
            theme_few(base_size = 40) +
            theme(legend.position="none", axis.text.x = element_text(size = 26, vjust = 0.5, angle = 90),
                 axis.text.y = element_text(size = 26), plot.margin = margin(1,unit = "cm"))

ggsave(filename = "Data/Figures/Figure2/GWAS.SaigeV2.Creatinine.V5.png", plot = g_scr, device = "png", 
       width = 21, height = 9, dpi = 300)

g_scr

scr_gwas_x %>% arrange(p.value) %>% head()

scr_gwas_x %>% filter(POS %in% c(45379909,45466280))

## Write overall GWAS results.
dim(scr_gwas)
head(scr_gwas)
write.table(x = scr_gwas, file = "Data/Saige.V2/Step2/extractedChrAll.SaigeGWAS.tsv", sep = "\t", row.names = F, quote = F)

-log10(0.0000000000000000228)

## Run fine mapping












