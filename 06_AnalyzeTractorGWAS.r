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
library('gridExtra')
library('ggrepel')
library('susieR')

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

#Implement a read function.
fread_df = function(x) {
    return(as.data.frame(fread(x))) 
}

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


read_tractor_chunks = function(tractor_result_files) { 
    # Read and store the files in a list
    tractor_files_read = lapply(tractor_result_files, fread_df)

    # Combine the list of data frames into one data frame
    tractor_results = do.call(rbind, tractor_files_read)

    #Reset index.
    rownames(tractor_results) = NULL
    
    ## Define the rows which need to be annotated.
    tractor_anno = tractor_results %>% filter(pval_anc1 <= 0.00000001)
    tractor_noanno = tractor_results %>% filter(pval_anc1 > 0.00000001) %>% 
                mutate(Gene = "None") 
    
    variant_lists <- split(tractor_anno, seq(nrow(tractor_anno)))
    GeneNames <- mclapply(variant_lists, find_gene, mane = mane, mc.cores = 8)

    # Combine the results back into DF1
    tractor_anno$Gene <- unlist(GeneNames)
    
    
    #Lable only the leading variant.
    tractor_anno = tractor_anno %>%    
                            group_by(Gene) %>%
                            mutate(Gene = ifelse(pval_anc1 == min(pval_anc1), Gene, "None")) %>%
                            ungroup()
    
    tractor_results = rbind(tractor_anno, tractor_noanno)
    
    
    return(tractor_results)
}

tractor_results = read_tractor_chunks(paste0('Data/Tractor/NewPaintings/Chunks/GWAS/TractorRunChunk_', 
                                                    c(1:16), 
                                                    '.tsv'))

dim(tractor_results)
head(tractor_results)
tractor_results %>% arrange(pval_anc0) %>% head()
tractor_results %>% arrange(pval_anc1) %>% head()
median(tractor_results$N)

tractor_results_0 = tractor_results %>% select(CHR, POS, ID, LAprop_anc0, 
                                               AF_anc0, se_anc0, beta_anc0, 
                                               REF, ALT,
                                               pval_anc0, Gene)
names(tractor_results_0) = c("Chr", "Pos", "ID", "LAprop_anc", 
                             "AF_anc", "se_anc", "beta_anc", 
                             "REF", "ALT",
                             "pval_anc", "Gene")
tractor_results_0$Ancestry = "European"

tractor_results_1 = tractor_results %>% select(CHR, POS, ID, LAprop_anc1, 
                                               AF_anc1, se_anc1, beta_anc1, 
                                               REF, ALT,
                                               pval_anc1, Gene)
names(tractor_results_1) = c("Chr", "Pos", "ID", "LAprop_anc", 
                             "AF_anc", "se_anc", "beta_anc",
                             "REF", "ALT",
                             "pval_anc", "Gene")
tractor_results_1$Ancestry = "African"

tractor_results_summary = rbind(tractor_results_0, tractor_results_1)
head(tractor_results_summary)

tractor_results_summary %>% filter(ID == "chr15:45379909:C:T")
tractor_results_summary %>% filter(ID == "chr15:45466280:A:G")

tractor_results_summary_african = tractor_results_summary %>% 
                                    filter(Ancestry == "African") %>% 
                                    filter(pval_anc < 0.000075) %>% 
                                    arrange(pval_anc) 

names(tractor_results_summary_african) = c("CHR", "POS", "MarkerID", "LAprop_anc",
                                          "AF_anc", "se_anc", "BETA", "Allele1",
                                          "Allele2", "p.value", "Gene", "Ancestry")
dim(tractor_results_summary_african)
head(tractor_results_summary_african)

write.table(x = tractor_results_summary_african, file = "Data/Tractor/NewPaintings/Chunks/TractorSummary.1E5.tsv",
           row.names = F, quote = F, sep = "\t")

g = ggplot(tractor_results_summary, aes(x=LAprop_anc, color=as.factor(Ancestry))) +
    facet_wrap(~as.factor(Ancestry),scales="free") +
    geom_density(linewidth = 1.5) +
    labs(x = "Proportional ancestry", y = "Density") +
    scale_color_manual(values = c("African" = "blue3", "European" = "darkorange1")) +
    theme_bw(base_size = 32) +
    theme(legend.position = "none")

ggsave(filename = "Data/Figures/GWAS.LocalAncestry.HaplotypeDistribution.png", plot = g, device = "png", 
       width = 16, height = 10, dpi = 300)
g

tractor_results_summary = tractor_results_summary %>% mutate(dbSNP = ifelse(Pos == 45379909, "rs2467850", 
                                                 ifelse(Pos == 45466280, "rs55673230",
                                                       ifelse(Pos == 30755467, "rs34830202", 
                                                              ifelse(Pos == 45403676, "rs1153848", "None")))))
                                                             
tractor_results_summary %>% arrange(pval_anc) %>% head()

ancestry_map <- list("African" = "African haplotypes", "European" = "European haplotypes")


ancestry_labeller <- function(variable,value){
  return(ancestry_map[value])
}

g1 = ggplot(tractor_results_summary, aes(x=Pos, y=-log10(pval_anc), color = as.factor(Ancestry))) +
    geom_hline(yintercept=7.3, linetype = "longdash", color = "darkred", linewidth = 1.5) +
    geom_point(alpha=0.75, size=2) +
    labs(x = "Position (Chromosome 15)", y = bquote(-log[10](p-value))) +
    scale_color_manual(values = c("African" = "blue3", "European" = "darkorange1")) +

    geom_label_repel(data=subset(tractor_results_summary, tractor_results_summary$dbSNP == "rs1153848"), 
                   aes(label=dbSNP), size=10, max.overlaps = 10, nudge_x = 20000000, nudge_y = 0.85, color = "black")  + 
    geom_label_repel(data=subset(tractor_results_summary, tractor_results_summary$dbSNP == "rs2467850"), 
                   aes(label=dbSNP), size=10, max.overlaps = 10, nudge_x = 20000000, nudge_y = -0.85, color = "black")  + 
    geom_label_repel(data=subset(tractor_results_summary, tractor_results_summary$dbSNP == "rs55673230" & Ancestry == "African"), 
                   aes(label=dbSNP), size=10, max.overlaps = 10, nudge_x = 25000000, nudge_y = 2.5, color = "black")  + 


    scale_x_continuous(breaks = c(20000000, 40000000, 60000000, 80000000, 100000000),
                       labels = c("20M", "40M", "60M", "80M", "100M")) +
    facet_wrap(~Ancestry, labeller = ancestry_labeller) +
    ylim(c(0,14)) + 
    theme_few(base_size = 40) +
    theme(legend.position="none", 
          legend.title = element_blank(), 
          axis.text = element_text(size = 26),
          plot.margin = margin(1,unit = "cm"))


ggsave(filename = "Data/Figures/Figure2/GWAS.Tractor.AncestrySpecificGWAS.png", plot = g1, device = "png", 
       width = 21, height = 9, dpi = 300)

g1

head(tractor_results_summary)
tractor_results_summary %>% arrange(pval_anc) %>% head()

write.table(file = "Data/Tractor/NewPaintings/Chunks/FineMapping/Loci.txt",
            x = tractor_results_summary %>% 
                        filter(Pos <= 45403676+500000 & 
                               Pos >= 45403676-500000 & 
                               Ancestry == "African" &
                               !is.na(beta_anc) &
                               AF_anc > 0.05 & 
                               AF_anc < 0.95) %>%
                        select(ID),
           row.names = F,
           col.names = F,
           quote = F,
           sep = "\t")

### Extract top loci.
left_end = 45403676-500000
right_end = 45403676+500000

tractor_results_summary_loci = tractor_results_summary %>% 
                                filter(Pos <= right_end & Pos >= left_end & Ancestry == "African")

dim(tractor_results_summary_loci)
head(tractor_results_summary_loci)

gene_tracks = mane %>% filter(CHR == 15 & 
                              ((left_end <= Start & Start <= right_end) | (left_end <= Stop & Stop <= right_end)))
gene_tracks = gene_tracks %>% mutate(YPOS = runif(n=nrow(gene_tracks), min=0.042, max=1))
gene_tracks

g4 = ggplot(tractor_results_summary_loci, 
            aes(x=Pos, y=-log10(pval_anc))) +
    geom_point(alpha=1, size=2.5, color = "blue3") +
    labs(x = "Chromosome", y = "-log10(p-value)") +
    geom_hline(yintercept=7.3) +
    geom_label_repel(data=subset(tractor_results_summary, tractor_results_summary$pval_anc < 0.0000000000005), 
                   aes(label=ID), size=4, max.overlaps = 20, nudge_x = 400000, color = "black")  + 
    ylim(c(0,14)) + xlim(c(44910000, 45910000)) +
    theme_few(base_size = 32) +
    theme(legend.position="bottom", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

g5 <- ggplot(gene_tracks) +
    geom_segment(aes(x = Start, xend = Stop, y = YPOS, yend = YPOS),
                 size = 10, color = "darkmagenta") +
    geom_text(data=gene_tracks,
            aes(Stop+1000,YPOS,label=Gene), vjust = 0.5, hjust = 0, size = 8) +     
    scale_y_continuous(limits = c(-0.05,1.05), breaks = c(-0.05,0.5,1.05), labels = c("00", "5", "0")) + 
    scale_x_continuous(limits = c(44910000, 45850000),
                       breaks = c(45000000, 45250000, 45500000, 45750000), 
                       labels = c("45.00M", "45.25M", "45.50M", "45.75M")) + 
    theme_bw(base_size = 32) +
    labs(x = "", y = "Tracks") +
    theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"))


combined_plot <- grid.arrange(g4, g5, ncol = 1, heights = c(3, 1.25),padding = unit(0, "cm"))
combined_plot

ggsave(filename = "Data/Figures/GWAS.Tractor.AfricanLocusZoom.png", plot = combined_plot, device = "png", 
       width = 18, height = 14, dpi = 300)

ld = as.data.frame(fread("Data/Tractor/NewPaintings/Chunks/FineMapping/extractedChr15.LD.ld"))
ld_variants = as.data.frame(fread("Data/Tractor/NewPaintings/Chunks/FineMapping/extractedChr15.LD.bim", header = F)) %>% pull(V2)
head(ld)
head(ld_variants)

ld.mat = matrix(as.vector(data.matrix(ld)), nrow=length(ld_variants), ncol=length(ld_variants))
tractor_fine_mapping <- tractor_results_summary %>%
                        filter(ID %in% ld_variants) %>%
                        mutate(ID = factor(ID, levels = ld_variants)) %>% 
                        arrange(ID) %>%
                        filter(!is.na(beta_anc) & Ancestry == "African")

dim(tractor_fine_mapping)
head(tractor_fine_mapping)

### Run SusieR
fitted_rss = susie_rss(bhat = tractor_fine_mapping$beta_anc, 
                       shat = tractor_fine_mapping$se_anc, 
                       R = ld.mat, 
                       max_iter = 300,
                       n = 18000,
                       L = 5)

print(fitted_rss$sets)

tractor_fine_mapping$variable = rownames(tractor_fine_mapping)
tractor_fine_mapping_result = merge(x=tractor_fine_mapping, 
                             y=summary(fitted_rss)$vars, 
             by="variable", all.x=TRUE)

head(tractor_fine_mapping_result)
tractor_fine_mapping_result %>% arrange(-cs) %>% head(20)

g4 = ggplot(tractor_fine_mapping_result, 
            aes(x=Pos, y=-log10(pval_anc))) +
    geom_point(alpha=1, size=2.5, color = "blue3") +
    geom_point(data=subset(tractor_fine_mapping_result, 
                           tractor_fine_mapping_result$cs!=-1), color = "limegreen", size=4) +
    labs(x = "Chromosome", y = "-log10(p-value)") +
    geom_hline(yintercept=7.3) +
    geom_label_repel(data=tractor_fine_mapping_result %>% filter(cs != -1), 
                   aes(label=ID), size=4, max.overlaps = 20, nudge_x = 400000, color = "black")  + 
    ylim(c(0,14)) + xlim(c(44910000, 45910000)) +
    theme_few(base_size = 32) +
    theme(legend.position="bottom", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

g5 <- ggplot(gene_tracks) +
    geom_segment(aes(x = Start, xend = Stop, y = YPOS, yend = YPOS),
                 size = 10, color = "darkmagenta") +
    geom_text(data=gene_tracks,
            aes(Stop+1000,YPOS,label=Gene), vjust = 0.5, hjust = 0, size = 8) +     
    scale_y_continuous(limits = c(-0.05,1.05), breaks = c(-0.05,0.5,1.05), labels = c("00", "5", "0")) + 
    scale_x_continuous(limits = c(44910000, 45850000),
                       breaks = c(45000000, 45250000, 45500000, 45750000), 
                       labels = c("45.00M", "45.25M", "45.50M", "45.75M")) + 
    theme_bw(base_size = 32) +
    labs(x = "", y = "Tracks") +
    theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"))


combined_plot <- grid.arrange(g4, g5, ncol = 1, heights = c(3, 1.25),padding = unit(0, "cm"))
combined_plot

ggsave(filename = "Data/Figures/GWAS.Tractor.AfricanLocusZoom.png", plot = combined_plot, device = "png", 
       width = 18, height = 14, dpi = 300)













