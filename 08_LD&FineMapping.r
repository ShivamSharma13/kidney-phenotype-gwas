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

### Extract top loci.
left_end = 45379909-250000
right_end = 45379909+250000

gene_tracks = mane %>% filter(CHR == 15 & 
                              ((left_end <= Start & Start <= right_end) | (left_end <= Stop & Stop <= right_end)))
gene_tracks = gene_tracks %>% mutate(YPOS = runif(n=nrow(gene_tracks), min=0.042, max=1))
gene_tracks

saige_gwas = as.data.frame(fread("Data/Saige.V2/Step2/extractedChr15.SaigeGWAS"))
dim(saige_gwas)
saige_gwas %>% arrange(p.value) %>% head()

lead_snp = "chr15:45379909:C:T"
output_ld_rsq = "Data/FineMapping/SaigeGWAS/extractedChr15.ACAF.RSqCorrelationsLeadSNP"

ld_rsq_command = paste("/usr/bin/plink2",
                       "--bfile",
                       "Data/Saige.V2/ACAF/extractedChr15.ACAF.QC.PlinkNative", 
                       "--r2-unphased", 
                       "--ld-snp",
                       lead_snp, 
                       "--ld-window-r2 0",
                       "--ld-window-kb 250",
                       "--out",
                       output_ld_rsq,
                       sep = " ")

s = system(ld_rsq_command, intern = TRUE)
s

ld_rsq = as.data.frame(fread("Data/FineMapping/SaigeGWAS/extractedChr15.ACAF.RSqCorrelationsLeadSNP.vcor"))
ld_rsq = merge(ld_rsq, saige_gwas, by.x = "ID_B", by.y = "MarkerID") %>% select(ID_B, CHR, POS, p.value, UNPHASED_R2)

#Add the lead variant back manually.
ld_rsq = rbind(ld_rsq %>% mutate(dbSNP = NA), 
               data.frame(ID_B = c("chr15:45379909:C:T"), CHR = c(15), 
                          POS = c(45379909), p.value = c(0.0000000000000000228), UNPHASED_R2 = c(1), dbSNP = "rs2467850"))

ld_rsq = ld_rsq %>%
            mutate(dbSNP = ifelse(POS == 45466280, "rs55673230", dbSNP))

# Add threshold cutoffs.
ld_rsq = ld_rsq %>% mutate(Tier = cut(UNPHASED_R2,
                               breaks = c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.01), 
                               labels = c("Tier5", "Tier4", "Tier3", "Tier2", "Tier1"),
                               right = FALSE))

dim(ld_rsq)
ld_rsq %>% arrange(-UNPHASED_R2) %>% head()

colors <- c(
  "Tier5" = "dodgerblue4",  # Dark Blue
  "Tier4" = "lightskyblue",  # Light Blue
  "Tier3" = "forestgreen",      # Yellow
  "Tier2" = "darkgoldenrod2",      # Orange
  "Tier1" = "red3"     # Dark Red
)

g4 = ggplot(ld_rsq, 
            aes(x=POS, y=-log10(p.value), fill = Tier)) +
    geom_hline(yintercept=7.3, linetype = "longdash", size = 1.5, color = "darkred") +
    geom_label_repel(data=ld_rsq %>% filter(!is.na(dbSNP)), 
                   aes(label=dbSNP), size=10, max.overlaps = 25, nudge_x = 120000, color = "black", fill = "white")  + 
    geom_point(size=6, pch = 21, color = "black", stroke = 0.5) +
    geom_point(data=subset(ld_rsq, 
                           ld_rsq$POS == 45379909), color = "black", 
               stroke = 1, fill = "#ffff3d", size=9, pch = 23) +
    labs(x = "Chromosome", y = bquote(-log[10](p-value))) +
    scale_fill_manual(values = colors) +  # Apply custom colors based on categories
    ylim(c(0,18)) + xlim(c(45129909, 45629909)) +
    theme_few(base_size = 40) +
    theme(legend.position="none", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, -0.35, 0), "cm"))

g5 <- ggplot(gene_tracks) +
    geom_segment(aes(x = Start, xend = Stop, y = YPOS, yend = YPOS),
                 linewidth = 10, color = "darkmagenta") +
    geom_text(data=gene_tracks,
            aes(Stop+1000,YPOS,label=Gene), vjust = 0.5, hjust = 0, size = 10) +     
    scale_y_continuous(limits = c(-0.05,1.05), breaks = c(-0.05,0.5,1.05), labels = c("00", "5", "0")) + 
    scale_x_continuous(limits = c(45129909, 45629909),
                       breaks = c(45100000, 45200000, 45300000, 45400000, 45500000, 45600000), 
                       labels = c("45.1M", "45.2M", "45.3M", "45.4M", "45.5M", "45.6M")) + 
    theme_bw(base_size = 40) +
    labs(x = "", y = "Tracks") +
    theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"), 
            axis.text.y = element_text(color = "white"),         # Change y-axis text color
            axis.ticks.y = element_line(color = "white"),        # Change y-axis ticks color
            axis.title.y = element_text(color = "white"))


combined_plot <- grid.arrange(g4, g5, ncol = 1, heights = c(3, 1.25),padding = unit(0, "cm"))
combined_plot

ggsave(filename = "Data/Figures/Figure3/GWAS.Saige.LocusZoom.LD.V2.pdf", plot = combined_plot, device = "pdf", 
       width = 18, height = 14, dpi = 300)


saige_gwas = as.data.frame(fread("Data/Saige.V2/Step2/extractedChr15.SaigeGWAS"))
dim(saige_gwas)
head(saige_gwas)

loci_file = "Data/FineMapping/SaigeGWAS/Loci.SaigeGWAS.txt"

write.table(file = loci_file,
            x = saige_gwas %>% 
                        filter(POS <= 45379909+250000 & 
                               POS >= 45379909-250000 & 
                               !is.na(BETA)) %>%
                        select(MarkerID),
           row.names = F,
           col.names = F,
           quote = F,
           sep = "\t")

output_ld_r = "Data/FineMapping/SaigeGWAS/extractedChr15.ACAF.LD"

ld_r_command = paste("/usr/bin/plink",
                       "--bfile",
                       "Data/Saige.V2/ACAF/extractedChr15.ACAF.QC.PlinkNative", 
                       "--keep-allele-order", 
                       "--extract",
                       loci_file, 
                       "--r square ",
                       "--make-bed",
                       "--out",
                       output_ld_r,
                       sep = " ")

s = system(ld_r_command, intern = TRUE)
s

#maybe delete the bed file?

ld = as.data.frame(fread("Data/FineMapping/SaigeGWAS/extractedChr15.ACAF.LD.ld"))
ld_variants = as.data.frame(fread("Data/FineMapping/SaigeGWAS/extractedChr15.ACAF.LD.bim", header = F)) %>% pull(V2)
head(ld)
head(ld_variants)

ld.mat = matrix(as.vector(data.matrix(ld)), nrow=length(ld_variants), ncol=length(ld_variants))
saige_fine_mapping <- saige_gwas %>%
                        filter(MarkerID %in% ld_variants) %>%
                        mutate(MarkerID = factor(MarkerID, levels = ld_variants)) %>% 
                        arrange(MarkerID) 

dim(saige_fine_mapping)
head(saige_fine_mapping)

### Run SusieR
fitted_rss_saige = susie_rss(bhat = saige_fine_mapping$BETA, 
                       shat = saige_fine_mapping$SE, 
                       R = ld.mat, 
                       n = 18979,
                       L = 10)
print(fitted_rss_saige$sets)

saige_fine_mapping$variable = rownames(saige_fine_mapping)
saige_fine_mapping_result = merge(x=saige_fine_mapping, 
                             y=summary(fitted_rss_saige)$vars, 
             by="variable", all.x=TRUE)

head(saige_fine_mapping_result)

# Annotate top hit variants.
saige_fine_mapping_result = saige_fine_mapping_result %>%
            mutate(dbSNP = NA) %>%
            mutate(dbSNP = ifelse(POS == 45379909, "rs2467850", dbSNP))

saige_fine_mapping_result %>% arrange(-cs) %>% head(20)

#Write susie results.
saige_fine_mapping_result %>% filter(cs != -1)
write.table(saige_fine_mapping_result %>% filter(cs != -1), 
            file = "Data/FineMapping/SaigeGWAS/SaigeGWAS.FineMappingResults.txt",
           sep = "\t", row.names = F, quote = F)

write.table(saige_fine_mapping_result %>% filter(cs != -1) %>% select(CHR, POS), 
            file = "Data/FineMapping/SaigeGWAS/SaigeGWAS.FineMappingResults.VariantIDS.txt",
           sep = "\t", row.names = F, quote = F, col.names = F)


g4 = ggplot(saige_fine_mapping_result, 
            aes(x=POS, y=-log10(p.value))) +
    geom_hline(yintercept=7.3, linetype = "longdash", size = 1.5, color = "darkred") +
    geom_label_repel(data=saige_fine_mapping_result %>% filter(!is.na(dbSNP)), 
                   aes(label=dbSNP), size=8, max.overlaps = 25, nudge_x = 120000, color = "black", fill = "white")  +
    geom_point(size=5, pch = 21, color = "black", stroke = 0.5, fill = "dodgerblue4") +
    geom_point(data=subset(saige_fine_mapping_result, 
                           saige_fine_mapping_result$cs!=-1), fill = "red3", size = 6, stroke = 0.5, pch = 24) +
    geom_point(data=subset(saige_fine_mapping_result, 
                           saige_fine_mapping_result$POS == 45379909), color = "black", 
               stroke = 1, fill = "#ffff3d", size=9, pch = 23) + 
    labs(x = "Chromosome", y = bquote(-log[10](p-value))) +
    ylim(c(0,18)) + xlim(c(45129909, 45629909)) +
    theme_few(base_size = 40) +
    theme(legend.position="bottom", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, -0.3, 0), "cm"))

g5 <- ggplot(gene_tracks) +
    geom_segment(aes(x = Start, xend = Stop, y = YPOS, yend = YPOS),
                 size = 10, color = "darkmagenta") +
    geom_text(data=gene_tracks,
            aes(Stop+1000,YPOS,label=Gene), vjust = 0.5, hjust = 0, size = 10) +     
    scale_y_continuous(limits = c(-0.05,1.05), breaks = c(-0.05,0.5,1.05), labels = c("00", "5", "0")) + 
    scale_x_continuous(limits = c(45129909, 45629909),
                       breaks = c(45100000, 45200000, 45300000, 45400000, 45500000, 45600000), 
                       labels = c("45.1M", "45.2M", "45.3M", "45.4M", "45.5M", "45.6M")) + 
    theme_bw(base_size = 40) +
    labs(x = "", y = "Tracks") +
    theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"), 
            axis.text.y = element_text(color = "white"),         # Change y-axis text color
            axis.ticks.y = element_line(color = "white"),        # Change y-axis ticks color
            axis.title.y = element_text(color = "white"))


combined_plot <- grid.arrange(g4, g5, ncol = 1, heights = c(3, 1.25),padding = unit(0, "cm"))
combined_plot

ggsave(filename = "Data/Figures/Figure3/GWAS.Saige.LocusZoom.Susie.V2.pdf", plot = combined_plot, device = "pdf", 
       width = 18, height = 14, dpi = 300)

######################################





















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

tractor_results_summary %>% filter(Pos %in% c(45466280,45379909))
tractor_results_summary %>% arrange(pval_anc) %>% head()

tractor_loci_file = "Data/FineMapping/Tractor/Loci.Tractor.txt"

write.table(file = tractor_loci_file,
            x = tractor_results_summary %>% 
                        filter(Pos <= 45403676+250000 & 
                               Pos >= 45403676-250000 & 
                               Ancestry == "African" &
                               !is.na(beta_anc)) %>%
                        select(ID),
           row.names = F,
           col.names = F,
           quote = F,
           sep = "\t")

### Extract top loci.
left_end = 45403676-250000
right_end = 45403676+250000

tractor_results_summary_loci = tractor_results_summary %>% 
                                filter(Pos <= right_end & Pos >= left_end & Ancestry == "African")

dim(tractor_results_summary_loci)
head(tractor_results_summary_loci)

gene_tracks = mane %>% filter(CHR == 15 & 
                              ((left_end <= Start & Start <= right_end) | (left_end <= Stop & Stop <= right_end)))
gene_tracks = gene_tracks %>% mutate(YPOS = runif(n=nrow(gene_tracks), min=0.042, max=1))
gene_tracks

output_ld_r = "Data/FineMapping/Tractor/extractedChr15.ACAF.LD"

ld_r_command = paste("/usr/bin/plink",
                       "--bfile",
                       "Data/Tractor/NewPaintings/ACAF/extractedChr15.Overlap", 
                       "--keep-allele-order", 
                       "--extract",
                       tractor_loci_file, 
                       "--r square ",
                       "--make-bed",
                       "--out",
                       output_ld_r,
                       sep = " ")

s = system(ld_r_command, intern = TRUE)
s

#maybe delete the bed file?

ld = as.data.frame(fread("Data/FineMapping/Tractor/extractedChr15.ACAF.LD.ld"))
ld_variants = as.data.frame(fread("Data/FineMapping/Tractor/extractedChr15.ACAF.LD.bim", header = F)) %>% pull(V2)
head(ld)
head(ld_variants)
tail(ld_variants)


ld.mat = matrix(as.vector(data.matrix(ld)), nrow=length(ld_variants), ncol=length(ld_variants))
tractor_fine_mapping <- tractor_results_summary %>%
                        filter(ID %in% ld_variants) %>%
                        mutate(ID = factor(ID, levels = ld_variants)) %>% 
                        arrange(ID) %>%
                        filter(Ancestry == "African")

dim(tractor_fine_mapping)
head(tractor_fine_mapping)
tail(tractor_fine_mapping)

### Run SusieR
fitted_rss_tractor = susie_rss(bhat = tractor_fine_mapping$beta_anc, 
                       shat = tractor_fine_mapping$se_anc, 
                       R = ld.mat, 
                       n = 18979,
                       L = 10)
print(fitted_rss_tractor$sets)

tractor_fine_mapping$variable = rownames(tractor_fine_mapping)
tractor_fine_mapping_result = merge(x=tractor_fine_mapping, 
                             y=summary(fitted_rss_tractor)$vars, 
             by="variable", all.x=TRUE)

head(tractor_fine_mapping_result)

# Annotate top hit variants.
tractor_fine_mapping_result = tractor_fine_mapping_result %>%
            mutate(dbSNP = NA) %>%
            mutate(dbSNP = ifelse(Pos == 45403676, "rs1153848", dbSNP))

tractor_fine_mapping_result %>% arrange(-cs) %>% head()

g4 = ggplot(tractor_fine_mapping_result, 
            aes(x=Pos, y=-log10(pval_anc))) +
    geom_label_repel(data=tractor_fine_mapping_result %>% filter(!is.na(dbSNP)), 
                   aes(label=dbSNP), size=8, max.overlaps = 25, nudge_x = 120000, color = "black", fill = "white")  +
    geom_hline(yintercept=7.3, linetype = "longdash", size = 1, color = "darkred") +
    geom_point(size=5, pch = 21, color = "black", stroke = 0.5, fill = "dodgerblue4") +
    geom_point(data=subset(tractor_fine_mapping_result, 
                           tractor_fine_mapping_result$cs!=-1), fill = "red3", size = 6, stroke = 0.5, pch = 24) +
    geom_point(data=subset(tractor_fine_mapping_result, 
                           !is.na(tractor_fine_mapping_result$dbSNP)), color = "black", 
               stroke = 1, fill = "#ffff3d", size=9, pch = 23) + 
    labs(x = "Chromosome", y = "-log10(p-value)") +
    ylim(c(0,18)) + xlim(c(45153676, 45653676)) +
    theme_few(base_size = 32) +
    theme(legend.position="bottom", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, -0.3, 0), "cm"))

g5 <- ggplot(gene_tracks) +
    geom_segment(aes(x = Start, xend = Stop, y = YPOS, yend = YPOS),
                 size = 10, color = "darkmagenta") +
    geom_text(data=gene_tracks,
            aes(Stop+1000,YPOS,label=Gene), vjust = 0.5, hjust = 0, size = 8) +     
    scale_y_continuous(limits = c(-0.05,1.05), breaks = c(-0.05,0.5,1.05), labels = c("00", "5", "0")) + 
    scale_x_continuous(limits = c(45153676, 45653676),
                       breaks = c(45100000, 45200000, 45300000, 45400000, 45500000, 45600000), 
                       labels = c("45.1M", "45.2M", "45.3M", "45.4M", "45.5M", "45.6M")) + 
    theme_bw(base_size = 32) +
    labs(x = "", y = "Tracks") +
    theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"), 
            axis.text.y = element_text(color = "white"),         # Change y-axis text color
            axis.ticks.y = element_line(color = "white"),        # Change y-axis ticks color
            axis.title.y = element_text(color = "white"))


combined_plot <- grid.arrange(g4, g5, ncol = 1, heights = c(3, 1.25),padding = unit(0, "cm"))
combined_plot

ggsave(filename = "Data/Figures/Figure3/GWAS.Tractor.LocusZoom.Susie.pdf", plot = combined_plot, device = "pdf", 
       width = 18, height = 14, dpi = 300)


