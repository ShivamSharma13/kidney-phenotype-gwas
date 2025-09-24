rm(list=ls())
set.seed(13)

options(scipen=100, digits=3)

library('dplyr')
library('data.table')
library('parallel')
library('tidyr')

options(repr.matrix.max.cols=50, repr.matrix.max.rows=500)

process_chr <- function(chr){

    ##Read consequence annotations.
    consequence_annotations <- as.data.frame(fread(paste0('Data/Annotation/African/extractedChr',
                                                chr,
                                                '.AfricanCohort.Annotated.vcf'), skip = 12))
    names(consequence_annotations)[1] = "CHR"
    dim(consequence_annotations)
    head(consequence_annotations)

    ##Read dbnsfp annotations.
    dbnsfp_annotations <- as.data.frame(fread(paste0('Data/Annotation/African/extractedChr',
                                                     chr,
                                                     '.AfricanCohort.DBNSFPAnnotated.vcf'), skip = 21))
    names(dbnsfp_annotations)[1] = "CHR"
    dbnsfp_annotations$INFO <- gsub(".*\\|;", "", dbnsfp_annotations$INFO)
    dim(dbnsfp_annotations)
    head(dbnsfp_annotations)

    consequence_annotations_extracted = consequence_annotations %>% 
                            separate(INFO, into = paste0("col", 1:16), sep = "\\|") %>%
                            select("ID", "col2","col3", "col4")

    names(consequence_annotations_extracted)[c(2,3,4)] = c("Consequence", "Impact", "Gene")
    head(consequence_annotations_extracted)

    ##Define the impact column names.
    dbnsfp_resources <- c("dbNSFP_CADD_phred", "dbNSFP_LRT_pred", 
                          "dbNSFP_MetaSVM_pred", "dbNSFP_MutationTaster_pred", 
                          "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_SIFT_pred")

    # Split the Column_with_info column based on the delimiter ";"
    split_info <- strsplit(dbnsfp_annotations$INFO, ";")

    # Extract values for each column
    for (i in seq_along(dbnsfp_resources)) {
      dbnsfp_annotations[[dbnsfp_resources[i]]] <- sapply(split_info, function(x) {
        val <- grep(paste0(dbnsfp_resources[i], "="), x, value = TRUE)
        if (length(val) == 0) 
            return("None") # If no matching value found, return NA
        else 
            return(sub(paste0(dbnsfp_resources[i], "="), "", val))
      })
    }

    #Remove the clumsy INFO column.
    dbnsfp_annotations = dbnsfp_annotations %>% select(-INFO)
    head(dbnsfp_annotations)

    #Assign annotations.
    dbnsfp_annotations = dbnsfp_annotations %>% 
                            mutate(DBNSFP_Consequence = case_when(
                                    (dbNSFP_LRT_pred == "D" & dbNSFP_MutationTaster_pred == "D" &
                                        dbNSFP_Polyphen2_HDIV_pred == "D" & dbNSFP_Polyphen2_HVAR_pred == "D" &
                                        dbNSFP_SIFT_pred == "D") ~ "LikelyDeleterious",
                                    (dbNSFP_LRT_pred == "D" | dbNSFP_MutationTaster_pred == "D" |
                                        dbNSFP_Polyphen2_HDIV_pred == "D" | dbNSFP_Polyphen2_HVAR_pred == "D" |
                                        dbNSFP_SIFT_pred == "D") ~ "PossiblyDeleterious",
                                    (dbNSFP_LRT_pred != "D" & dbNSFP_MutationTaster_pred != "D" &
                                        dbNSFP_Polyphen2_HDIV_pred != "D" & dbNSFP_Polyphen2_HVAR_pred != "D" &
                                        dbNSFP_SIFT_pred != "D") ~ "LikelyBenign"))

    dbnsfp_annotations = dbnsfp_annotations %>% select(ID, DBNSFP_Consequence)
    head(dbnsfp_annotations)

    #Merge the annotations but make sure to bring all from others.
    dim(consequence_annotations_extracted)
    merged_annotations = merge(consequence_annotations_extracted, dbnsfp_annotations, by = "ID", all.x = T)
    dim(merged_annotations)
    head(merged_annotations)

    #Make final annotation column and remove unnecessary variants.
    merged_annotations = merged_annotations %>% 
            filter(Impact == "HIGH" | 
                   DBNSFP_Consequence %in% c("PossiblyDeleterious", "LikelyDeleterious", "LikelyBenign")) %>%
            mutate(FinalAnnotation = ifelse(Impact == "HIGH", "pLOF", DBNSFP_Consequence)) %>%
            select(ID, Gene, FinalAnnotation)

    #And finally process one gene at a time.
    for (gene_name in unique(merged_annotations$Gene)) {
        this_gene_data = merged_annotations %>% filter(Gene == gene_name)

        #Extract the vector.
        id_vector = c(gene_name, "var", this_gene_data %>% pull(ID))
        impact_vector = c(gene_name, "anno", this_gene_data %>% pull(FinalAnnotation))

        #Check the lengths of these vectors.
        if (length(id_vector) != length(impact_vector)) {
            print(paste("Something is wrong with gene", gene_name))
        }

        #Convert vectors to a string.
        id_vector = paste(id_vector,collapse="\t")
        impact_vector = paste(impact_vector,collapse="\t")

        group_file_path = paste0("Data/Group/African/Group.chr",
                                 chr,
                                 ".tsv")

        #Append it to the file.
        write(id_vector,file=group_file_path,append=TRUE)
        write(impact_vector,file=group_file_path,append=TRUE)

        #Finishing touch!
        #print(paste("Finished with gene", gene_name))

    }
    
}

chrs = seq(1,22)
#chrs = c("X")

ress = mclapply(chrs, function(i) process_chr(i), mc.cores = 22)

process_chr <- function(chr){

    
    ##Read consequence annotations.
    consequence_annotations <- as.data.frame(fread(paste0('Data/Annotation/European/extractedChr',
                                                chr,
                                                '.EuropeanCohort.Annotated.vcf'), skip = 12))
    names(consequence_annotations)[1] = "CHR"
    dim(consequence_annotations)
    head(consequence_annotations)

    ##Read dbnsfp annotations.
    dbnsfp_annotations <- as.data.frame(fread(paste0('Data/Annotation/European/extractedChr',
                                                     chr,
                                                     '.EuropeanCohort.DBNSFPAnnotated.vcf'), skip = 21))
    names(dbnsfp_annotations)[1] = "CHR"
    dbnsfp_annotations$INFO <- gsub(".*\\|;", "", dbnsfp_annotations$INFO)
    dim(dbnsfp_annotations)
    head(dbnsfp_annotations)

    consequence_annotations_extracted = consequence_annotations %>% 
                            separate(INFO, into = paste0("col", 1:16), sep = "\\|") %>%
                            select("ID", "col2","col3", "col4")

    names(consequence_annotations_extracted)[c(2,3,4)] = c("Consequence", "Impact", "Gene")
    head(consequence_annotations_extracted)

    ##Define the impact column names.
    dbnsfp_resources <- c("dbNSFP_CADD_phred", "dbNSFP_LRT_pred", 
                          "dbNSFP_MetaSVM_pred", "dbNSFP_MutationTaster_pred", 
                          "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_SIFT_pred")

    # Split the Column_with_info column based on the delimiter ";"
    split_info <- strsplit(dbnsfp_annotations$INFO, ";")

    # Extract values for each column
    for (i in seq_along(dbnsfp_resources)) {
      dbnsfp_annotations[[dbnsfp_resources[i]]] <- sapply(split_info, function(x) {
        val <- grep(paste0(dbnsfp_resources[i], "="), x, value = TRUE)
        if (length(val) == 0) 
            return("None") # If no matching value found, return NA
        else 
            return(sub(paste0(dbnsfp_resources[i], "="), "", val))
      })
    }

    #Remove the clumsy INFO column.
    dbnsfp_annotations = dbnsfp_annotations %>% select(-INFO)
    head(dbnsfp_annotations)

    #Assign annotations.
    dbnsfp_annotations = dbnsfp_annotations %>% 
                            mutate(DBNSFP_Consequence = case_when(
                                    (dbNSFP_LRT_pred == "D" & dbNSFP_MutationTaster_pred == "D" &
                                        dbNSFP_Polyphen2_HDIV_pred == "D" & dbNSFP_Polyphen2_HVAR_pred == "D" &
                                        dbNSFP_SIFT_pred == "D") ~ "LikelyDeleterious",
                                    (dbNSFP_LRT_pred == "D" | dbNSFP_MutationTaster_pred == "D" |
                                        dbNSFP_Polyphen2_HDIV_pred == "D" | dbNSFP_Polyphen2_HVAR_pred == "D" |
                                        dbNSFP_SIFT_pred == "D") ~ "PossiblyDeleterious",
                                    (dbNSFP_LRT_pred != "D" & dbNSFP_MutationTaster_pred != "D" &
                                        dbNSFP_Polyphen2_HDIV_pred != "D" & dbNSFP_Polyphen2_HVAR_pred != "D" &
                                        dbNSFP_SIFT_pred != "D") ~ "LikelyBenign"))

    dbnsfp_annotations = dbnsfp_annotations %>% select(ID, DBNSFP_Consequence)
    head(dbnsfp_annotations)

    #Merge the annotations but make sure to bring all from others.
    dim(consequence_annotations_extracted)
    merged_annotations = merge(consequence_annotations_extracted, dbnsfp_annotations, by = "ID", all.x = T)
    dim(merged_annotations)
    head(merged_annotations)

    #Make final annotation column and remove unnecessary variants.
    merged_annotations = merged_annotations %>% 
            filter(Impact == "HIGH" | 
                   DBNSFP_Consequence %in% c("PossiblyDeleterious", "LikelyDeleterious", "LikelyBenign")) %>%
            mutate(FinalAnnotation = ifelse(Impact == "HIGH", "pLOF", DBNSFP_Consequence)) %>%
            select(ID, Gene, FinalAnnotation)

    #And finally process one gene at a time.
    for (gene_name in unique(merged_annotations$Gene)) {
        this_gene_data = merged_annotations %>% filter(Gene == gene_name)

        #Extract the vector.
        id_vector = c(gene_name, "var", this_gene_data %>% pull(ID))
        impact_vector = c(gene_name, "anno", this_gene_data %>% pull(FinalAnnotation))

        #Check the lengths of these vectors.
        if (length(id_vector) != length(impact_vector)) {
            print(paste("Something is wrong with gene", gene_name))
        }

        #Convert vectors to a string.
        id_vector = paste(id_vector,collapse="\t")
        impact_vector = paste(impact_vector,collapse="\t")

        group_file_path = paste0("Data/Group/European/Group.chr",
                                 chr,
                                 ".tsv")

        #Append it to the file.
        write(id_vector,file=group_file_path,append=TRUE)
        write(impact_vector,file=group_file_path,append=TRUE)

        #Finishing touch!
        #print(paste("Finished with gene", gene_name))

    }
    
}

chrs = seq(1,22)
#chrs = c("X")

ress = mclapply(chrs, function(i) process_chr(i), mc.cores = 22)


