version 1.0
 

workflow process_acaf {
    input {
        File acaf_bgen_file_list
        File acaf_sample_file_list
        File saige_sample_file
        File plink_sample_file
        
        File step1_model_file
        File step1_variance_ratio_file
        File step0_sparse_grm_file
        File step0_sparse_grm_sample_file
        
        File space_list
        File cores_list
        File memory_list
        File plinkmemory_list
        
        File chromosome_list
        File output_prefix_list       
    }
    
    Array[File] acaf_bgen_files = read_lines(acaf_bgen_file_list)
    Array[File] acaf_sample_files = read_lines(acaf_sample_file_list)
    Array[File] chromosomes = read_lines(chromosome_list)

    Array[File] spaces = read_lines(space_list)
    Array[File] cores = read_lines(cores_list)
    Array[File] memorys = read_lines(memory_list)
    Array[File] plinkmemorys = read_lines(plinkmemory_list)

    Array[File] output_prefixes = read_lines(output_prefix_list)
    
    parameter_meta {
        acaf_bgen_files: {localization_optional: true}
        acaf_sample_files: {localization_optional: true}
        saige_sample_file: {localization_optional: true}
        plink_sample_file: {localization_optional: true}
        
        step1_model_file: {localization_optional: true}
        step1_variance_ratio_file: {localization_optional: true}
        step0_sparse_grm_file: {localization_optional: true}
        step0_sparse_grm_sample_file: {localization_optional: true}
    }
 
    scatter (idx in range(length(acaf_bgen_files))) {
        call downsample {
            input:
                bgen_file = acaf_bgen_files[idx],
                sample_file = acaf_sample_files[idx],
                output_prefix = output_prefixes[idx],
                chromosome = chromosomes[idx],

                step1_model_file = step1_model_file,
                step1_variance_ratio_file = step1_variance_ratio_file,
                step0_sparse_grm_file = step0_sparse_grm_file,
                step0_sparse_grm_sample_file = step0_sparse_grm_sample_file,

                space = spaces[idx],
                core = cores[idx],
                memory = memorys[idx],
                plinkmemory = plinkmemorys[idx],

                saige_sample_file = saige_sample_file,
                plink_sample_file = plink_sample_file
                
        }
    }
}
 
task downsample {
    input{
        File bgen_file
        File sample_file

        File saige_sample_file
        File plink_sample_file
        
        File step1_model_file
        File step1_variance_ratio_file
        File step0_sparse_grm_file
        File step0_sparse_grm_sample_file
        
        String space
        String core
        String memory
        String plinkmemory
        
        String output_prefix
        String chromosome

    }
 
    command <<<
        set -e
        echo "Working on "
        echo ~{bgen_file}
        
        # Make output directory.
        mkdir -p Plink
        mkdir -p Saige
        
        # Do the required QC for missingness, hwe and mac.
        plink2 --bgen ~{bgen_file} ref-first \
        --sample ~{sample_file} \
        --keep ~{plink_sample_file} \
        --geno 0.1 \
        --hwe 1e-15 \
        --maf 0.01 --mac 20 \
        --threads ~{core} --memory ~{plinkmemory} \
        --export bgen-1.3 ref-first \
        --out Plink/~{output_prefix}.ACAF.QC
        
        # Make BED files.
        plink2 --bgen Plink/~{output_prefix}.ACAF.QC.bgen ref-first \
        --sample Plink/~{output_prefix}.ACAF.QC.sample \
        --threads ~{core} --memory ~{plinkmemory} \
        --make-bed \
        --out Plink/~{output_prefix}.ACAF.QC.PlinkNative 
        
        ## Run GWAS step 2.
        step2_SPAtests.R --bedFile=Plink/~{output_prefix}.ACAF.QC.PlinkNative.bed --famFile=Plink/~{output_prefix}.ACAF.QC.PlinkNative.fam --bimFile=Plink/~{output_prefix}.ACAF.QC.PlinkNative.bim --AlleleOrder=alt-first --SAIGEOutputFile=Saige/~{output_prefix}.SaigeGWAS --chrom=~{chromosome} --LOCO=FALSE --minMAC=20 --sampleFile=~{saige_sample_file} --GMMATmodelFile=~{step1_model_file} --varianceRatioFile=~{step1_variance_ratio_file} --sparseGRMFile=~{step0_sparse_grm_file} --sparseGRMSampleIDFile=~{step0_sparse_grm_sample_file} --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_fastTest=TRUE --is_output_moreDetails=TRUE
    
    >>>
    
    output{
        File qced_bgen = "Plink/~{output_prefix}.ACAF.QC.bgen"
        File qced_sample = "Plink/~{output_prefix}.ACAF.QC.sample"
        File qced_log = "Plink/~{output_prefix}.ACAF.QC.log"
        
        File final_bed = "Plink/~{output_prefix}.ACAF.QC.PlinkNative.bed"
        File final_bim = "Plink/~{output_prefix}.ACAF.QC.PlinkNative.bim"
        File final_fam = "Plink/~{output_prefix}.ACAF.QC.PlinkNative.fam"
        File final_log = "Plink/~{output_prefix}.ACAF.QC.PlinkNative.log"
        
        File saige_step2 = "Saige/~{output_prefix}.SaigeGWAS"
        File saige_step2_index = "Saige/~{output_prefix}.SaigeGWAS.index"

        
    }
 
    runtime {
        docker:"gcr.io/docker-to-gcr-saige/lonely_grass:v4"
        memory: memory+"GB"
        cpu: core
        disks: "local-disk " + space + " HDD"
        bootDiskSizeGb: 75
    }
}


#*********Step 1**********#
step1_fitNULLGLMM.R --sparseGRMFile=Data/Saige/GRM/extractedChrAllPruned.CreatinineCohort.SparseMatrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=Data/Saige/GRM/extractedChrAllPruned.CreatinineCohort.SparseMatrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --plinkFile=Data/Saige/GRM/Genotype/extractedChrAllPruned.CreatinineCohort.Genotype. \
    --useSparseGRMtoFitNULL=TRUE \
    --phenoFile=Data/Saige/Phenotype/PhenotypeFile.tsv \
    --phenoCol=MedianCreatinine \
    --covarColList=SexCoded,Age,AgeSq,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --sexCol=SexCoded \
    --qCovarColList=SexCoded \
    --sampleIDColinphenoFile=SampleID \
    --traitType=quantitative \
    --invNormalize=TRUE \
    --SampleIDIncludeFile=Data/Saige/Cohort.Creatinine.tsv \
    --nThreads=8 \
    --outputPrefix=Data/Saige//Step1/SaigeResults.Creatinine
