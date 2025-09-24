version 1.0
 

workflow perform_annotation {
    input {
        File vcf_file_list
        File output_prefix_list
        File dbnsfp_txt_file
        File dbnsfp_tbi_file
    }
    
    Array[File] vcf_files = read_lines(vcf_file_list)
    Array[File] output_prefixes = read_lines(output_prefix_list)

    
    parameter_meta {
        vcf_files: {localization_optional: true}
        dbnsfp_txt_file: {localization_optional: true}
        dbnsfp_tbi_file: {localization_optional: true} 
    }
 
    scatter (idx in range(length(vcf_files))) {
        call annotation {
            input:
                vcf_file = vcf_files[idx],
                output_prefix = output_prefixes[idx],
                dbnsfp_txt_file = dbnsfp_txt_file,
                dbnsfp_tbi_file = dbnsfp_tbi_file
                

        }
    }
}
 
task annotation {
    input {
        File vcf_file
        File dbnsfp_txt_file
        File dbnsfp_tbi_file
        String output_prefix
    }

    command <<<
        set -e
        echo "Working on ~{vcf_file}"
        
        ### Make output directory.
        mkdir -p Annotation
        
        ### Annotate using SnpEff.
        java -Xmx16g -jar /usr/local/bin/snpEff.jar GRCh38.mane.1.2.refseq         ~{vcf_file} > Annotation/~{output_prefix}.Annotated.vcf
        
        echo "Finished annotating using SnpEff"
        
        ### Extract only missense variants. These variants will be annotated using SnpSift
        bcftools view -i 'INFO/ANN~"missense_variant"'         Annotation/~{output_prefix}.Annotated.vcf >         Annotation/~{output_prefix}.AnnotatedMissenseOnly.vcf
        
        echo "Created missense only VCF file"
    
        ### Annotate using DBNSFP
        java -Xmx16g -jar /usr/local/bin/SnpSift.jar         dbnsfp -noLog         -f SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,CADD_phred,MetaSVM_pred,MutationAssessor_pred,FATHMM_pred,MetaLR_pred         -v -db ~{dbnsfp_txt_file}         Annotation/~{output_prefix}.AnnotatedMissenseOnly.vcf >         Annotation/~{output_prefix}.DBNSFPAnnotated.vcf
    
    >>>
    
    output{
        File annotated_vcf = "Annotation/~{output_prefix}.Annotated.vcf"
        File annotated_missense_only_vcf = "Annotation/~{output_prefix}.AnnotatedMissenseOnly.vcf"
        File dbnsfp_annotated_vcf = "Annotation/~{output_prefix}.DBNSFPAnnotated.vcf"

    }
     
    runtime {
        docker:"gcr.io/docker-to-gcr-saige/lonely_grass:v3"
        memory: "160 GB"
        cpu: "24"
        disks: "local-disk 300 HDD"
        bootDiskSizeGb: 75
    }
}
 
