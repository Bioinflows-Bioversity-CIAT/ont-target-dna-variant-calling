rule haplotype_caller_gvcf:
    input:
        # single or list of bam files
        bam = rules.rehead_bam_file.output.bam,
        bai = rules.rehead_bam_file.output.index,
        ref = rules.copy_reference.output,
        genome_dict = rules.create_dict.output,
        interval_list = rules.get_interval_list.output
    output:
        gvcf = f"{base_dir}/variant_calling/GATK/{{ref}}/HaplotypeCaller/intervals/{{chrom}}/{{sample_id}}_{{chrom}}.g.vcf.gz"
    log:
        f"{base_dir}/log/variant_calling/GATK/{{ref}}/HaplotypeCaller/intervals/{{chrom}}/{{sample_id}}_{{chrom}}.log"
    params:
        extra=get_GATK_HaplotypeCaller_params(),
        intervals = lambda wildcards, input: input.interval_list
    threads: 4
    resources:
        mem_mb=7000,
    wrapper:
        "v4.7.2/bio/gatk/haplotypecaller"

rule combine_by_sample_gvcfs:
    input:
        gvcfs = get_gvcfs_by_sample,
        ref = rules.copy_reference.output,
    output:
        gvcf = f"{base_dir}/variant_calling/GATK/{{ref}}/CombineGVCFs/{{sample_id}}.g.vcf.gz"
    log:
        f"{base_dir}/log/variant_calling/GATK/{{ref}}/CombineGVCFs/{{sample_id}}.log"
    params:
        extra = get_GATK_CombineGVCFs_params(),  
    resources:
        mem_mb=10240,
    wrapper:
        "v4.7.2/bio/gatk/combinegvcfs"


rule genomics_db_import:
    input:
        gvcfs=get_gvcfs_DB,
        interval_list = rules.get_interval_list.output
    output:
        db = directory(f"{base_dir}/variant_calling/GATK/{{ref}}/DB/{{chrom}}")
    log:
        f"{base_dir}/log/variant_calling/GATK/{{ref}}/DB/{{chrom}}.log"
    params:
        extra= get_GenomicsDBImport_params(),  # optional
        intervals = lambda wildcards, input: input.interval_list
    threads: 10
    resources:
        mem_mb=34000,
    wrapper:
        "v4.7.2/bio/gatk/genomicsdbimport"


rule genotype_gvcfs:
    input:
        genomicsdb = rules.genomics_db_import.output.db,
        ref=rules.copy_reference.output,
    output:
        vcf = f"{base_dir}/variant_calling/GATK/{{ref}}/DB/GenotypeGVCFs/{{chrom}}/{{interval_i}}-{{interval_e}}.vcf.gz"
    log:
        f"{base_dir}/log/variant_calling/GATK/{{ref}}/DB/GenotypeGVCFs/{{chrom}}/{{interval_i}}-{{interval_e}}.log"
    params:
        extra=get_GenotypeGVCFs_params(),
        intervals = lambda wildcards: "{chrom}:{interval_i}-{interval_e}".format(chrom = wildcards.chrom,
            interval_i = wildcards.interval_i,
            interval_e = wildcards.interval_e)
    resources:
        mem_mb=5120
    wrapper:
        "v4.7.2/bio/gatk/genotypegvcfs"


rule bcftools_concat:
    input:
        calls = get_interval_raw_vcfs,
        fai = f"{base_dir}/resources/{{ref}}/{{ref}}.fasta.fai"
    output:
        vcf = f"{base_dir}/variant_calling/GATK/{{ref}}/DB/GenotypeGVCFs/merged.vcf.gz"
    log:
        f"{base_dir}/log/variant_calling/GATK/{{ref}}/DB/GenotypeGVCFs/merged.log"
    params:
        uncompressed_bcf=False,
        extra="-a -Oz",  # optional parameters for bcftools concat (except -o)
    threads: 30
    resources:
        mem_mb=10240,
    wrapper:
        "v4.7.2/bio/bcftools/concat"
