rule variant_calling_per_sample:
    input:
        bam = rules.rehead_bam_file.output.bam,
        ref = rules.copy_reference.output,
    output:
        gvcf = f"{base_dir}/variant_calling/Clair3/{{ref}}/{{sample_id}}/merge_output.gvcf.gz"
    conda:
        "envs/clair3.yml"
    threads:
        5
    shell:
        """
        run_clair3.sh --bam_fn={input.bam} \
        --ref_fn={input.ref} \
        --output={base_dir}/variant_calling/Clair3/{wildcards.ref}/{wildcards.sample_id} \
        --sample_name="{wildcards.sample_id}" \
        --threads={threads} \
        --gvcf \
        --platform="ont" \
        --print_ref_calls \
        --include_all_ctgs \
        --model_path={config[clair3][model]}
        """

rule filter_target_roi_clair3:
    input:
        gvcf = rules.variant_calling_per_sample.output,
        bed = rules.get_interval_bed.output
    output:
        gvcf = f"{base_dir}/variant_calling/Clair3/{{ref}}/{{sample_id}}/merge_output_roi.gvcf.gz"
    conda:
        "envs/ngs.yml"
    threads:
        1
    shell:
        """
        bcftools view \
        -R {input.bed} \
        -Oz -o {output.gvcf} \
        {input.gvcf} && \
        tabix -p vcf {output.gvcf}
        """

rule genomics_db_import_clair3:
    input:
        gvcfs=get_clair3_gvcfs,
        interval_list = rules.get_interval_list_clair3.output
    output:
        db = directory(f"{base_dir}/variant_calling/Clair3/{{ref}}/DB")
    log:
        f"{base_dir}/log/variant_calling/Clair3/{{ref}}/DB.log"
    params:
        extra= get_GenomicsDBImport_params(),  # optional
        intervals = lambda wildcards, input: input.interval_list
    threads: 10
    resources:
        mem_mb=34000,
    wrapper:
        "v4.7.2/bio/gatk/genomicsdbimport"

rule genotype_Clair3_gvcfs:
    input:
        genomicsdb = rules.genomics_db_import_clair3.output.db,
        ref=rules.copy_reference.output,
        interval_list = rules.get_interval_list_clair3.output
    output:
        vcf = f"{base_dir}/variant_calling/Clair3/{{ref}}/merged.vcf.gz"
    log:
        f"{base_dir}/log/variant_calling/Clair3/{{ref}}/merge.log"
    params:
        extra=get_GenotypeGVCFs_params(),
        intervals = lambda wildcards, input: input.interval_list
    resources:
        mem_mb=5120
    wrapper:
        "v4.7.2/bio/gatk/genotypegvcfs"

# rule merge_clair3_gvcfs:
#     input:
#         gvcfs = get_clair3_gvcfs,
#         bed = rules.get_interval_bed.output,
#     output:
#         temp(f"{base_dir}/variant_calling/Clair3/{{ref}}/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}/merged.bcf")
#     params:
#         config = "~/workflows/ont-target-dna-variant-calling/workflow/scripts/dna_nexus_clair3.yml"
#     threads:
#         1
#     conda:
#         "ngs"
#     shell:
#         """
#         glnexus_cli --config {params.config} \
#         -dir {base_dir}/variant_calling/Clair3/{wildcards.ref}/{wildcards.exp_name}/{wildcards.device_id}/{wildcards.flowcell_type}/{wildcards.flowcell_id}/{wildcards.run_id}/Gl.Nexus \
#         {input.gvcfs} > {output}
#         """

# rule compress_clair3_vcf:
#     input:
#         rules.merge_clair3_gvcfs.output
#     output:
#         f"{base_dir}/variant_calling/Clair3/{{ref}}/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}/merged.vcf.gz"
#     conda:
#         "ngs"
#     threads:
#         5
#     shell:
#         """
#         bcftools view {input} | bgzip -@ {threads} -c > {output}
#         """
