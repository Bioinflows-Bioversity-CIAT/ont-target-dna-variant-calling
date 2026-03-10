rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

checkpoint genome_faidx:
    input:
        path = rules.copy_reference.output
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.fasta.fai"
    cache: True
    wrapper:
        "v4.7.2/bio/samtools/faidx"

rule get_interval_list:
    output:
        f"{base_dir}/resources/{{ref}}/intervals_{{chrom}}.list"
    run:
        i_interval = interval_units[(interval_units["ref_name"] == wildcards.ref) & (interval_units["chrom"] == wildcards.chrom)]
        i_interval["locus"].to_csv(output[0], index=False, header=False)

rule get_interval_list_clair3:
    output:
        f"{base_dir}/resources/{{ref}}/intervals.list"
    run:
        i_interval = interval_units[interval_units["ref_name"] == wildcards.ref]
        i_interval["locus"].to_csv(output[0], index=False, header=False)

rule get_interval_bed:
    output:
        f"{base_dir}/resources/{{ref}}/intervals.bed"
    run:
        i_interval = interval_units[interval_units["ref_name"] == wildcards.ref]
        i_interval[["chrom","pos_i","pos_e"]].to_csv(output[0], sep="\t", index=False, header=False)

rule create_dict:
    input:
        rules.copy_reference.output
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.dict"
    log:
        f"{base_dir}/log/reference/{{ref}}_dict.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v4.7.2/bio/picard/createsequencedictionary"
