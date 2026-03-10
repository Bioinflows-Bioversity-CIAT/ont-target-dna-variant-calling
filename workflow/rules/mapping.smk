rule minimap:
    input:
        fastq = get_fastq_per_sample,
        ref = f"{base_dir}/resources/{{ref}}/{{ref}}.fasta"
    output:
        sam = temp(f"{base_dir}/mapping/{{ref}}/{{sample_id}}.sam")
    conda:
        "envs/ngs.yml"
    shell:
        """
        minimap2 -ax map-ont -Y {input.ref} {input.fastq} > {output.sam}
        """

rule sam_bam:
    input:
        sam = rules.minimap.output
    output:
        bam = temp(f"{base_dir}/mapping/{{ref}}/{{sample_id}}.bam"),
        bam_sort = temp(f"{base_dir}/mapping/{{ref}}/{{sample_id}}_sort.bam"),
        bam_index = temp(f"{base_dir}/mapping/{{ref}}/{{sample_id}}_sort.bam.bai")
    conda:
        "envs/ngs.yml"
    shell:
        """
        samtools view -S -F 0x100 -F 0x800 -b {input.sam} > {output.bam} && \
        samtools sort -o {output.bam_sort} {output.bam} && \
        samtools index {output.bam_sort}
        """ 

rule rehead_bam_file:
    input:
        bam = rules.sam_bam.output.bam_sort
    output:
        bam = f"{base_dir}/mapping/{{ref}}/{{sample_id}}_rehead.bam",
        index = f"{base_dir}/mapping/{{ref}}/{{sample_id}}_rehead.bam.bai"
    conda:
        "envs/ngs.yml"
    shell:
        """
        samtools addreplacerg \
        -r "@RG\tID:{wildcards.sample_id}\tSM:{wildcards.sample_id}\tPL:ONT" \
        --output-fmt BAM \
        -o {output.bam} {input.bam} && \
        samtools index {output.bam}
        """
