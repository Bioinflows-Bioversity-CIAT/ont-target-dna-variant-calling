rule trim_reads_se:
    input:
        get_sample_fastq,
    output:
        f"{base_dir}/trimming/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}/{{sample_id}}.fastq.gz"
    params:
        **config["trimmomatic"]["pe"]
    log:
        f"{base_dir}/log/trimming/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}/{{sample_id}}.log"
    threads:
        5
    resources:
        mem_mb=4096
    wrapper:
        "v4.7.2/bio/trimmomatic/se"
