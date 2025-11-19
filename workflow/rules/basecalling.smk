rule dorado_multiplex_basecalling:
    input:
        run = get_run_pod5_dir
    output:
        fastq = f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}.fastq.gz"
    resources:
        gpu = 1
    log:
        f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}_basecalling.log"
    conda:
        "ngs"
    shell:
        """
        {config[dorado][path]} basecaller \
        --emit-fastq \
        --kit-name {config[dorado][barcode_kit]} \
        -x {config[dorado][device]} \
        -r \
        --min-qscore {config[dorado][min-q]} \
        --trim {config[dorado][trim]} \
        {config[dorado][model]} \
        {input.run} 2> {log} | bgzip > {output.fastq}
        """

rule generate_sample_sheet:
    output:
        f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}.sample_sheet.csv"
    run:
        exp_samples = sample_units[sample_units["exp_name"] == wildcards.exp_name]
        exp_report = reports[wildcards.exp_name][wildcards.run_id]["protocol_run_info"]
        exp_samples["experiment_id"] = wildcards.exp_name
        exp_samples["kit"] = exp_report["user_info"]["kit_info"]["sequencing_kit"]
        exp_samples["flow_cell_id"] = wildcards.flowcell_id
        exp_samples["position_id"] = wildcards.device_id
        exp_samples["protocol_run_id"] = wildcards.run_id
        exp_samples["flow_cell_product_code"] = exp_report["user_info"]["user_specified_product_code"]
        exp_samples["alias"] = exp_samples["sample_id"]
        out_cols = ["experiment_id", "kit", "flow_cell_id",
                    "position_id", "protocol_run_id", "flow_cell_product_code",
                    "alias", "barcode"]
        exp_samples[out_cols].to_csv(output[0], index=False)

checkpoint demultiplex_basecalling:
    input:
        fastq = rules.dorado_multiplex_basecalling.output.fastq,
        sample_sheet = rules.generate_sample_sheet.output
    output:
        directory(f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}")
    conda:
        "ngs"
    log:
        f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}_demux.log"
    shell:
        """
        {config[dorado][path]} demux \
        --output-dir {output} \
        --kit-name {config[dorado][barcode_kit]} \
        --emit-fastq \
        --emit-summary \
        --sample-sheet {input.sample_sheet} \
        --threads {threads} \
        {input.fastq} 2> {log}
        """
