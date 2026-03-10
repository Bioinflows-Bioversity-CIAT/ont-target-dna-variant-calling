# import basic packages
import pandas as pd
from snakemake.utils import validate
from collections import defaultdict
import glob
import json

def get_fastq_per_sample(wildcards):
    return sample_units.loc[wildcards.sample_id, "fastq_path"]

def get_reference_fasta(wildcards):
    return(reference_units.loc[wildcards.ref, "ref_path"])

def get_GATK_HaplotypeCaller_params():
    # Annotation params
    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['G']])
    kmers = ' '.join(["-kmer-size {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['kmer-size']])
    extra = ' '.join(["{p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['extra']])
    min_base_qual = "--min-base-quality-score {p}".format(p=config["GATK"]['HaplotypeCaller']['min_base-qual'])
    params = ' '.join([min_base_qual, annot, kmers, extra])
    return params

def get_GATK_CombineGVCFs_params():
    # Annotation params
    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['CombineGVCFs']['G']])
    return annot

def get_gvcfs_by_sample(wildcards):
    ref_intervals = interval_units[interval_units["ref_name"] == wildcards.ref]
    chroms = list(ref_intervals.chrom.unique())
    gvcfs = expand("{base_dir}/variant_calling/GATK/{ref}/HaplotypeCaller/intervals/{chrom}/{sample_id}_{chrom}.g.vcf.gz",
                   base_dir = base_dir,
                   chrom = chroms,
                   **wildcards)
    return gvcfs

def get_gvcfs_DB(wildcards):
    gvcfs_list = list()
    for isample, row in sample_units.iterrows():
        gvcf = "{base_dir}/variant_calling/GATK/{ref}/CombineGVCFs/{sample_id}.g.vcf.gz".format(
            base_dir = base_dir,
            sample_id = isample,
            **wildcards)
        gvcfs_list.append(gvcf)
    return gvcfs_list

def get_GenomicsDBImport_params():
    params = ""
    for param in config['GATK']['GenomicsDBImport'].keys():
        if type(config['GATK']['GenomicsDBImport'][param]) != list:
            params += "--{param} {value} ".format(param = param,
                value = config['GATK']['GenomicsDBImport'][param]
            )
        else:
            for option in config['GATK']['GenomicsDBImport'][param]:
                params += "{option} ".format(option = option)
    return params

def get_GenotypeGVCFs_params():
    params = ""
    for param in config['GATK']['GenotypeGVCFs'].keys():
        if param != "interval_length":
            if type(config['GATK']['GenotypeGVCFs'][param]) != list:
                params += "--{param} {value} ".format(param = param,
                    value = config['GATK']['GenotypeGVCFs'][param])
            else:
                for option in config['GATK']['GenotypeGVCFs'][param]:
                    params += "{option} ".format(option = option)
    return params

def get_interval_raw_vcfs(wildcards):
    vcfs = list()
    ref_intervals = interval_units[interval_units["ref_name"] == wildcards.ref]
    print(wildcards)
    for index, row in ref_intervals.iterrows():
        i_vcf = "{base_dir}/variant_calling/GATK/{ref}/DB/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz".format(
            base_dir = base_dir,
            chrom = row.chrom,
            interval_i = str(row.pos_i),
            interval_e = str(row.pos_e),
            **wildcards)
        vcfs.append(i_vcf)
    return vcfs

def get_final_GATK_vcfs():
    vcfs = [f"{base_dir}/variant_calling/GATK/{ref}/DB/GenotypeGVCFs/merged.vcf.gz" for ref in reference_units["ref_name"].tolist()]
    return vcfs
def get_final_Clair3_vcfs():
    vcfs = [f"{base_dir}/variant_calling/Clair3/{ref}/merged.vcf.gz" for ref in reference_units["ref_name"].tolist()]
    return vcfs

def get_clair3_gvcfs(wildcards):
    gvcfs_list = list()
    for isample, row in sample_units.iterrows():
        gvcf = "{base_dir}/variant_calling/Clair3/{ref}/{sample_id}/merge_output_roi.gvcf.gz".format(
            base_dir = base_dir,
            sample_id = isample,
            **wildcards)
        gvcfs_list.append(gvcf)
    return gvcfs_list

sample_units = (
    pd.read_csv(config["sample_units_path"], sep="\t", dtype={"sample_id": str})
    .set_index("sample_id", drop=False)
    .sort_index()
)

reference_units = (
    pd.read_table(config['reference_units_path'], sep = '\t')
    .set_index("ref_name", drop=False)
)

interval_units = pd.read_table(config["interval_units_path"], sep="\t")
interval_units["locus"] = interval_units.apply(lambda r: f"{r.chrom}:{r.pos_i}-{r.pos_e}", axis=1)
# validate sample sheet and config file
# validate(config, schema="../../config/schemas/config.schema.yaml")

validate(reference_units, schema="../../config/schemas/references.schema.yaml")


base_dir = config['base_dir']
print(sample_units)
print(reference_units)
print(interval_units)


wildcard_constraints:
    sample_id="|".join(sample_units['sample_id'].unique()),
    ref = "|".join(reference_units.index),
    chrom = "Chr[0-9]+|scaffold_[0-9]+",
    interval_i = "\d+",
    interval_e = "\d+",
