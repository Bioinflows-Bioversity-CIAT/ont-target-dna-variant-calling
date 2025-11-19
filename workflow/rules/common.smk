# import basic packages
import pandas as pd
from snakemake.utils import validate
from collections import defaultdict
import glob
import json

# read sequencing_units sheet
sequencing_units = (
    pd.read_csv(config["sequencing_units_path"], sep="\t", dtype={"exp_name": str})
    .set_index("exp_name", drop=False)
    .sort_index()
)

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
validate(sequencing_units, schema="../../config/schemas/sequencing_units.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")
validate(reference_units, schema="../../config/schemas/references.schema.yaml")


base_dir = config['base_dir']
print(sequencing_units)
print(sample_units)
print(reference_units)
print(interval_units)

def get_report_file(directory):
    rep_dict = glob.glob(directory + "/report_*.json")
    if len(rep_dict) == 1:
        return rep_dict[0]
    elif len(rep_dict) > 1:
        print("Multiple report json files in folder: ", directory)
    else:
        rep_table = glob.glob(directory + "/final*.txt")
        if len(rep_table) == 1:
            return rep_table[0]
        else:
            print("Any report file found in folder:", directory)

def read_ont_report(directory):
    report_path = get_report_file(directory)
    if report_path[-4:] == 'json':
        with open(report_path, 'r') as file:
            # Load the JSON data
            rep_dict = json.load(file)
    else:
        rep_cust_dict = dict()
        with open(report_path, 'r') as file:
            for line in file:
                data = line.split('=')
                if len(data) > 1:
                    rep_cust_dict[data[0]] = data[1].replace('\n', '')
                else:
                    rep_cust_dict[data[0]] = ''
        flowcell_type = rep_cust_dict['protocol'].split(':')[2]
        # Interface when no json but final report
        rep_dict = {
            'protocol_run_info': {
                'run_id': rep_cust_dict['protocol_run_id'],
                'args': ['--pod5=off'],
                'device': {'device_id': rep_cust_dict['instrument']},
                'flow_cell': {'flow_cell_id': rep_cust_dict['flow_cell_id'], 'product_code': flowcell_type},
            }
        }

    return rep_dict

# Read reports using sequencing units datasheet
run_ids_list = list()
devices_list = list()
flowcell_types_list = list()
flowcell_ids_list = list()

reports = defaultdict(dict)
for i, row in sequencing_units.iterrows():
    i_report = read_ont_report(row['run_path'])
    reports[row['exp_name']][i_report['protocol_run_info']['run_id']] = i_report
    reports[row['exp_name']][i_report['protocol_run_info']['run_id']]['run_path'] = row['run_path']

    run_ids_list.append(i_report['protocol_run_info']['run_id'])
    devices_list.append(i_report['protocol_run_info']['device']['device_id'])
    flowcell_types_list.append(i_report['protocol_run_info']['flow_cell']['product_code'])
    flowcell_ids_list.append(i_report['protocol_run_info']['flow_cell']['flow_cell_id'])

def get_run_pod5_dir(wildcards):
    rep_dict = reports[wildcards.exp_name][wildcards.run_id]
    run_dir = rep_dict['run_path']
    if "--pod5=on" in rep_dict['protocol_run_info']['args']:
        return run_dir + '/pod5/'
    else:
            return f"{base_dir}/basecalling/{{exp_name}}/{{device_id}}/{{flowcell_type}}/{{flowcell_id}}/{{run_id}}_pod5".format(**wildcards)

def get_fastqs(base_dir):
    fastqs = list()
    for exp_name in reports.keys():
        for run_id in reports[exp_name].keys():
            rep_dict = reports[exp_name][run_id]
            device_id = rep_dict['protocol_run_info']['device']['device_id']
            flowcell_id = rep_dict['protocol_run_info']['flow_cell']['flow_cell_id']
            flowcell_type = rep_dict['protocol_run_info']['flow_cell']['product_code']
            output = f"{base_dir}/basecalling/{exp_name}/{device_id}/{flowcell_type}/{flowcell_id}/{run_id}/"
            fastqs.append(output)
    return fastqs

def get_sample_fastq(wildcards):
    checkpoint_output = checkpoints.demultiplex_basecalling.get(**wildcards).output
    return  f"{checkpoint_output}/{{run_id}}_{{sample_id}}.fastq".format(**wildcards)

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
    gvcfs = expand("{base_dir}/variant_calling/GATK/{ref}/HaplotypeCaller/intervals/{chrom}/{exp_name}/{device_id}/{flowcell_type}/{flowcell_id}/{run_id}/{sample_id}_{chrom}.g.vcf.gz",
                   base_dir = base_dir,
                   chrom = chroms,
                   **wildcards)
    return gvcfs

def get_gvcfs_DB(wildcards):
    checkpoint_output = checkpoints.demultiplex_basecalling.get(**wildcards).output[0]
    sample_list = glob.glob(checkpoint_output + "/*.fastq")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))

    gvcfs_list = list()
    for isample, row in sample_units.iterrows():
        gvcf = "{base_dir}/variant_calling/GATK/{ref}/CombineGVCFs/{exp_name}/{device_id}/{flowcell_type}/{flowcell_id}/{run_id}/{sample_id}.g.vcf.gz".format(
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
        i_vcf = "{base_dir}/variant_calling/GATK/{ref}/DB/GenotypeGVCFs/{exp_name}/{device_id}/{flowcell_type}/{flowcell_id}/{run_id}/{chrom}/{interval_i}-{interval_e}.vcf.gz".format(
            base_dir = base_dir,
            chrom = row.chrom,
            interval_i = str(row.pos_i),
            interval_e = str(row.pos_e),
            **wildcards)
        vcfs.append(i_vcf)
    return vcfs

def get_final_GATK_vcfs():
    vcfs = list()
    for exp_name in reports.keys():
        for run_id in reports[exp_name].keys():
            rep_dict = reports[exp_name][run_id]
            device_id = rep_dict['protocol_run_info']['device']['device_id']
            flowcell_id = rep_dict['protocol_run_info']['flow_cell']['flow_cell_id']
            flowcell_type = rep_dict['protocol_run_info']['flow_cell']['product_code']
            vcfs.extend([f"{base_dir}/variant_calling/GATK/{ref}/DB/GenotypeGVCFs/{exp_name}/{device_id}/{flowcell_type}/{flowcell_id}/{run_id}/merged.vcf.gz" for ref in reference_units["ref_name"].tolist()])
    return vcfs

wildcard_constraints:
    exp_name = "|".join(sequencing_units['exp_name'].unique()),
    device_id = "|".join(devices_list),
    flowcell_type = "|".join(flowcell_types_list),
    flowcell_id = "|".join(flowcell_ids_list),
    run_id = "|".join(run_ids_list),
    sample_id="|".join(sample_units['sample_id'].unique()),
    ref = "|".join(reference_units.index),
    chrom = "Chr[0-9]+|scaffold_[0-9]+",
    interval_i = "\d+",
    interval_e = "\d+",
