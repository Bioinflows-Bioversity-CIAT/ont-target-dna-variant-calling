## Workflow overview

This workflow is designed to perform variant calling on target ONT-sequenced regions using the GATK HaplotypeCaller and Clair3. The input data for this workflow includes:
- A reference genomes in FASTA format
- A text file containing the regions of interest (rois)
- A text file containing the paths to the input FASTQ files

### Step 1: Reference genome directory

Is requiered create a tab-delimited text file indicating the assembly name (`ref_name`) and the path to the reference genome (`ref_path`), as shown below:

| ref_name | ref_path |
| --- | --- |
| hg38 | /path/to/hg38.fa |

Is important to note that the reference name should be unique and should not contain spaces. You can find an example in `.test/config/reference_units.txt`.

### Step 2: Regions of interest (rois) directory

You need to create a tab-delimited text file indicating the reference genome where the region is relative to (`ref_name`), the chromosome name (`chrom`), the start position of the region (`pos_i`), the end position of the region (`pos_e`).
| ref_name | chrom | pos_i | pos_e |
| --- | --- | --- | --- |
| hg38 | chr1 | 1000000 | 2000000 |

### Step 3: Fastq files directory

A tab-delimited text file is required to indicate the sample name (`sample_id`) and the path to the input FASTQ files (`fastq_path`), as shown below:
| sample_id | fastq_path |
| --- | --- |
| sample1 | /path/to/sample1.fastq.gz |
| sample2 | /path/to/sample2.fastq.gz |


### Step 4: Modify config.yaml

In `config/config.yaml`, you can find the parameters that can be modified according to your needs. The parameters include:
- `base_dir`: The base directory where the workflow will be executed.
- `sample_units_path`: The path to the text file containing the sample names and paths to the input FASTQ files.
- `reference_units_path`: The path to the text file containing the reference genome names and paths.
- `interval_units_path`: The path to the text file containing the regions of interest.
- `GATK` parameters
- `Clair3` parameters

