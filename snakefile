import numpy as np
import pandas as pd

configfile: "config.yaml"

# Define the list of samples, either directly or through config
output_dir = config["output_dir"]

# dataframe of isolates to run. First column is the sample ID, second column is a comma-separated string of the sequencing run IDs
df_samples_runs = pd.read_csv(config['isolates_to_run'], sep='\t', header=None)
sample_run_dict = dict(zip(df_samples_runs[0], df_samples_runs[1].str.split(',')))

include: "rules.smk"

# Define helper functions to construct paths dynamically
def sample_out_dir(sample_ID):
    return f"{output_dir}/{sample_ID}"

def run_out_dir(sample_ID, run_ID):
    return f"{output_dir}/{sample_ID}/{run_ID}"

# Define rules in the Snakefile

rule all:
    input:
        [f"{output_dir}/{sample_ID}/{run_ID}/{run_ID}.fastq.gz" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        [f"{output_dir}/{sample_ID}/{run_ID}/read_QC/{run_ID}.filtered.fastq" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        # [f"{output_dir}/{sample_ID}/variant_calling/{sample_ID}.variants.vcf" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/variant_calling/{sample_ID}.variants.tsv" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/{run_ID}/bam/{run_ID}.bam.bai" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        # [f"{output_dir}/{sample_ID}/bam/{sample_ID}.depth.tsv.gz" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/{run_ID}/fastlin/output.txt" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        # [f"{output_dir}/{sample_ID}/{run_ID}/fastplong/{run_ID}.trimmed.fastq" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]]
        # [f"{output_dir}/{sample_ID}/{run_ID}/kraken/kraken_report" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]]
        # [f"{output_dir}/{sample_ID}/{run_ID}/fastlin/primary_lineage.txt" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        # [f"{output_dir}/{sample_ID}/bam/{sample_ID}.bam" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/lineage/F2_Coll2014.txt" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/WHO_resistance/{sample_ID}_pred_AF_thresh_75.csv" for sample_ID in sample_run_dict.keys()],