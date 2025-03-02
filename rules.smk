import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

primary_directory = "/Users/skulkarni/Desktop/Research/LR_data"


rule get_input_FASTQ_files:
    group:
        "sequential"
    output:
        fastq = f"{run_out_dir}/{{run_ID}}.fastq.gz",
    params:
        run_out_dir = run_out_dir,
    shell:       
        """
        # Download the FASTQ file
        fastq-dump {wildcards.run_ID} --outdir {params.run_out_dir} --gzip
        """


# rule trim_adapters_HiFi:
#     input:
#         fastq = f"{run_out_dir}/{{run_ID}}.fastq.gz",
#     output:
#         fastqc_directory = directory(f"{run_out_dir}/read_QC"), # need to just specify an output directory for fastqc
#         hifi_filt_fastq = f"{run_out_dir}/{{run_ID}}.filt.fastq.gz", # filtered FASTQ file without adapter sequences
#         hifi_contam = temp(f"{run_out_dir}/{{run_ID}}.contaminant.blastout"), # contaminated reads
#         hifi_stats = temp(f"{run_out_dir}/{{run_ID}}.stats"), # log file
#         hifi_blocklist = temp(f"{run_out_dir}/{{run_ID}}.blocklist"), # headers of reads to be removed
#     conda:
#         f"{primary_directory}/envs/read_QC.yaml"
#     params:
#         min_read_length = config["min_read_length"],
#         hifi_adapter_filt_script = "/Users/skulkarni/Desktop/git/HiFiAdapterFilt/hifiadapterfilt.sh",
#         hifi_adapter_filt_prefix = f"{run_out_dir}/{{run_ID}}",
#     shell:
#         """
#         # remove adapters. This creates some output files
#         bash {params.hifi_adapter_filt_script} -p {params.hifi_adapter_filt_prefix}
#         """



# rule trim_adapters_ONT:
#     input:
#         fastq = f"{run_out_dir}/{{run_ID}}.fastq.gz",
#     output:
#         fastqc_directory = directory(f"{run_out_dir}/read_QC"), # need to just specify an output directory for fastqc
#         hifi_filt_fastq = f"{run_out_dir}/{{run_ID}}.filt.fastq.gz", # filtered FASTQ file without adapter sequences
#         hifi_contam = temp(f"{run_out_dir}/{{run_ID}}.contaminant.blastout"), # contaminated reads
#         hifi_stats = temp(f"{run_out_dir}/{{run_ID}}.stats"), # log file
#         hifi_blocklist = temp(f"{run_out_dir}/{{run_ID}}.blocklist"), # headers of reads to be removed
#     params:
#         min_read_length = config["min_read_length"],
#         hifi_adapter_filt_script = "/Users/skulkarni/Desktop/git/HiFiAdapterFilt/hifiadapterfilt.sh",
#         hifi_adapter_filt_prefix = f"{run_out_dir}/{{run_ID}}",
#     shell:
#         """
#         # remove adapters with dorado (downloaded executable and added to path)
#         """


rule filter_low_quality_reads:
    input:
        fastq = f"{run_out_dir}/{{run_ID}}.fastq.gz",
    output:
        fastq_filtered = f"{run_out_dir}/read_QC/{{run_ID}}.filtered.fastq",
    conda:
        f"{primary_directory}/envs/read_QC.yaml"
    params:
        min_read_length = config["min_read_length"],
        min_mean_quality = config["min_mean_quality"],
        read_QC_dir = f"{run_out_dir}/read_QC",
    shell:
        """
        # remove short reads with filtlong
        # filtlong --min_length {params.min_read_length} --keep_percent 90 {input.fastq} > {output.fastq_filtered}
        
        # drop reads with an average quality score less than 10, which is a 10% percent error rate
        fastplong --length_required {params.min_read_length} --mean_qual {params.min_mean_quality} -i {input.fastq} -o {output.fastq_filtered}

        # run FASTQC to inspect stats after QC above
        fastqc -o {params.read_QC_dir} {output.fastq_filtered}

        # rm {input.fastq}
        """



# rule kraken_classification:
#     input:
#         fastq_filtered = f"{run_out_dir}/read_QC/{{run_ID}}.filtered.fastq",
#     output:
#         fastq_classified_filtered = f"{run_out_dir}/kraken/{{run_ID}}.kraken.filtered.fastq",
#         kraken_report = f"{run_out_dir}/kraken/kraken_report",
#         kraken_classifications = temp(f"{run_out_dir}/kraken/kraken_classifications"),
#     conda:
#         f"{primary_directory}/envs/variant_calling.yaml"
#     params:
#         kraken_db = f"{primary_directory}/{config['kraken_db']}",
#         output_dir = output_dir,
#     shell:
#         """
#         kraken2 --db {params.kraken_db} --threads 8 {input.fastq_filtered} --report {output.kraken_report} --classified-out {output.fastq_classified_filtered} > {output.kraken_classifications}

#         # rm {input.fastq_filtered}
#         """


# rule align_reads:
#     input:
#         fastq_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}.kraken.filtered.fastq",
#     output:
#         sam_file = f"{run_out_dir}/bam/{{run_ID}}.sam",
#         bam_file = f"{run_out_dir}/bam/{{run_ID}}.bam",
#         bam_index_file = f"{run_out_dir}/bam/{{run_ID}}.bam.bai",
#     conda:
#         f"{primary_directory}/envs/variant_calling.yaml"
#     params:
#         ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
#         long_read_source = config['long_read_source'], # hifi or ont, depending on the type of long-read sequencing
#     shell:
#         """
#         # index the reference FASTA file
#         minimap2 -d "{params.ref_genome}.mmi" {params.ref_genome}

#         # align reads to the reference genome using minimap2. The -a flag means to generate alignments in SAM format.
#         # the default is to output in PAF format.
#         minimap2 -x "map-{params.long_read_source}" -a "{params.ref_genome}.mmi" {input.fastq_trimmed_classified} -o {output.sam_file}

#         # sort alignment and convert to bam file
#         samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

#         # index alignment, which creates a .bai index file
#         samtools index {output.bam_file}
#         """


# # check depth
# rule get_BAM_file_depths_merge:
#     input:
#         bam_files = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
#     params:
#         ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
#         sample_out_dir = sample_out_dir,
#     output:
#         depth_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv"),
#         depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz",
#         merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
#         merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam.bai",
#     conda:
#         f"{primary_directory}/envs/variant_calling.yaml"
#     shell:
#         """
#         # get all runs associated with this sample_ID and compute depth
#         # -a computes depth at all positions, not just those with non-zero depth
#         # -Q is for minimum mapping quality: use 1, so that multiply mapped reads aren't counted. These have mapping quality of 0
#         samtools depth -a -Q 1 {input.bam_files} > {output.depth_file}

#         # get the length of the reference genome
#         genome_length=$(tail -n +2 {params.ref_genome} | tr -d '\n' | wc -c) # remove first line (FASTA header) and newline characters, then count characters to get ref genome length

#         # when there are multiple bam files, each one is its own column in the depth file.
#         num_sites_H37Rv=$(wc -l {output.depth_file} | awk '{{print $1}}')
    
#         if [ ! "$num_sites_H37Rv" -eq "$genome_length" ]; then
#             echo "Check that all $genome_length sites in the H37Rv reference genome are in {output.depth_file}, which currently has $num_sites_H37Rv sites"
#             exit 1
#         fi

#         # compress for space
#         gzip -c {output.depth_file} > {output.depth_file_gzip}

#         # merge them using samtools. works because the original bam files were sorted after converting from SAM files
#         samtools merge {input.bam_files} -o {output.merged_bam_file}
        
#         # index the merged BAM file for variant calling
#         samtools index {output.merged_bam_file}

#         # delete the original BAM files for each run and the index files
#         for bam_file in {input.bam_files}; do
#             if [ -e "$bam_file" ]; then
#                 rm "$bam_file"
#                 rm "$bam_file.bai"
#             fi
#         done
#         """


# # variant calling with pilon
# rule pilon_variant_calling:
#     input:
#         merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
#         merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam.bai",
#     output:
#         vcf_file = temp(f"{sample_out_dir}/variant_calling/{{sample_ID}}.vcf"),
#         vcf_file_gzip = f"{sample_out_dir}/variant_calling/{{sample_ID}}.vcf.gz",
#         vcf_file_variants_only = f"{sample_out_dir}/variant_calling/{{sample_ID}}.variants.vcf",
#         fasta_file = temp(f"{sample_out_dir}/variant_calling/{{sample_ID}}.fasta"),        
#     params:
#         ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
#         sample_pilon_dir = f"{sample_out_dir}/variant_calling",
#     conda:
#         f"{primary_directory}/envs/variant_calling.yaml"
#     shell:
#         """
#         pilon -Xmx30g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant --threads 8
            
#         # then gzip the full VCF file and delete the unzipped version. Also delete the FASTA file because it's not needed
#         gzip -c {output.vcf_file} > {output.vcf_file_gzip}

#         # save the variants only (non-REF calls) to another VCF file. Left-align indels and deduplicate variants with the same POS, REF, and ALT.
#         # this affects those cases where the position of the indel is ambiguous
#         # however, because of the shifting positions, the position of the indel can change, so need to sort it again
#         bcftools norm --rm-dup none --fasta-ref {params.ref_genome} {output.vcf_file_gzip} | bcftools sort | bcftools view --types snps,indels,mnps,other > {output.vcf_file_variants_only}
#         """


# # structural variant detection with sniffles
# rule sniffles_variant_calling:
#     input:
#         merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
#         merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam.bai",
#     output:
#         vcf_file = f"{sample_out_dir}/variant_calling/{{sample_ID}}.SV.vcf",  
#     params:
#         ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
#     conda:
#         f"{primary_directory}/envs/variant_calling.yaml"
#     shell:
#         """
#         sniffles --input {input.merged_bam_file} --vcf {output.vcf_file} --reference {params.ref_genome} --threads 8
#         """



# rule write_variants_to_TSV:
#     input:
#         pilon_vcf = f"{sample_out_dir}/variant_calling/{{sample_ID}}.variants.vcf",
#         sniffles_vcf = f"{sample_out_dir}/variant_calling/{{sample_ID}}.SV.vcf",
#     output:
#         pilon_TSV = f"{sample_out_dir}/variant_calling/{{sample_ID}}.variants.tsv",  
#         sniffles_TSV = f"{sample_out_dir}/variant_calling/{{sample_ID}}.SV.tsv"
#     conda:
#         f"{primary_directory}/envs/variant_annotation.yaml"
#     shell:
#         """
#         # YAY, SnpSift adds headers to the TSV files
#         # remove ANN because haven't run annotation on it yet
#         SnpSift extractFields {input.pilon_vcf} POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC -e "" > {output.pilon_TSV}

#         # SUPPORT is read depth supporting the SV
#         # COVERAGE field gives depth at the upstream, start, center, end, downstream of the SV
#         # STRAND is the strand senses (+ and -) that support the SV. Useful to determine if there is strand bias if only one strand supports the SV
#         SnpSift extractFields {input.sniffles_vcf} POS REF ALT FILTER QUAL IMPRECISE VAF SUPPORT COVERAGE STRAND -e "" > {output.sniffles_TSV}
#         """