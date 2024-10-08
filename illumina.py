#!/usr/bin/env python

import os
import pandas as pd
import concurrent.futures

_author_ = 'Khadim GUEYE'

workflow_path = "/hps/nobackup/cochrane/ena/users/khadim/project/ena-covid-systematic-analysis-workflow"

params = {
    "SARS2_FA": os.path.join(workflow_path, "data/NC_063383.1.fasta"),
    "SARS2_FA_FAI": os.path.join(workflow_path, "data/NC_063383.1.fasta.fai"),
    "SECRETS": os.path.join(workflow_path, "prepro/projects_accounts.csv"),
    "INDEX": os.path.join(workflow_path, "prepro/illumina.index.tsv"),
    "STOREDIR": os.path.join(workflow_path, "prepro/storeDir"),
    "OUTDIR": os.path.join(workflow_path, "prepro/results"),
    "STUDY": 'PRJEB55834'
}

def map_to_reference(run_accession, sample_accession, input_file_1, input_file_2, sars2_fasta, sars2_fasta_fai, projects_accounts_csv, study_accession , num_threads , STOREDIR , OUTDIR , workflow_path):
    # Retrieve FTP credentials from projects_accounts_csv
    line = os.popen(f"grep {study_accession} {projects_accounts_csv}").read()
    ftp_id = line.split(',')[2]
    ftp_password = line.split(',')[4]
    
    # Download FASTQ files using wget
    if ftp_id == 'public':
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_1.fastq.gz {input_file_1}")
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_2.fastq.gz {input_file_2}")
    else:
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_1.fastq.gz {input_file_1} --user={ftp_id} --password={ftp_password}")
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_2.fastq.gz {input_file_2} --user={ftp_id} --password={ftp_password}")
    
    # Execute Trimmomatic to perform read trimming
    os.system(f"cd {STOREDIR}; trimmomatic PE {run_accession}_1.fastq.gz {run_accession}_2.fastq.gz {run_accession}_trim_1.fq {run_accession}_trim_1_un.fq {run_accession}_trim_2.fq {run_accession}_trim_2_un.fq -summary {run_accession}_trim_summary -threads {num_threads} SLIDINGWINDOW:5:30 MINLEN:50")

    # Index the reference genome using BWA
    os.system(f"cd {STOREDIR}; bwa index {sars2_fasta}")
    
    # Map reads to the reference genome using BWA-MEM
    os.system(f"cd {STOREDIR}; bwa mem -t {num_threads} {sars2_fasta} {run_accession}_trim_1.fq {run_accession}_trim_2.fq | samtools view -bF 4 - | samtools sort - > {run_accession}_paired.bam")
    os.system(f"cd {STOREDIR}; bwa mem -t {num_threads} {sars2_fasta} {run_accession}_trim_1_un.fq {run_accession}_trim_2_un.fq | samtools view -bF 4 - | samtools sort - > {run_accession}_unpaired.bam")
    os.system(f"cd {STOREDIR}; samtools merge {run_accession}.bam {run_accession}_paired.bam {run_accession}_unpaired.bam")
    os.system(f"cd {STOREDIR}; rm {run_accession}_paired.bam {run_accession}_unpaired.bam")
    
    # Generate pileup file
    os.system(f"cd {STOREDIR}; samtools mpileup -a -A -Q 30 -d 8000 -f {sars2_fasta} {run_accession}.bam > {run_accession}.pileup")
    
    # Generate coverage file
    os.system(f"cd {STOREDIR}; cat {run_accession}.pileup | awk '{{print $2,\",\"$3,\",\"$4}}' > {run_accession}.coverage")
    
    # Index the BAM file
    os.system(f"cd {STOREDIR}; samtools index {run_accession}.bam")
    
    # Perform indel quality calibration using LoFreq
    os.system(f"cd {STOREDIR}; lofreq indelqual --dindel {run_accession}.bam -f {sars2_fasta} -o {run_accession}_fixed.bam")
    os.system(f"cd {STOREDIR}; samtools index {run_accession}_fixed.bam")
    
    # Call variants using LoFreq
    os.system(f"cd {STOREDIR}; lofreq call-parallel --no-default-filter --call-indels --pp-threads {num_threads} -f {sars2_fasta} -o {run_accession}.vcf {run_accession}_fixed.bam")
    
    # Filter variants
    os.system(f"cd {STOREDIR}; lofreq filter --af-min 0.25 -i {run_accession}.vcf -o {run_accession}_filtered.vcf")
    
    # Compress VCF files
    os.system(f"cd {STOREDIR}; bgzip -f {run_accession}.vcf")
    os.system(f"cd {STOREDIR}; bgzip -f {run_accession}_filtered.vcf")
    
    # Index compressed VCF files
    os.system(f"cd {STOREDIR}; tabix {run_accession}.vcf.gz")
    
    # Generate statistics for VCF files
    os.system(f"cd {STOREDIR}; bcftools stats {run_accession}.vcf.gz > {run_accession}.stat")
    
    # Perform SNP annotation using SnpEff
    os.system(f"cd {STOREDIR}; snpEff -q -no-downstream -no-upstream -noStats NC_063383.1 {run_accession}.vcf > {run_accession}.annot.vcf")

    
    # Generate consensus sequence from VCF
    os.system(f"{workflow_path}/bin/vcf_to_consensus.py -dp 10 -af 0.25 -v {STOREDIR}/{run_accession}.vcf.gz -d {STOREDIR}/{run_accession}.coverage -o {STOREDIR}/{run_accession}_headless_consensus.fasta -n {run_accession} -r {sars2_fasta}")
    os.system(f"{workflow_path}/bin/fix_consensus_header.py {STOREDIR}/{run_accession}_headless_consensus.fasta > {STOREDIR}/{run_accession}_consensus.fasta")
    os.system(f"cd {STOREDIR}; bgzip -f {run_accession}_consensus.fasta")
    
    # Prepare output directory and move files
    os.system(f"cd {STOREDIR}; mkdir -p {run_accession}_output")
    os.system(f"cd {STOREDIR}; mv {run_accession}_trim_summary {run_accession}.annot.vcf {run_accession}.bam {run_accession}.coverage {run_accession}.stat {run_accession}.vcf.gz {run_accession}_output")
    
    # Create output tarball
    os.system(f"cd {STOREDIR}; tar -zcvf {run_accession}_output.tar.gz {run_accession}_output")

    # Create directory and copy files
    os.system(f"mkdir -p {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")
    os.system(f"cp {workflow_path}/bin/config.yaml {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")
    os.system(f"mv {STOREDIR}/{run_accession}_output.tar.gz {STOREDIR}/{run_accession}_filtered.vcf.gz {STOREDIR}/{run_accession}_consensus.fasta.gz {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")

    # Cleanup
    os.system(f"cd {OUTDIR}; rm -rf {run_accession}*")


def ena_analysis_submit(run_accession, sample_accession, projects_accounts_csv, study_accession, workflowPath):
    # Read the line containing Webin ID and password from the projects_accounts_csv file
    line = os.popen(f"grep {study_accession} {projects_accounts_csv}").read()
    webin_id = line.split(',')[3]
    webin_password = line.split(',')[4]

    # Build paths for the output files
    output_tgz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_output.tar.gz"
    filtered_vcf_gz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_filtered.vcf.gz"
    consensus_fasta_gz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_consensus.fasta.gz"

    # Submit files using ena-analysis-submitter.py
    if study_accession == 'PRJEB55834':
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB55825 -s {sample_accession} -r {run_accession} -f {output_tgz_path} -a PATHOGEN_ANALYSIS -au {webin_id} -ap {webin_password} > {output_tgz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB55823 -s {sample_accession} -r {run_accession} -f {filtered_vcf_gz_path} -a FILTERED_VARIATION -au {webin_id} -ap {webin_password} > {filtered_vcf_gz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB55824 -s {sample_accession} -r {run_accession} -f {consensus_fasta_gz_path} -a SEQUENCE_CONSENSUS -au {webin_id} -ap {webin_password} > {consensus_fasta_gz_path}_submit.txt")
    else:
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {output_tgz_path} -a PATHOGEN_ANALYSIS -au {webin_id} -ap {webin_password} > {output_tgz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {filtered_vcf_gz_path} -a FILTERED_VARIATION -au {webin_id} -ap {webin_password} > {filtered_vcf_gz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {consensus_fasta_gz_path} -a SEQUENCE_CONSENSUS -au {webin_id} -ap {webin_password} > {consensus_fasta_gz_path}_submit.txt")


def process_data(row):
    try:
        study_accession = params['STUDY']
        run_accession = row['run_accession']
        sample_accession = row['sample_accession']
        fastq_urls = row['fastq_ftp'].split(';')
        output_folder = f"{workflow_path}/prepro/results/{study_accession}/{run_accession}_output"
        
        with open('accession_submitted.txt', 'r') as file:
            accession_submitted = file.read()

        # If output files don't exist, run map_to_reference
        if run_accession not in accession_submitted :
            map_to_reference(run_accession, sample_accession, fastq_urls[-2], fastq_urls[-1], params['SARS2_FA'], params['SARS2_FA_FAI'], params['SECRETS'], study_accession, num_threads=8, STOREDIR=params['STOREDIR'], OUTDIR=params['OUTDIR'], workflow_path=workflow_path)
            ena_analysis_submit(run_accession, sample_accession, params['SECRETS'], study_accession, workflow_path)
            os.system(f'cat {output_folder}/*txt > {output_folder}/submit.txt')
            with open(f'{output_folder}/submit.txt', 'r') as submit_file:
                submit = submit_file.read()
            if 'error' not in submit.lower():
                 with open('accession_submitted.txt', 'a') as f:
                     f.write(f'{run_accession}\n')
            os.system (f'rm -rf {output_folder} ')
    except Exception as e:
        os.system(f'echo {e} >> getion_error_illumina.err')



def main():
    index = pd.read_csv(params['INDEX'], sep='\t')
    
    # Number of tasks to run in parallel
    max_workers = 100
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_data, row) for _, row in index.iterrows()]
        
        # Wait for all tasks to finish
        for future in concurrent.futures.as_completed(futures):
            pass  # Can be extended for handling results or exceptions

if __name__ == "__main__":
    main()
