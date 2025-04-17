#!/usr/bin/env python

#Author: Sarah Schmedes
#Email: sarah.schmedes@flhealth.gov
#Uses local containers, VADR for error flags, reports pangolin output in main report
'''
This program takes in Illumina paired-end fastqs using the ARTIC primer schemes for SARS-CoV-2
and generates a consensus assembly by using a reference-based assembly method by mapping read to
a reference.
'''

import os
import sys
import subprocess
import argparse
import datetime
import pandas as pd
import re
import os.path

#Parse arguments, get path for fastqs, primer version
parser = argparse.ArgumentParser(usage='flaq_sc2.py <input_dir> [options]')
parser.add_argument('input', help='path to dir with fastqs')
parser.add_argument('--primer_bed', help='path to ARTIC SC2 primer bed', required=True)
parser.add_argument('--lib_frag', default='frag', choices=['no_frag', 'frag'], help='specify if input amplicons were fragmented, (default: no_frag)', required=True)
parser.add_argument('--threads', default=8, dest='threads', help='specify number of threads, (default: %(default)s)', required=True)
parser.add_argument('--ref_fasta', help='path to reference fasta', required=True)
parser.add_argument('--ref_gff', help='path to reference gff', required=True)
parser.add_argument('--sotc', help='comma separated list of SOTCs to screen (e.g., S:L452R,S:E484K', required=True)
parser.add_argument('--pango_path', help='path to pangolin container', required=False)
parser.add_argument('--pangolin', help='pangolin version (e.g., v2.3)', required=False)
parser.add_argument('--pangolin_data', help='pangolin-data version (e.g., v1.3', required=False)
parser.add_argument('--version', action='version', version='This is flaq_sc2: Version 2.6.7', help='print version')

if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

input_dir = os.path.abspath(args.input) + '/'
primers = os.path.abspath(args.primer_bed)
frag = args.lib_frag
threads = str(args.threads)
ref = os.path.abspath(args.ref_fasta)
gff = os.path.abspath(args.ref_gff)
sotc_v = args.sotc
sotc_v = sotc_v.split(',')
pango_path = args.pango_path
pango_v = args.pangolin
pdata = args.pangolin_data
if pango_path is not None and pango_v is not None and pdata is not None:
    pangolin = pango_v + '_pdata-' + pdata
else:
    pangolin = "auto-pulled"

cwd = os.getcwd() + '/'

output_dir = cwd + datetime.date.today().strftime('%Y-%m-%d') + '_flaq_run'
subprocess.run('mkdir -p ' + output_dir, shell=True, check=True) #make output directory date_flaq_legion_run
subprocess.run('mkdir -p assemblies', shell=True, check=True) #make final assemblies folder
subprocess.run('mkdir -p variants', shell=True, check=True) #make final variants folder
subprocess.run('mkdir -p vadr_error_reports', shell=True, check=True) #folder for vadr error report for easier review
if os.path.exists('nextclade_latest.sif'):
    subprocess.run('rm nextclade_latest.sif', shell=True, check=True)
subprocess.run('singularity pull nextclade_latest.sif docker://nextstrain/nextclade:latest', shell=True, check=True) #pull latest nextclade container
subprocess.run('mkdir -p data/nextclade', shell=True, check=True) #create nextclade dataset folder
subprocess.run('singularity exec nextclade_latest.sif nextclade dataset get --name sars-cov-2 --output-dir data/nextclade', shell=True, check=True) #get latest nextclade dataset

# Pull the latest pangolin container if path not provided
if not args.pango_path:
    if os.path.exists('pangolin_latest.sif'):
        subprocess.run('rm pangolin_latest.sif', shell=True, check=True)
    subprocess.run('singularity pull pangolin_latest.sif docker://staphb/pangolin:latest', shell=True, check=True)
    # Get the container version info
    try:
        proc = subprocess.run('singularity exec pangolin_latest.sif pangolin --all-versions', shell=True, capture_output=True, text=True, check=True)
        version_info = proc.stdout.strip()
        # Extract pangolin version
        pango_match = re.search(r'pangolin: ([\d.]+)', version_info)
        # Extract pangolin-data version
        pdata_match = re.search(r'pangolin-data: ([\d.]+)', version_info)
        if pango_match and pdata_match:
            pango_v = "v" + pango_match.group(1)
            pdata = pdata_match.group(1)
            pangolin = pango_v + '_pdata-' + pdata
        else:
            pango_v = "latest"
            pdata = "latest"
            pangolin = "latest_auto-pulled"
    except Exception as e:
        print(f"Warning: Error extracting pangolin version: {e}")
        pango_v = "latest"
        pdata = "latest"
        pangolin = "latest_auto-pulled"
    pango_path = "pangolin_latest.sif"
else:
    if pango_v is not None and pdata is not None:
        pangolin = pango_v + '_pdata-' + pdata
    else:
        pangolin = "user-specified-container"

adapters = '/bbmap/resources/adapters.fa'
phix = '/bbmap/resources/phix174_ill.ref.fa.gz'

#Get sample names
samples = []
fastqs = []

#Look at some code examples to get fastq names R1, _1 or R1_001 (make work for more sample types later)
for f in os.listdir(input_dir):
    if f.endswith('.fastq.gz'):
        fastqs.append(f)
        sn = f.split("_")
        sn = sn[0]
        samples.append(sn)
unique = set(samples)
samples = list(unique)
samples.sort()

#Create output file
report = open(output_dir + '/report.txt', 'w')
header = ['sampleID', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag', 'pangolin_version', 'lineage', 'SOTC']
report.write('\t'.join(map(str,header)) + '\n')

#Run pipeline for each sample
for s in samples:
    sample_dir = output_dir + '/' + s + '/'
    subprocess.run('mkdir -p ' + sample_dir, shell=True, check=True) #mkdir for each sample name
    subprocess.run('cp ' + input_dir + s + '*.fastq.gz ' + sample_dir, shell=True, check=True) #cp fastqs to new dir

    out_log = open(sample_dir + s + '.out', 'w')
    err_log = open(sample_dir + s + '.err', 'w')

    #Get number of raw reads
    proc_1 = subprocess.run('zcat ' + sample_dir + s + '_1.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_1.stdout.rstrip()
    reads_1 = int(wc_out_1) / 4
    proc_2 = subprocess.run('zcat ' + sample_dir + s + '_2.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_2 = proc_2.stdout.rstrip()
    reads_2 = int(wc_out_2) / 4
    raw_reads = reads_1 + reads_2
    raw_reads = int(raw_reads)

    #Run fastqc on original reads
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fastq.gz --threads ' + threads + ' --outdir ' + sample_dir, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_original_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_original_fastqc.zip', shell=True, check=True)

    #Run trimmomatic
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/trimmomatic_0.39.sif trimmomatic PE -threads ' + threads + ' -phred33 -trimlog ' + sample_dir + s + '.log ' + sample_dir + s + '_1.fastq.gz ' + sample_dir + s + '_2.fastq.gz ' + sample_dir + s + '_1.trim.fq.gz ' + sample_dir + s + '_unpaired_1.trim.fq.gz ' + sample_dir + s + '_2.trim.fq.gz ' + sample_dir + s + '_unpaired_2.trim.fq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ' + sample_dir + s + '_trimstats.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    #rm unpaired reads
    subprocess.run('rm ' + sample_dir + s + '*_unpaired*.trim.fq.gz', shell=True, check=True)
    #rm fastq files copied from previous dir
    subprocess.run('rm ' + sample_dir + s + '*.fastq.gz', shell=True, check=True)

    #Run bbduk to remove Illumina adapter sequences and any PhiX contamination
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=' + sample_dir + s + '_1.trim.fq.gz in2=' + sample_dir + s + '_2.trim.fq.gz out1=' + sample_dir + s + '_1.rmadpt.fq.gz out2=' +sample_dir + s + '_2.rmadpt.fq.gz ref=' + adapters + ' stats=' + sample_dir + s + '.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=' + sample_dir + s + '_1.rmadpt.fq.gz in2=' + sample_dir + s + '_2.rmadpt.fq.gz out1=' + sample_dir + s + '_1.fq.gz out2=' + sample_dir + s + '_2.fq.gz outm=' + sample_dir + s + '_matchedphix.fq ref=' + phix + ' k=31 hdist=1 stats=' + sample_dir + s + '_phixstats.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('rm ' + sample_dir + '*.trim.fq.gz', shell=True, check=True)
    subprocess.run('rm ' + sample_dir + '*.rmadpt.fq.gz', shell=True, check=True)

    #Run fastqc on clean forward and reverse reads
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fq.gz --threads ' + threads, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_clean_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_clean_fastqc.zip', shell=True, check=True)

    #Run multiqc
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/multiqc_1.8.sif multiqc ' + sample_dir + '*_fastqc.zip -o ' + sample_dir, shell=True, stdout=out_log, stderr=err_log, check=True)

    #Get number of clean reads
    proc_c1 = subprocess.run('zcat ' + sample_dir + s + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c1 = proc_c1.stdout.rstrip()
    reads_c1 = int(wc_out_c1) / 4
    proc_c2 = subprocess.run('zcat ' + sample_dir + s + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c2 = proc_c2.stdout.rstrip()
    reads_c2 = int(wc_out_c2) / 4
    clean_reads = reads_c1 + reads_c2
    clean_reads = int(clean_reads)

    #Map reads to reference
    align_dir = sample_dir + 'alignment/'
    subprocess.run('mkdir ' + align_dir, shell=True, check=True)

    #If frag == 'no_frag', do not remove PCR duplicates
    if frag == 'no_frag':
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bwa_0.7.17.sif bwa mem -t ' + threads + ' ' + ref + ' ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools view - -F 4 -u -h | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort > ' + align_dir + s + '.sorted.bam', shell=True, check=True)

    #If frag == 'frag', remove PCR duplicates
    elif frag == 'frag':
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/bwa_0.7.17.sif bwa mem -t ' + threads + ' ' + ref + ' ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools view - -F 4 -u -h | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -n > ' + align_dir + s + '.namesorted.bam', shell=True, check=True) #output name sorted bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools fixmate -m ' + align_dir + s + '.namesorted.bam ' + align_dir + s + '.fixmate.bam', shell=True, check=True) #fixmate
        #Create positional sorted bam from fixmate.bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -o ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.fixmate.bam', shell=True, check=True)
        #Mark duplicate reads
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools markdup ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.markdup.bam', shell=True, check=True)
        #Remove duplicate reads
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools markdup -r ' + align_dir + s + '.positionsort.bam ' + align_dir + s + '.dedup.bam', shell=True, check=True)
        #Sort dedup.bam and rename to .sorted.bam
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort -o ' + align_dir + s + '.sorted.bam ' + align_dir + s + '.dedup.bam', shell=True, check=True)

    #Index final sorted bam from either no_frag or frag paths
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools index ' + align_dir + s + '.sorted.bam', shell=True, check=True)

    #Trim primers with iVar
    subprocess.run('cd ' + align_dir + ' && ivar trim -i ' + s + '.sorted.bam  -b ' + primers + ' -p ' + s + '.primertrim -e', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort ' + align_dir + s + '.primertrim.bam -o ' + align_dir + s + '.primertrim.sorted.bam', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools index ' + align_dir + s + '.primertrim.sorted.bam', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools coverage ' + align_dir  + s + '.primertrim.sorted.bam -o ' + align_dir + s + '.coverage.txt', shell=True, stdout=out_log, stderr=err_log, check=True)

    #Get map stats
    with open(align_dir + s + '.coverage.txt', 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        mapq = stats[8]

    #Get percentage of mapped reads/clean reads
    percent_map = "%0.4f"%(((int(reads_mapped)/int(clean_reads)))*100)

    #Call variants
    var_path = sample_dir + 'variants/'
    subprocess.run('mkdir ' + var_path, shell=True, check=True)
    subprocess.run('cd ' + var_path + ' && samtools mpileup -A -d 8000 --reference ' + ref + ' -B -Q 0 ' + align_dir + s + '.primertrim.sorted.bam | ivar variants -r ' + ref + ' -m 10 -p ' + s + '.variants -q 20 -t 0.25 -g ' + gff, shell=True, stdout=out_log, stderr=err_log, check=True)

    #Generate consensus assembly
    assem_path = sample_dir + 'assembly/'
    subprocess.run('mkdir ' + assem_path, shell=True, check=True)
    subprocess.run('cd ' + assem_path + ' && samtools mpileup -A -B -d 8000 --reference ' + ref + ' -Q 0 ' + align_dir + s + '.primertrim.sorted.bam | ivar consensus -t 0 -m 10 -n N -p ' + s + '.consensus', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Gather QC metrics for consensus assembly
    with open(assem_path + s + '.consensus.fa', 'r') as assem:
        header = assem.readline()
        header = header.rstrip()
        bases = assem.readline()
        bases = bases.rstrip()
        num_bases = len(bases)
        ns = bases.count('N')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)

    #Rename header in fasta to just sample name
    subprocess.run('sed -i \'s/^>.*/>' + s + '/\' ' + assem_path + s + '.consensus.fa', shell=True, check=True)

    #QC flag
    pg_flag = ''
    dp_flag = ''
    qc_flag = ''
    if float(pg) < 79.5:
        pg_flag = 'FAIL: Percent genome < 80%'
        qc_flag = qc_flag + pg_flag
    else:
        if float(depth) < 100:
            dp_flag = 'FAIL: Mean read depth < 100x'
            qc_flag = qc_flag + dp_flag
        if qc_flag == '':
            qc_flag = qc_flag + 'PASS'

    #Copy passing ivar assemblies to assemblies folder and variant files to variants folder
    if qc_flag == 'PASS':
        subprocess.run('cp ' + assem_path + s + '.consensus.fa assemblies/', shell=True, check=True)
        subprocess.run('cp ' + var_path + s + '.variants.tsv variants/', shell=True, check=True)

        #Run VADR
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/vadr_1.3.sif /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 ' + assem_path + s + '.consensus.fa > ' + assem_path + s + '.trimmed.fasta', shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/vadr_1.3.sif v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ ' + assem_path + s + '.trimmed.fasta -f ' + assem_path + 'vadr_results', shell=True, stdout=out_log, stderr=err_log, check=True)

        #Parse through VADR outputs to get PASS or REVIEW flag
        vadr_flag = ''

        with open(assem_path + 'vadr_results/vadr_results.vadr.pass.list', 'r') as p_list:
            result = p_list.readline()
            result = result.rstrip()
            if result == s:
                vadr_flag = 'PASS'

        with open(assem_path + 'vadr_results/vadr_results.vadr.fail.list', 'r') as f_list:
            result = f_list.readline()
            result = result.rstrip()
            if result == s:
                vadr_flag = 'REVIEW'

        #Copy VADR error report to main analysis folder for easier review
        if vadr_flag == 'REVIEW':
            subprocess.run('cp ' + assem_path + 'vadr_results/vadr_results.vadr.alt.list vadr_error_reports/', shell=True, check=True)
            subprocess.run('mv vadr_error_reports/vadr_results.vadr.alt.list vadr_error_reports/' + s + '.vadr.alt.list', shell=True, check=True)

        #Run pangolin
        subprocess.run('singularity exec -B $(pwd):/data ' + pango_path + ' pangolin ' + assem_path + s + '.consensus.fa', shell=True, check=True)
        #Get lineage
        proc = subprocess.run('tail -n 1 lineage_report.csv | cut -d \',\' -f 2', shell=True, check=True, capture_output=True, text=True)
        lineage = proc.stdout.rstrip()
        subprocess.run('mv lineage_report.csv ' + assem_path, shell=True, check=True)

        #Run nextclade
        subprocess.run('singularity exec -B $(pwd):/data nextclade_latest.sif nextclade run --output-csv ' + assem_path + 'nextclade_report.csv --jobs ' + threads + ' --input-dataset data/nextclade ' + assem_path + s + '.consensus.fa', shell=True, check=True)

        #Parse nextclade output and screen for sotc
        with open(assem_path + 'nextclade_report.csv', 'r') as nc:
            header = nc.readline()
            c_results = nc.readline()
            c_results = c_results.rstrip()
            data = c_results.split(',')
            sotc = []
            for v in sotc_v:
                if v in data:
                    sotc.append(v)
            sotc_out = (',').join(sotc)
    else:
        vadr_flag = 'NA'
        lineage = 'NA'
        sotc_out = 'NA'

    #Write to output file
    results = [s, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag, pangolin, lineage, sotc_out]
    report.write('\t'.join(map(str,results)) + '\n')
    out_log.close()
    err_log.close()

report.close()

#Create multi-fasta file of only passing assemblise
subprocess.run('cat assemblies/*.fa > assemblies.fasta', shell=True, check=True)

# Get nextclade version
try:
    proc = subprocess.run('singularity exec nextclade_latest.sif nextclade --version', shell=True, capture_output=True, text=True, check=True)
    nextclade_version_output = proc.stdout.strip()
    version_match = re.search(r'nextclade\s+([\d.]+)', nextclade_version_output)
    if version_match:
        nextclade_version = version_match.group(1)
    else:
        nextclade_version = "latest"
except:
    nextclade_version = "latest"

#Run nextclade
subprocess.run('singularity exec -B $(pwd):/data nextclade_latest.sif nextclade run --output-tsv nextclade_report.tsv --jobs ' + threads + ' --input-dataset data/nextclade assemblies.fasta', shell=True, check=True)

if os.path.exists('nextclade_report.tsv'):
    df = pd.read_csv('nextclade_report.tsv', sep='\t')
    df.insert(df.columns.get_loc("clade"), "nextclade_version", nextclade_version)
    df.to_csv('nextclade_report.tsv', sep='\t', index=False)
