#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=flaq_sc2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ENTER EMAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=3-00
#SBATCH --output=flaq_sc2.%j.out
#SBATCH --error=flaq_sc2.%j.err

#Run script/command and use $SLURM_CPUS_ON_NODE
module load singularity

primer_version="4.1"

python flaq_sc2.py fastqs/ --primer_bed primers/ARTIC-V${primer_version}.bed --lib_frag frag --threads $SLURM_CPUS_ON_NODE --ref_fasta reference/nCoV-2019.reference.fasta --ref_gff reference/GCF_009858895.2_ASM985889v3_genomic.gff --sotc S:L452R,S:E484K  --pango_path /apps/staphb-toolkit/containers/pangolin_4.0.5-pdata-1.3.sif --pangolin v4.0.5 --pangolin_data v1.3 
