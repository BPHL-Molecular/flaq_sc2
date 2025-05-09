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
module load apptainer

primer_version="4.1"

#Pangolin-related arguments/parameters are now optional. By default, the pipeline will now pull the latest pangolin and pangolin-data versions, as well as nextclade.
python flaq_sc2.py fastqs/ --primer_bed primers/ARTIC-V${primer_version}.bed --lib_frag frag --threads $SLURM_CPUS_ON_NODE --ref_fasta reference/nCoV-2019.reference.fasta --ref_gff reference/GCF_009858895.2_ASM985889v3_genomic.gff --sotc S:L452R,S:E484K
