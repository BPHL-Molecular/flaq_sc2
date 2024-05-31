# FLAQ-SC2 (Florida Assembly Quality - SARS-CoV-2)
FL BPHL's SARS-CoV-2 (SC2) analysis pipeline for Illumina paired-end, whole-genome tiled-amplicon data. 

## About
FLAQ-SC2 was developed to analyze Illumina paired-end, whole-genome tiled-amplicon data (i.e., [ARTIC protocol](https://artic.network/ncov-2019)) from SC2-positive clinical specimens. The pipeline generates annotated, consensus assemblies along with reports including read/mapping/assembly quality metrics, Pango lineage, and quality flags to support screening samples prior to submissions to public repositories (e.g., GISAID, NCBI Genbank). The current version will run only on [HiPerGator](https://www.rc.ufl.edu/about/hipergator/)(HPG) using local Singularity containers for each pipeline process.

Stay tuned for FLAQ-SC2's upgrade to [Daytona](https://github.com/BPHL-Molecular/Daytona), a platform agnostic [Nextflow](https://www.nextflow.io/) workflow. Daytona is currently under active development.

## Dependencies
- Python3
- Singularity/Apptainer
- [iVar](https://github.com/andersen-lab/ivar)
- Git

To load python3 into your current environment on HiPerGator, either use `module load python` to get the lastest version of python or activate your base conda environment. For more information on how to set up your base conda environment on HPG, see the [HiPerGator Analysis Reference Guide](https://github.com/StaPH-B/southeast-region/tree/master/hipergator)).

Singularity/Apptainer will be loaded as a module during your job execution on HPG using the sbatch job script in this repository. 

The iVar Singularity container on HiPerGator does not run consistently. Instead, iVar can be installed via [conda](https://bioconda.github.io/recipes/ivar/README.html?highlight=ivar#package-package%20&#x27;ivar&#x27;). Activate your ivar conda environment with `conda activate ivar`.

Git is already installed in your HPG environment upon login.

## Primers
The default primer in the pipeline is ARTIC-V4.1.bed. If your SARS-CoV-2 data use different ARTIC primer, you need change the line "primers="4.1"" in sbatch_flaq_sc2.sh. For example, if ARTIC-V5.3.2.bed is used, primers="4.1" should be repalced with primers="5.3.2".
## Usage

For first time use, clone this repository to a directory in blue on HPG, such as in /blue/bphl-\<state\>/\<user\>/repos/bphl-molecular/.
```
cd /blue/bphl-<state>/<user>/repos/bphl-molecular/
git clone https://github.com/BPHL-Molecular/flaq_sc2.git
```
For future use, update any changes to your local repository on HPG by navigating to the flaq_sc2 repository and pulling any changes.
```
cd flaq_sc2/
git pull
```
To run the FLAQ-SC2 pipeline, copy all files from the flaq_sc2 local repository to your analysis folder. Make an input directory and copy your fastqs.
```
mkdir <analysis_dir>
cd <analysis_dir>
cp /blue/bphl-<state>/<user>/repos/bphl-molecular/flaq_sc2/* .
mkdir fastqs/
cp /path/to/fastqs/*.fastq.gz fastqs/
```
Rename your fastq files to the following format: sample_1.fastq.gz, sample_2.fastq.gz. See below for a helpful command to rename your R1 and R2 files.
```
cd fastqs/
for i in *_R1_001.fastq.gz; do mv -- "$i" "${i%[PATTERN to REMOVE]}_1.fastq.gz"; done
for i in *_R2_001.fastq.gz; do mv -- "$i" "${i%[PATTERN to REMOVE]}_2.fastq.gz"; done
```
Edit your sbatch job submission script to include your email to receive an email notification upon job END or FAIL. Replace ENTER EMAIL in `#SBATCH --mail-user=ENTER EMAIL` with your email address. Make sure there is no space between = and your email address. Edit additional sbatch parameters as needed to run your job succesfully, such as the length of time the job will run.

Submit your job.
```
sbatch sbatch_flaq_sc2.sh
```

## Main processes
- [Fastqc](https://github.com/s-andrews/FastQC)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
- [BWA](https://github.com/lh3/bwa)
- [Samtools](https://github.com/samtools/samtools)
- [iVar](https://github.com/andersen-lab/ivar)
- [Pangolin](https://github.com/cov-lineages/pangolin)
- [VADR](https://github.com/ncbi/vadr)

## Primary outputs

Outputs from each process for each individual sample can be found in a sample-specific subdirectory within the FLAQ-SC2 analysis directory. Report.txt contains the main summary report with read/mapping/assembly quality metrics, Pango lineage, and quality flags to support screening samples prior to submissions to public repositories (e.g., GISAID, NCBI Genbank). Additional details can be found in the report outputs from each process. Final assemblies (.fasta), variant files (.variant.tsv), and annotation alerts are copied into the run directory for easier access for use in downstream analyses or for samples that require manual review.

```
analysis_dir/
|__ <date>_flaq_run/
     |__ report.txt
     |__ sample1/
     |__ sample2/
|__ assemblies/
|__ vadr_error_reports/
|__ variants/
```

## Developed by:
[@SESchmedes](https://www.github.com/SESchmedes)<br />

Please email bphl16BioInformatics@flhealth.gov for questions and support.
