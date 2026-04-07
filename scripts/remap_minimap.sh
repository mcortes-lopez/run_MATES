#!/bin/bash
#SBATCH --mail-user=mcorteslopez@nygenome.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=mapping_multimap
#SBATCH --mem=48G
#SBATCH --partition=pe2
#SBATCH --output=logs/%x_%a.log
#SBATCH --error=logs/%x_%a.e
#SBATCH --cpus-per-task=12
#SBATCH --array=0-11

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate sicelore2.1
module load samtools

max_cpus=8
juncbed="/gpfs/commons/groups/landau_lab/mariela/tools/genomes/Homo_sapiens/refdata-gex-GRCh38-2020-A/genes/genes.bed"
minimapfasta="/gpfs/commons/groups/landau_lab/mariela/tools/genomes/Homo_sapiens/refdata-gex-GRCh38-2020-A/fasta/genome.fa"

declare -a SAMPLES=($(find /gpfs/commons/groups/landau_lab/mariela/CRC_project/output/oarfish/fastq_demux -type f | xargs -I{} basename {} ".fastq" ))

SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}

fastq="/gpfs/commons/groups/landau_lab/mariela/CRC_project/output/oarfish/fastq_demux/${SAMPLE}.fastq"
sed 's/\$/?/g' ${fastq} > ${SAMPLE}_tmp_qc_renewed_fq.fq
fastqfile="${SAMPLE}_tmp_qc_renewed_fq.fq"

minimap2 -ax splice -uf -t $max_cpus --junc-bed $juncbed $minimapfasta $fastqfile | samtools view -bS -@ $max_cpus - | samtools sort -m 4G -@ $max_cpus -o bams/${SAMPLE}_passed_TE.bam - && samtools index bams/${SAMPLE}_passed_TE.bam
rm $fastqfile
