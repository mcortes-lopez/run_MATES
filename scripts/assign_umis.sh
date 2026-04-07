#!/bin/bash
#SBATCH --mail-user=mcorteslopez@nygenome.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=umiassign
#SBATCH --partition=pe2
#SBATCH --output=logs/%x_%a.log
#SBATCH --error=logs/%x_%a.e
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --array=0-11

module purge
module load java 

export PATH="/nfs/sw/java/jdk-17.0.4/bin/:$PATH"
export JAVA_CMD="/nfs/sw/java/jdk-17.0.4/bin/java"

nanopore="/gpfs/commons/groups/landau_lab/mariela/tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar"
refflat="/gpfs/commons/groups/nygcfaculty/LandauONT/CH259/extra_reference_files/gencode.v31.refFlat"
max_cpus=4


declare -a SAMPLES=($(find /gpfs/commons/groups/landau_lab/mariela/CRC_project/output/oarfish/fastq_demux -type f | xargs -I{} basename {} ".fastq" ))
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}

# Input and output BAM files
mappingbam="bams/${SAMPLE}_passed_TE.bam"
outputbam="bams/${SAMPLE}_retag_TE.bam"

# Run Java command
sicelore="/gpfs/commons/groups/landau_lab/mariela/tools/sicelore-2.1/Jar/Sicelore-2.1.jar"
module load java 
javaXmx='-Xmx20g'
#java -jar $javaXmx $sicelore AddBamMoleculeTags -I ${mappingbam} -O ${outputbam} -CELLTAG "CB" -UMITAG "UB"
#samtools index -@4 ${outputbam} 

if java -jar $javaXmx $sicelore AddBamMoleculeTags -I ${mappingbam} -O ${outputbam} -CELLTAG "CB" -UMITAG "UB"; then
    echo "Java command completed successfully."
    
    # Index the output BAM
    if samtools index -@4 ${outputbam}; then
        echo "Indexing completed successfully."

        # Safely remove the input file
        echo "Removing input BAM file: ${mappingbam}"
        rm -f ${mappingbam}
    else
        echo "Error: Failed to index the output BAM. Input BAM will not be removed."
        exit 1
    fi
else
    echo "Error: Java command failed. Input BAM will not be removed."
    exit 1
fi