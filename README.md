# Running MATES

## Preparation

To run MATES, we need 3 input files, which are lists in a text format 

1. samples - Contains the sample ids that we will analyze 
2. barcodes - Contains the paths to each barcode file. Make sure to remove any additional symbol to the 16nt barcode - Only letters.
3. bams - Contains the paths to each bam file mapped to the genome. 

It is recommended that the bams contain the cell barcode tag as CB and the UMI tag as UB. We can use `sinto` for this.

Before starting, make sure of soft linking the MATES installation folder (`/gpfs/commons/groups/landau_lab/mariela/tools/MATES`) to your running directory, as well as the references `*.bed` and `*.csv` files. More on how to create these references is described in the official documentation.

Also, for running the MATES functions, you can use the environment in this path: `/gpfs/commons/home/mcorteslopez/.conda/envs/mates_env` or create a new environment by cloning it. 

### ONT

For ONT we need to modify the fastq read qualities and then re-map following by re-tagging of the bam file. This is because MATES does a quality filter and the default quality of the reads is a low quality score. The score is not meaningful per se because it is added after collapsing the read molecules, briefly, some useful commands would be: 

```bash
max_cpus=8
juncbed="/gpfs/commons/groups/landau_lab/mariela/tools/genomes/Homo_sapiens/refdata-gex-GRCh38-2020-A/genes/genes.bed"
minimapfasta="/gpfs/commons/groups/landau_lab/mariela/tools/genomes/Homo_sapiens/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
nanopore="/gpfs/commons/groups/landau_lab/mariela/tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar"
refflat="/gpfs/commons/groups/nygcfaculty/LandauONT/CH259/extra_reference_files/gencode.v31.refFlat"

# Modify fastq
fastq="${SAMPLE}.fastq"
sed 's/\$/?/g' ${fastq} > ${SAMPLE}_tmp_qc_renewed_fq.fq
fastqfile="${SAMPLE}_tmp_qc_renewed_fq.fq"

# Re-map fastq 
minimap2 -ax splice -uf -t $max_cpus --junc-bed $juncbed $minimapfasta $fastqfile | samtools view -bS -@ $max_cpus - | samtools sort -m 4G -@ $max_cpus -o bams/${SAMPLE}_passed_TE.bam - && samtools index bams/${SAMPLE}_passed_TE.bam
rm $fastqfile

## Input and output BAM files
mappingbam="bams/${SAMPLE}_passed_TE.bam"
outputbam="bams/${SAMPLE}_retag_TE.bam"

#Re-tag aligned file
sicelore="/gpfs/commons/groups/landau_lab/mariela/tools/sicelore-2.1/Jar/Sicelore-2.1.jar"
module load java 
javaXmx='-Xmx20g'
java -jar $javaXmx $sicelore AddBamMoleculeTags -I ${mappingbam} -O ${outputbam} -CELLTAG "CB" -UMITAG "UB"
samtools index ${outputbam} 

```

### PacBio

For PacBio, we can use the bam output corrected contained in the `mapped_corrected/` folder, `${sample}_me_corrected_modified.bam`. The tag conversion will only be required for the UMIs: 

```bash
sinto tagtotag --from XM --to UB --bam ${sample}_me_corrected_modified.bam -o ${sample}_me_corrected_modified_retag_full.bam
samtools index ${sample}_me_corrected_modified_retag_full.bam
```

These re-tagged bams will be listed in the bam list file. 

## Running MATES

We can execute the commands interactively or create a script to run them one by one: 


First we need to split the bams per barcode: 

```python
import MATES
from MATES import bam_processor
from MATES import data_processor
from MATES import MATES_model
from MATES import TE_quantifier


bam_processor.split_bam_files(data_mode="10X",sample_list_file='data/scRNA_samples.txt', process_num=4,
                              threads_num = 12, bam_path_file='data/scRNA_bams.txt',bc_ind='CB', long_read=True, 
                              bc_path_file='data/scRNA_barcodes.txt')

```

This part will generate a folder called `long_read` which will contain subfolders with the sample id listed in the samples file and inside of each of these folders, there will be a barcode folder for each cell barcode listed in our barcode list file. It might be recommended to remove these intermediate files after successful runs. 

After splitting the bams, we need to quantify the reads mapping to TEs: 

```python
# bam_processor.count_long_reads(TE_mode="exclusive", data_mode="10X", threads_num=6, <- Use this if keeping only intragenic TEs
bam_processor.count_long_reads(TE_mode="inclusive", data_mode="10X", threads_num=6, 
                               sample_list_file='data/scRNA_samples.txt',
                               bc_path_file='data/scRNA_barcodes.txt')
```

This process takes a long time to run and might require additional resources, depending of the sample. If successful, you will see a `count_long_reads` folder with individual barcode named folders.
The  last step allow us to quantify the TEs and generate the count matrices that we can use in downstream tools. 

```python
from MATES.scripts import TE_locus_quantifier
#TE_locus_quantifier.unique_locus_TE_MTX(TE_mode="exclusive", data_mode="10X", <- Use this if keeping only intragenic TEs
TE_locus_quantifier.unique_locus_TE_MTX(TE_mode="inclusive", data_mode="10X", 
                    sample_list_file='data/scRNA_samples.txt', long_read = True)
```

This part might use the most memory, because it has to merge all the counts. It generates a folder  called `long_read_output/Unique_exclusive` or `long_read_output/Unique_inclusive` (depending of your settings) where subsets with the sample names contain the matrix, barcode and feature files used for downstream analysis. 

I include a script if you prefer to run these steps in the cluster at [scripts/run_mates.sh](scripts/run_mates.sh).


More MATES documentation in [their github](https://github.com/mcgilldinglab/MATES)