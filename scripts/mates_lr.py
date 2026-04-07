#!/usr/bin/env python

import MATES
from MATES import bam_processor
from MATES import data_processor
from MATES import MATES_model
from MATES import TE_quantifier

# Commented out IPython magic to ensure Python compatibility.
# %cd /gpfs/commons/groups/landau_lab/mariela/SON_project/MATES/scripts/

#bam_processor.split_bam_files(data_mode="10X",sample_list_file='data/scRNA_samples.txt', process_num=4,
#                              threads_num = 12, bam_path_file='data/scRNA_bams.txt',bc_ind='CB', long_read=True, 
#                              bc_path_file='data/scRNA_barcodes.txt')



#bam_processor.count_long_reads(TE_mode="exclusive", data_mode="10X", threads_num=6, 
bam_processor.count_long_reads(TE_mode="inclusive", data_mode="10X", threads_num=6, 
                               sample_list_file='data/scRNA_samples.txt',
                               bc_path_file='data/scRNA_barcodes.txt')
# Parameters
### count_long_reads(TE_mode, data_mode, threads_num, sample_list_file, ref_path='Default', bc_path_file=None):

from MATES.scripts import TE_locus_quantifier
#TE_locus_quantifier.unique_locus_TE_MTX(TE_mode="exclusive", data_mode="10X", 
TE_locus_quantifier.unique_locus_TE_MTX(TE_mode="inclusive", data_mode="10X", 
                    sample_list_file='data/scRNA_samples.txt', long_read = True)
# Parameters
#
#unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False)

#from MATES import TE_locus_quantifier
#TE_locus_quantifier.unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False, bc_path_file=None)
