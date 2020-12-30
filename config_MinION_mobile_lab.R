#
# Copyright 2019 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

####################################################################################################
##Note: rows starting with '#' are notes for the user, and are ignored by the software
#if do_subsampling_flag <- 1, subsampling of num_fast5_files fast5 files is performed; otherwise set do_subsampling_flag <- 0
do_subsampling_flag <- 0
#num_fast5_files is the number of fast5 files to be subsampled/analysed (if do_subsampling_flag <- 1)                          
num_fast5_files <- 25
#BC_int are the barcodes used in the experiment
#BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
#barcode kits
#barcode_kits <- c("EXP-NBD103", "EXP-NBD114", "EXP-PBC001", "EXP-PBC096", "SQK-16S024", "SQK-LWB001", "SQK-RAB201", "SQK-RBK001", "SQK-RBK004", "SQK-RLB001")
barcode_kits <- c("SQK-RAB204")
#kit (1D/1D^2 reads/rapid 16S)
kit <- "SQK-RAB204"
#flowcell chemistry (R9.4/R9.5 chemistry)
flowcell <- "FLO-MIN106"
#fast_basecalling_flag <- 1 if you want to use the fast basecalling algorithm (FLO-MIN106 only); otherwise set fast_basecalling_flag <- 0 if you want to use the accurate but slow one
fast_basecalling_flag <- 1
#pair_strands_flag <- 1 if, in case a 1d2 kit and FLO-MIN107 flow-cell have been used, you want to perform 1d2 basecalling; otherwise set pair_strands_flag <- 0
pair_strands_flag <- 0
#require_two_barcodes_flag <- 1 if you want to keep only reads with a barcode (tag) at both ends of the read; otherwise set require_two_barcodes_flag <- 0
require_two_barcodes_flag <- 0
#save_space_flag <- 1 if you want temporary files to be automatically deleted; otherwise set save_space_flag <- 0
save_space_flag <- 0
#set the maximum number of threads to be used
num_threads <- 30
#set a mean amplicon length [bp]: for amplicon length I refer to the length of the biological sequence after adapters and primers trimming
amplicon_length <- 1400
#fixed_lenfil_flag <- 1 if you want to keep reads in the range [amplicon_length - lenfil_tol/2; amplicon_length + lenfil_tol/2]; otherwise set fixed_lenfil_flag <- 0 if you want to keep reads in the range [mean_length -2*sd; mean_length + 2*sd] where mean_length and sd are evaluated on a sample basis
fixed_lenfil_flag <- 1
#if fixed_lenfil_flag <- 1, lenfil_tol [bp] is the size of the window centered in amplicon_length for reads to be kept
lenfil_tol <- 300
#set primers length [bp]; if the kit you used already contains PCR primers as part of the adapters, you can set this value to 0
primers_length <- 25
#min read quality value
min_qual <- 7
#Choose taxonomic classifier between Blast and Vsearch
CLASSIFIER <- "Vsearch"
#MAX_ACCEPTS is the maximum number of hits for each query; if a value > 1 is used, a consensus taxonomy for the MAX_ACCEPTS hits is retrieved
MAX_ACCEPTS <- 3
#QUERY_COV is the minimum fraction of a query sequence that should be aligned to a sequence in the database
QUERY_COV <- 0.8
#ID_THR is the minimum alignment identity threshold
ID_THR <- 0.85
########################################################################################################
#PIPELINE DIR
PIPELINE_DIR <- "/path/to/MetONTIIME"
#MINICONDA DIR
MINICONDA_DIR <- "/path/to/miniconda3"
#basecaller_dir
BASECALLER_DIR <- "/path/to/ont-guppy-cpu/bin/"
#NCBI-downloaded sequences (QIIME2 artifact)
DB <- "/path/to/PRJNA33175_Bacterial_sequences.qza"
#Taxonomy of NCBI-downloaded sequences (QIIME2 artifact)
TAXONOMY <- "/path/to/PRJNA33175_taxonomy.qza"
#sample-metadata file describing samples metadata; it is created automatically if it doesn't exist
SAMPLE_METADATA <- "/path/to/sample-metadata.tsv"
########## End of user editable region #################################################################
#load BioStrings package
suppressMessages(library(Biostrings))
#path to MetONTIIME.sh
MetONTIIME <- paste0(PIPELINE_DIR, "/MetONTIIME.sh")
#path to subsample fast5
subsample_fast5 <- paste0(PIPELINE_DIR, "/subsample_fast5.sh")
#SEQTK
SEQTK <- paste0(MINICONDA_DIR, "/envs/MetONTIIME_env/bin/seqtk")
#PYCOQC
PYCOQC <- paste0(MINICONDA_DIR, "/envs/MetONTIIME_env/bin/pycoQC")
#NANOFILT
NANOFILT <- paste0(MINICONDA_DIR, "/envs/MetONTIIME_env/bin/NanoFilt")

