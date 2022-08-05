#
# Copyright 2022 Simone Maestri. All rights reserved.
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
#conf_basecalling_flag <- 1 if you want to specify a configuration file for base-calling (and additional parameters) insted of choosing the default by specifying kit and flowcell
conf_basecalling_flag <- 0
#conf_par_basecalling is the name of the config file (and additional parameters, such as the device for GPU-accelerated basecalling) in case config_basecalling_flag <- 1
conf_par_basecalling <- "-c dna_r10.4_e8.1_hac.cfg --device 'auto' --gpu_runners_per_device 1"
#if skip_demultiplexing_flag <- 1 demultiplexing is skipped; otherwise set skip_demultiplexing_flag <- 0
skip_demultiplexing_flag <- 0
#require_two_barcodes_flag <- 1 if you want to keep only reads with a barcode (tag) at both ends of the read; otherwise set require_two_barcodes_flag <- 0
require_two_barcodes_flag <- 0
#save_space_flag <- 1 if you want temporary files to be automatically deleted; otherwise set save_space_flag <- 0
save_space_flag <- 0
#set the maximum number of threads to be used
num_threads <- 30
#min_seq_length is the minimum sequence length (bp) to be retained
min_seq_length <- 1200
#max_seq_length is the maximum sequence length (bp) to be retained
max_seq_length <- 1600
#set primers length [bp]; if the kit you used already contains PCR primers as part of the adapters, you can set this value to 0
primers_length <- 25
#min read quality value
min_qual <- 9
#Choose taxonomic classifier between Blast and Vsearch
CLASSIFIER <- "Vsearch"
#MAX_ACCEPTS is the maximum number of hits for each query; if a value > 1 is used, a consensus taxonomy for the MAX_ACCEPTS hits is retrieved
MAX_ACCEPTS <- 3
#QUERY_COV is the minimum fraction of a query sequence that should be aligned to a sequence in the database
QUERY_COV <- 0.8
#ID_THR is the minimum alignment identity threshold
ID_THR <- 0.90
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

