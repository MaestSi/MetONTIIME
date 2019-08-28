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

args = commandArgs(trailingOnly=TRUE)

config_file_ind <- grep(x = args, pattern = "\\.R")
raw_reads_ind <- setdiff(c(1, 2), config_file_ind)

config_file <- args[config_file_ind]
d1_tmp <- args[raw_reads_ind]

#load variables
source(config_file)

if (!exists("do_subsampling_flag")) {
  do_subsampling_flag <- 0
}

if (!exists("num_fast5_files")) {
  num_fast5_files <- 100
}

if (!exists("fast_basecalling_flag") || flowcell == "FLO-MIN107") {
  fast_basecalling_flag <- 0
}

if (!exists("amplicon_length")) {
  amplicon_length <- 710
}

if (!exists("fixed_lenfil_flag")) {
  fixed_lenfil_flag <- 0
}

if (!exists("lenfil_tol")) {
  lenfil_tol <- 300
}

if (!exists("pair_strands_flag") || flowcell == "FLO-MIN106") {
  pair_strands_flag <- 0
}

if (do_subsampling_flag == 1) {
    #d2 is the directory which is going to include processed reads
    d2 <- paste0(dirname(d1_tmp), "/", basename(d1_tmp), "_", num_fast5_files, "_subsampled_fast5_files_analysis")
} else {
    d2 <- paste0(dirname(d1_tmp), "/", basename(d1_tmp), "_analysis")
}

d2_basecalling <- paste0(d2, "/basecalling")
d2_preprocessing <- paste0(d2, "/preprocessing")

#d3 is the directory which is going to include results
d3 <- paste0(d2, "/analysis")

#logfile is the file which is going to include the log
logfile <- paste0(d3, "/logfile.txt")     

if (!exists("save_space_flag")) {
  save_space_flag <- 0
}

if (pair_strands_flag == 1) {
  basecaller <- paste0(BASECALLER_DIR, "/guppy_basecaller")
  basecaller_1d2 <- paste0(BASECALLER_DIR, "/guppy_basecaller_1d2")
} else {
  basecaller <- paste0(BASECALLER_DIR, "/guppy_basecaller")
}

demultiplexer <- paste0(BASECALLER_DIR, "/guppy_barcoder")

basecaller_version <- system(paste0(basecaller, " --version"), intern = TRUE)

if (!dir.exists(d2)) {
  dir.create(d2)
  dir.create(d2_basecalling)
  dir.create(d2_preprocessing)
  dir.create(d3)
  cat(text = "Welcome to MinION meta-barcoding mobile laboratory!", file = logfile, sep = "\n", append = TRUE)
  cat(text = "Welcome to MinION meta-barcoding mobile laboratory!", sep = "\n")
  cat(text = "\n", file = logfile, append = TRUE)
  cat(text = "\n")
  
  if (do_subsampling_flag == 1) {
    cat(text = paste0("Subsampling ", num_fast5_files, " fast5 files from directory ", d1_tmp), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Subsampling ", num_fast5_files, " fast5 files from directory ", d1_tmp), sep = "\n")
    system(command = paste0(subsample_fast5, " ", d1_tmp, " ", num_fast5_files))
    d1 <- paste0(gsub("\\/$", "", d1_tmp), "_", num_fast5_files, "_subsampled_fast5_files")
  } else {
    d1 <- d1_tmp
  }   

  cat(text = paste0("Raw reads directory: ", d1), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Raw reads directory: ", d1), sep = "\n")  
  cat(text = paste0("Basecalled reads directory: ", d2_basecalling), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Basecalled reads directory: ", d2_basecalling), sep = "\n")
  cat(text = paste0("Preprocessing directory: ", d2_preprocessing), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Preprocessing directory: ", d2_preprocessing), sep = "\n")
  cat(text = paste0("Analysis and results directory: ", d3), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Analysis and results directory: ", d3), sep = "\n")
  cat(text = paste0("Flow-cell: ", flowcell), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Flow-cell: ", flowcell), sep = "\n")
  cat(text = paste0("Kit: ", kit), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Kit: ", kit), sep = "\n")
  if (exists("amplicon_length")) {
    cat(text = paste0("Expected amplicon length [bp]: ", amplicon_length), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Expected amplicon length [bp]: ", amplicon_length), sep = "\n")
  }
  cat(text = paste0("Barcodes used in this experiment: ", paste0(BC_int, collapse = ", ")), file = logfile, sep = ", ", append = TRUE)
  cat(text = paste0("Barcodes used in this experiment: ", paste0(BC_int, collapse = ", ")), sep = ", ")
  cat(text = "\n", file = logfile, append = TRUE)
  cat(text = "\n")
  cat(text = paste0("Basecalling is going to be performed by ", basecaller_version), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Basecalling is going to be performed by ", basecaller_version), sep = "\n")
  if (fast_basecalling_flag == 1 && flowcell == "FLO-MIN106") {
    cat(text = "Basecalling model: fast", file = logfile, sep = "\n", append = TRUE)
    cat(text = "Basecalling model: fast", sep = "\n")
  } else if (fast_basecalling_flag != 1 && flowcell == "FLO-MIN106") {
    cat(text = "Basecalling model: high-accuracy", file = logfile, sep = "\n", append = TRUE)
    cat(text = "Basecalling model: high-accuracy", sep = "\n")
  }
  cat(text = paste0("Demultiplexing is going to be performed by guppy_barcoder after basecalling"), file = logfile, sep = ", ", append = TRUE)
  cat(text = paste0("Demultiplexing is going to be performed by guppy_barcoder after basecalling"), sep = ", ")
  cat(text = "\n", file = logfile, append = TRUE)
  cat(text = "\n")
} else {
  cat(text = paste0(d2, " directory already exists; delete or rename it and then restart the analysis!"), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0(d2, " directory already exists; delete or rename it and then restart the analysis!"), sep = "\n")
  return()
}

BC_tot <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_tot_full <- c(BC_tot, paste0("BC", 13:99))
BC_trash <- setdiff(BC_tot_full, BC_int)
BC_trash <- paste0("^", BC_trash)

cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

cat(text = paste0("Basecalling started at ", date()), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Basecalling started at ", date()), sep = "\n")

num_threads_caller <- round(num_threads/4)
if (fast_basecalling_flag == 1) {
  system(paste0(basecaller, " -r -i ", d1, " --cpu_threads_per_caller ", num_threads_caller, " --num_callers 4", " -c dna_r9.4.1_450bps_fast.cfg --hp_correct TRUE --fast5_out -s ", d2_basecalling, " --disable_pings"))
} else {
  system(paste0(basecaller, " -r -i ", d1, " --cpu_threads_per_caller ", num_threads_caller, " --num_callers 4", " --flowcell ", flowcell, " --kit ", kit, " --hp_correct TRUE --fast5_out -s ", d2_basecalling, " --disable_pings"))
}

if (pair_strands_flag == 1) {
  system(paste0(basecaller_1d2, " -r -i ", d2_basecalling, "/workspace --cpu_threads_per_caller ", num_threads_caller, " --num_callers 4", " --config dna_r9.5_450bps_1d2_raw.cfg -f ", d2_basecalling, "/sequencing_summary.txt -s ", d2, "/basecalling_1d2 --disable_pings"))
  d2_basecalling <- paste0(d2, "/basecalling_1d2")
}

cat(text = paste0("Basecalling finished at ", date()), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Basecalling finished at ", date()), sep = "\n")
cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

cat(text = paste0("Demultiplexing started at ", date()), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Demultiplexing started at ", date()), sep = "\n")
system(paste0(demultiplexer, " -r -i ", d2_basecalling, " -t ", num_threads, " -s ", d2_preprocessing, " --barcode_kits \"", paste0(barcode_kits, collapse = " "), "\"", " --kit ", kit))
cat(text = paste0("Demultiplexing finished at ", date()), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Demultiplexing finished at ", date()), sep = "\n")
cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

fastq_files <- list.files(path = d2_preprocessing, pattern = ".*\\.fastq", full.names = TRUE, recursive = TRUE)
num_reads_tot <- 0
for (i in 1:length(fastq_files)) {
  num_reads_tot <- num_reads_tot + length(grep(x = readLines(fastq_files[i]), pattern = "^\\+$"))
}

barcode_dirs <- grep(x = list.dirs(d2_preprocessing), pattern = paste0("barcode", substr(x = BC_int, start = 3, stop = 5), collapse = "|"), value = TRUE)
for (i in 1:length(barcode_dirs)) {
  fastq_files_curr_barcode <- list.files(path = barcode_dirs[i], pattern = ".*\\.fastq", full.names = TRUE)
  system(paste0("cat ", paste0(fastq_files_curr_barcode, collapse = " "), " > ", d2_preprocessing, "/BC", substr(x = basename(barcode_dirs[i]), start = 8, stop = 9), "_tmp1.fastq"))
}

cat(text = paste0("Number of basecalled reads: ", num_reads_tot), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Number of basecalled reads: ", num_reads_tot), sep = "\n")  
cat(text = "", file = logfile, sep = "\n", append = TRUE)
cat(text = "", sep = "\n")
cat(text = paste0("Creating folder: ", d2, "/qc, which is going to include stats about the sequencing run"), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Creating folder: ", d2, "/qc, which is going to include stats about the sequencing run"), sep = "\n")
cat(text = "", file = logfile, sep = "\n", append = TRUE)
cat(text = "", sep = "\n")
dir.create(paste0(d2, "/qc"))
cat(text = "Now performing quality control with PycoQC", file = logfile, sep = "\n", append = TRUE)
cat(text = "Now performing quality control with PycoQC", sep = "\n")
cat(text = "", file = logfile, sep = "\n", append = TRUE)
cat(text = "", sep = "\n")

if (pair_strands_flag == 1) {
  system(paste0(PYCOQC, " -f ", d2_basecalling, "/sequencing_summary.txt -b ", d2_preprocessing, "/barcoding_summary.txt -o ", d2, "/qc/pycoQC_report.html"))
} else {
  system(paste0(PYCOQC, " -f ", d2_basecalling, "/sequencing_summary.txt -b ", d2_preprocessing, "/barcoding_summary.txt -o ", d2, "/qc/pycoQC_report.html"))
}
demu_files <- list.files(path = d2_preprocessing, pattern = "BC", full.names = TRUE)
for (i in 1:length(demu_files)) {
  BC_val_curr <- substr(x = basename(demu_files[i]), start = 3, stop = 4)
  if (paste0("BC", BC_val_curr) %in% BC_int) {
    cat(text = paste0("Now trimming adapters with Porechop for sample BC", BC_val_curr), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Now trimming adapters with Porechop for sample BC", BC_val_curr), sep = "\n")
    system(paste0(PORECHOP, " -i ", d2_preprocessing, "/BC", BC_val_curr, "_tmp1.fastq -b ", d2_preprocessing, "/BC", BC_val_curr, "_porechop_dir_tmp --require_two_barcodes --extra_end_trim ", primers_length))
    fastq_file_curr <- list.files(path = paste0(d2_preprocessing, "/BC", BC_val_curr, "_porechop_dir_tmp"), pattern = paste0("BC", BC_val_curr, "\\.fastq"), full.names = TRUE)
    if (length(fastq_file_curr) == 0) {
      BC_int <- setdiff(BC_int, paste0("BC", BC_val_curr))
      BC_trash <- c(BC_trash, paste0("BC", BC_val_curr))
      next
    }
    system(paste0("cp ", d2_preprocessing, "/BC", BC_val_curr, "_porechop_dir_tmp/BC", BC_val_curr, ".fastq ", d2_preprocessing, "/BC", BC_val_curr, "_tmp2.fastq"))
    system(paste0(SEQTK, " seq -A ", d2_preprocessing, "/BC", BC_val_curr, "_tmp2.fastq > ", d2_preprocessing, "/BC", BC_val_curr, "_tmp1.fasta"))
    sequences <- readDNAStringSet(paste0(d2_preprocessing, "/BC", BC_val_curr, "_tmp1.fasta"), "fasta")
    ws <- width(sequences)
    read_length <- ws
    cat(text = paste0("Mean read length (stdev) for sample BC", BC_val_curr, ": ", sprintf("%.0f", mean(ws)), " (", sprintf("%.0f", sd(ws)), ")"), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Mean read length (stdev) for sample BC", BC_val_curr, ": ", sprintf("%.0f", mean(ws)), " (", sprintf("%.0f", sd(ws)), ")"), sep = "\n")
    if (fixed_lenfil_flag == 1) {
      lb <- amplicon_length - lenfil_tol/2
      ub <- amplicon_length + lenfil_tol/2
    } else {
      lb <- mean(ws) - 2*sd(ws)
      ub <- mean(ws) + 2*sd(ws)
    }
    ws_ok <- ws[intersect(which(ws > lb), which(ws < ub))]
    read_length_ok <- ws_ok
    cat(text = paste0("Now filtering out reads shorter than ", sprintf("%.0f", lb), " and longer than ", sprintf("%.0f", ub), " bp for sample BC", BC_val_curr), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Now filtering out reads shorter than ", sprintf("%.0f", lb), " and longer than ", sprintf("%.0f", ub), " bp for sample BC", BC_val_curr), sep = "\n")
    system(paste0("cat ", d2_preprocessing, "/BC", BC_val_curr, "_tmp2.fastq | ", remove_long_short, " ", lb, " ", ub, " > ", d3, "/BC", BC_val_curr, ".fastq"))
    system(paste0(SEQTK, " seq -A ", d3, "/BC", BC_val_curr, ".fastq > ", d3, "/BC", BC_val_curr, ".fasta"))
    cat(text = paste0("Mean read length for sample BC", BC_val_curr, " after filtering: ", sprintf("%.0f", mean(ws_ok)), " (", sprintf("%.0f", sd(ws_ok)), ")"), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Mean read length for sample BC", BC_val_curr, " after filtering: ", sprintf("%.0f", mean(ws_ok)), " (", sprintf("%.0f", sd(ws_ok)), ")"), sep = "\n")
    png(paste0(d2, "/qc/hist_BC", BC_val_curr, "_unfiltered.png"))
    hist(read_length, main = paste("Read length before filtering - sample BC", BC_val_curr), col = "blue", breaks = 100, xlab = "Read length", ylab = "Number of reads")
    dev.off()
    png(paste0(d2, "/qc/hist_BC", BC_val_curr, ".png"))
    hist(read_length_ok, main = paste("Read length after filtering - sample BC", BC_val_curr), col = "blue", breaks = 100, xlab = "Read length", ylab = "Number of reads")
    dev.off()
    cat(text = "\n", file = logfile, append = TRUE)
  }
}

BC_files_bn <-  list.files(path = d3, pattern = "^BC\\d+\\.fasta$")
BC_files_bn_fq <- list.files(path = d3, pattern = "^BC\\d+\\.fastq$")
BC_files <- paste0(d3, "/", BC_files_bn)
BC_files_fq <- paste0(d3, "/", BC_files_bn_fq)

num_reads_BC_int <- 0
num_reads_vec <- vector(length = length(BC_files))

for (i in 1:length(BC_files)) {
  num_reads_vec[i] <- length(grep(x = readLines(BC_files[i]), pattern = "^>"))
  cat(text = paste0("Number of reads assigned to ", gsub("\\.fasta$", "", basename(BC_files_bn[i])), ": ", num_reads_vec[i]), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Number of reads assigned to ", gsub("\\.fasta$", "", basename(BC_files_bn[i])), ": ", num_reads_vec[i]), sep = "\n")
  num_reads_BC_int <- num_reads_BC_int + num_reads_vec[i]
  system(command = paste0("rm  ", BC_files[i]))
}

cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

for (i in 1:length(BC_files_fq)) {
  cat(text = paste0("Compressing ", basename(BC_files_fq[i]), " file"), file = logfile, sep = "\n", append = TRUE)
  cat(text = paste0("Compressing ", basename(BC_files_fq[i]), " file"), sep = "\n")
  system(command = paste0("gzip ", BC_files_fq[i]))
}

if (save_space_flag == 1) {
    cat(text = "\n", file = logfile, append = TRUE)
    cat(text = "\n")
    cat(text = paste0("Deleting ", d2_preprocessing, " directory"), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Deleting ", d2_preprocessing, " directory"), sep = "\n")
    system(command = paste0("rm -r ", d2_preprocessing))
    cat(text = "\n", file = logfile, append = TRUE)
    cat(text = "\n")
    cat(text = paste0("Deleting ", d2_basecalling, " directory"), file = logfile, sep = "\n", append = TRUE)
    cat(text = paste0("Deleting ", d2_basecalling, " directory"), sep = "\n")
    system(command = paste0("rm -r ", d2_basecalling))
}

cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

cat(text = paste0("Running the MetONTIIME pipeline"), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Running the MetONTIIME pipeline"), sep = "\n")
cat(text = "\n", file = logfile, append = TRUE)
cat(text = "\n")

system(paste0(MetONTIIME, " ", d3, " ", SAMPLE_METADATA, " ", DB, " ", TAXONOMY))

cat(text = paste0("Workflow ended at ", date(), "!"), file = logfile, sep = "\n", append = TRUE)  
cat(text = paste0("Workflow ended at ", date(), "!"), sep = "\n")
cat(text = paste0("Look at the results in directory ", d3), file = logfile, sep = "\n", append = TRUE)
cat(text = paste0("Look at the results in directory ", d3), sep = "\n")
