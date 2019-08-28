# MetONTIIME

**MetONTIIME** is a Meta-barcoding pipeline for analysing ONT data in QIIME2 framework. The whole bioinformatic workflow consists of a preprocessing pipeline and a script emulating EPI2ME 16S workflow, so to make the whole bioinformatic analysis from raw fast5 files to taxonomy assignments straightforward and simple.

## Getting started

**Prerequisites**

* Miniconda3.
Tested with conda 4.6.11.
```which conda``` should return the path to the executable.
If you don't have Miniconda3 installed, you could download and install it with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Then, after completing _MetONTIIME_ installation, set the _MINICONDA_DIR_ variable in **config_MinION_mobile_lab.R** to the full path to miniconda3 directory.

* Guppy, the software for basecalling and demultiplexing provided by ONT. Tested with Guppy v3.2.
If you don't have [Guppy](https://community.nanoporetech.com/downloads) installed, choose an appropriate version and install it.
For example, you could download and unpack the archive with:
```
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_version_of_interest.tar.gz
tar -xf ont-guppy-cpu_version_of_interest.tar.gz
```
A directory _ont-guppy-cpu_ should have been created in your current directory.
Then, after completing _MetONTIIME_ installation, set the _BASECALLER_DIR_ variable in **config_MinION_mobile_lab.R** to the full path to _ont-guppy-cpu/bin_ directory.

* A fasta file downloaded from NCBI that you want to use as a reference database. For example, if you want to use the same database used by the EPI2ME 16S workflow for bacterial 16S gene, you can go to [BioProject 33175](https://www.ncbi.nlm.nih.gov/nuccore?term=33175%5BBioProject%5D), click _send to_ and select _Complete Record_ and _File_.

**Installation**

```
git clone https://github.com/MaestSi/MetONTIIME.git
cd MetONTIIME
chmod 755 *
./install.sh
```

A conda environment named _MetONTIIME_env_ is created, where seqtk, porechop, pycoQC and qiime2-2019.7 are installed.
Then, you can open the **config_MinION_mobile_lab.R** file with a text editor and set the variables _PIPELINE_DIR_ and _MINICONDA_DIR_ to the value suggested by the installation step.

# Usage

The first time you run the _MetONTIIME_ pipeline on a new database, you can use the **Import_database.sh** script for importing a fasta file as a pair of _QIIME2_ artifacts. This script downloads some taxonomy files from NCBI (~9.4 GB) and uses [entrez qiime](https://github.com/bakerccm/entrez_qiime) and _QIIME2_ to generate a _DNAFASTAFormat_ and a _HeaderlessTSVTaxonomyFormat_ artifacts, containing sequences and corresponding taxonomy. _Entrez_qiime_ is installed to a new conda environment named _entrez_qiime_env_.
After this step, you can open the **config_MinION_mobile_lab.R** file with a text editor and set the variables _DB_ and _TAXONOMY_ to the newly generated _QIIME2_ artifacts. After that, you can run the full _MetONTIIME_ pipeline using the wrapper script **Launch_MinION_mobile_lab.sh**.

**Import_database.sh**

Usage: Import_database.sh \<"sample_name".fasta\>

Input:
* \<"sample_name".fasta\>: a fasta file downloaded from NCBI containing sequences that you want to use as a reference database

Outputs:
* \<"sample_name"\_sequence.qza\>: _QIIME2_ artifact of type _DNAFASTAFormat_ containing reference sequences
* \<"sample_name"\_taxonomy.qza\>: _QIIME2_ artifact of type _HeaderlessTSVTaxonomyFormat_ containing taxonomy of reference sequences

**Launch_MinION_mobile_lab.sh**

Usage: Launch_MinION_mobile_lab.sh \<fast5_dir\>

Input
* \<fast5_dir\>: directory containing raw fast5 files

Outputs (saved in \<fast5_dir\>_analysis/analysis):
* feature-table_absfreq.tsv: file containing the number of reads assigned to each taxa for each sample
* feature-table_relfreq.tsv: file containing the proportion of reads assigned to each taxa for each sample
* species_counts.txt: subset of feature-table_absfreq.tsv file, containing only species names and counts
* taxa-bar-plots.qzv: _QIIME2_ visualization artifact of barplots with taxonomy abundances
* demux_summary.qzv: _QIIME2_ visualization artifact with summary of sequences assigned to each sample after demultiplexing
* logfile.txt, manifest.txt, sequences.qza, table.qz*, rep-seqs.qz*, taxonomy.qz*, table_collapsed.qza, feature-table_absfreq.biom, table_collapsed_relfreq.qz*, feature-table_relfreq.biom: temporary files useful for debugging or for further analyses

Outputs (saved in \<fast5_dir\>_analysis/qc):
* Read length distributions and pycoQC report

Outputs (saved in \<fast5_dir\>_analysis/basecalling):
* Temporary files for basecalling

Outputs (saved in \<fast5_dir\>_analysis/preprocessing):
* Temporary files for demultiplexing, filtering based on read length and adapters trimming

## Auxiliary scripts

In the following, auxiliary scripts run by **Launch_MinION_mobile_lab.sh** are listed. These scripts should not be called directly.

**MinION_mobile_lab.R**

Note: script run by _Launch_MinION_mobile_lab.sh_.

**config_MinION_mobile_lab.R**

Note: configuration script, must be modified before running _Launch_MinION_mobile_lab.sh_.

**subsample_fast5.sh**

Note: script run by _MinION_mobile_lab.R_ if _do_subsampling_flag_ variable is set to 1 in _config_MinION_mobile_lab.R_.

**remove_long_short.pl**

Note: script run by _MinION_mobile_lab.R_ for removing reads of abnormal length.

# Results visualization

All .qzv and .qza artifacts can be visualized either importing them to [QIIME2 View](https://view.qiime2.org/) or with command:

```
source activate MetONTIIIME_env
qiime tools view <file.qz*>
```
