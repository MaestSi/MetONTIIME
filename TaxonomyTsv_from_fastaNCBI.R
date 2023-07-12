#
# Copyright 2023 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@iit.it>
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

args <- commandArgs(trailingOnly = TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

suppressMessages(library("Biostrings"))
suppressMessages(library("taxize"))

sequences <- readDNAStringSet(dbSequencesFasta)
cat(sprintf("Imported %d fasta sequences\n", length(sequences)))

genbank_ids <- gsub(x = names(sequences), pattern = " .*", replacement = "")

max_attempts <- 10
num_uids_chunk <- 100
#process NBCI taxonomy UID in chunks
if (length(genbank_ids) < num_uids_chunk) {
  chunks_list <- list(1:length(genbank_ids))
} else {
  chunks_list <- split(1:length(genbank_ids), ceiling(seq(from = 1, to = length(genbank_ids))/num_uids_chunk))
}

if (exists("ENTREZ_KEY")) {
  cat("Converting GenBankID to NCBI taxonomy UID using the provided ENTREZ_KEY\n")
} else {
  cat("Converting GenBankID to NCBI taxonomy UID\n")
}

NCBI_taxa_uids_tmp <- c()
for (i in 1:length(chunks_list)) {
  cat(sprintf("Processing GenBankID chunk %d out of %d (%.2f%%)\n", i, length(chunks_list), i*100/length(chunks_list)))
  if (exists("ENTREZ_KEY")) {
    NCBI_taxa_uids_tmp_curr <- NULL
    attempt <- 1
    while(length(which(!is.na(NCBI_taxa_uids_tmp_curr))) != lapply(chunks_list, length)[[i]] && attempt <= max_attempts) {
      NCBI_taxa_uids_tmp_curr <- try(suppressMessages(genbank2uid(id = genbank_ids[chunks_list[[i]]], key = ENTREZ_KEY)))
      if (attempt > 1) cat(sprintf("Running attempt %d\n", attempt))
      attempt <- attempt + 1
    } 
  } else {
    NCBI_taxa_uids_tmp_curr <- NULL
    attempt <- 1
    while(length(which(!is.na(NCBI_taxa_uids_tmp_curr))) != lapply(chunks_list, length)[[i]] && attempt <= max_attempts) {
      NCBI_taxa_uids_tmp_curr <- try(suppressMessages(genbank2uid(id = genbank_ids[chunks_list[[i]]])))
      if (attempt > 1) cat(sprintf("Running attempt %d\n", attempt))
      attempt <- attempt + 1
    }
  }
  NCBI_taxa_uids_tmp <- c(NCBI_taxa_uids_tmp, NCBI_taxa_uids_tmp_curr)
}

NCBI_taxa_uids <- unlist(lapply(NCBI_taxa_uids_tmp, '[[', 1))
ind_fail <- which(is.na(NCBI_taxa_uids))
ind_ok <- which(!is.na(NCBI_taxa_uids))
if (length(ind_fail) > 0) {
  NCBI_taxa_uids_failed <- NCBI_taxa_uids_tmp[ind_fail]
}

NCBI_taxa_uids <- NCBI_taxa_uids[ind_ok]
NCBI_taxa_name_ok <- names(sequences)[ind_ok]
NCBI_taxa_name_fail <- names(sequences)[ind_fail]

if (exists("ENTREZ_KEY")) {
  cat("Retrieving taxonomy from NCBI taxonomy UID using the provided ENTREZ_KEY\n")
} else {
  cat("Retrieving taxonomy from NCBI taxonomy UID\n")
}

raw_classification <- c()
for (i in 1:length(chunks_list)) {
  cat(sprintf("Processing NCBI taxonomy UID chunk %d out of %d (%.2f%%)\n", i, length(chunks_list), i*100/length(chunks_list)))
  if (exists("ENTREZ_KEY")) {
    raw_classification_curr <- NULL
    attempt <- 1
    while(length(which(!is.na(raw_classification_curr))) != lapply(chunks_list, length)[[i]] && attempt <= max_attempts) {
      raw_classification_curr <- try(suppressMessages(classification(NCBI_taxa_uids[chunks_list[[i]]], db = 'ncbi', key = ENTREZ_KEY)))
      if (attempt > 1) cat(sprintf("Running attempt %d\n", attempt))
      attempt <- attempt + 1
    } 
  } else {
    raw_classification_curr <- NULL
    attempt <- 1
    while(length(which(!is.na(raw_classification_curr))) != lapply(chunks_list, length)[[i]] && attempt <= max_attempts) {
      raw_classification_curr <- try(suppressMessages(classification(NCBI_taxa_uids[chunks_list[[i]]], db = 'ncbi')))
      if (attempt > 1)  cat(sprintf("Running attempt %d\n", attempt))
      attempt <- attempt + 1
    }
  }
  raw_classification <- c(raw_classification, raw_classification_curr)
}

cat("Reformatting taxonomy\n")
full_taxonomy <- c()
for (i in 1:length(raw_classification)) {
  if (i %% 100 == 0) {
    cat(sprintf("Reformatted taxonomy for %d entries out of %d (%.2f%%)\n", i, length(sequences), 100*i/length(sequences)))
  }
  if (!is.null(ncol(raw_classification[[i]]))) {
    if (nrow(raw_classification[[i]]) > 1) {
      subspecies_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "subspecies"), 1]
      species_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "species"), 1]
      genus_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "genus"), 1]
      family_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "family"), 1]
      order_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "order"), 1]
      class_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "class"), 1]
      phylum_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "phylum"), 1]
      kingdom_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "superkingdom"), 1]
      full_taxonomy[i] <- paste(kingdom_name_curr, phylum_name_curr, class_name_curr, order_name_curr, family_name_curr, genus_name_curr, species_name_curr, subspecies_name_curr, sep = ";")
    } else {
      full_taxonomy[i] <- "Unclassified;;;;;;;"
    }
  } else {
    full_taxonomy[i] <- "Unclassified;;;;;;;"
  }
}
names(full_taxonomy) <- NCBI_taxa_name_ok

full_taxonomy_final <- vector(mode = "character", length = length(sequences))
names(full_taxonomy_final) <- names(sequences)
full_taxonomy_final[names(sequences)[ind_ok]] <- full_taxonomy[NCBI_taxa_name_ok]

if (length(ind_fail) > 0) {
  cat(sprintf("Failed to retrieve classification for %d entries out of %d (%.2f%%)\n", length(NCBI_taxa_uids_failed), length(sequences), 100*length(NCBI_taxa_uids_failed)/length(sequences)))
  full_taxonomy_final[names(sequences)[ind_fail]] <- NCBI_taxa_name_fail
}

seqnames_reduced <- unlist(lapply(strsplit(x = names(sequences), split = " "), '[[', 1))
full_taxonomy_df <- data.frame(seqnames=seqnames_reduced, taxonomy=unname(full_taxonomy_final))
cat("Writing dbTaxonomyTsv to file\n")
write.table(x = full_taxonomy_df, file = dbTaxonomyTsv, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
