#!/bin/bash

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

WORKING_DIR=$1
SAMPLE_METADATA=$2
DB=$3
TAXONOMY=$4
THREADS=$5
CLASSIFIER=$6

FASTQ_FILES=$(realpath $(find $WORKING_DIR | grep "\.fastq\.gz"))
MANIFEST=$WORKING_DIR"/manifest.txt"

if [ ! -f "$MANIFEST" ]; then
  echo -e sample-id"\t"absolute-filepath > $MANIFEST
  for f in $FASTQ_FILES; do
    s=$(echo $(basename $f) | sed 's/\.fastq\.gz//g')
    echo -e $s"\t"$f >> $MANIFEST;
  done
fi

if [ ! -f "$SAMPLE_METADATA" ]; then
  echo -e sample-id"\t"sample-name > $SAMPLE_METADATA
  for f in $FASTQ_FILES; do
    s=$(echo $(basename $f) | sed 's/\.fastq\.gz//g')
    echo -e $s"\t"$s >> $SAMPLE_METADATA;
  done
fi

cd $WORKING_DIR

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path $MANIFEST \
--input-format 'SingleEndFastqManifestPhred33V2' \
--output-path sequences.qza

qiime vsearch dereplicate-sequences \
--i-sequences sequences.qza \
--o-dereplicated-table table_tmp.qza \
--o-dereplicated-sequences rep-seqs_tmp.qza

qiime vsearch cluster-features-de-novo \
--i-sequences rep-seqs_tmp.qza \
--i-table table_tmp.qza \
--p-perc-identity 1 \
--o-clustered-table table.qza \
--o-clustered-sequences rep-seqs.qza \
--p-threads $THREADS

rm table_tmp.qza rep-seqs_tmp.qza

qiime demux summarize \
--i-data sequences.qza \
--o-visualization demux_summary.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file $SAMPLE_METADATA

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

CLASSIFIER_UC=$(awk '{print toupper($0)'} <<< $CLASSIFIER)

if [ "$CLASSIFIER_UC" == "BLAST" ]; then
  qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza  \
  --i-reference-reads $DB \
  --i-reference-taxonomy $TAXONOMY \
  --p-perc-identity 0.77 \
  --p-query-cov 0.3 \
  --p-maxaccepts 1 \
  --o-classification taxonomy.qza
elif [ "$CLASSIFIER_UC" == "VSEARCH" ]; then
  qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs.qza  \
    --i-reference-reads $DB \
    --i-reference-taxonomy $TAXONOMY \
    --p-perc-identity 0.77 \
    --p-query-cov 0.3 \
    --p-top-hits-only \
    --p-maxaccepts 1 \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --p-threads $THREADS \
    --o-classification taxonomy.qza
else
  echo "Classifier $CLASSIFIER is not supported (choose between Blast and Vsearch); running default classifier Blast"
  qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza  \
  --i-reference-reads $DB \
  --i-reference-taxonomy $TAXONOMY \
  --p-perc-identity 0.77 \
  --p-query-cov 0.3 \
  --p-maxaccepts 1 \
  --o-classification taxonomy.qza
fi

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file $SAMPLE_METADATA \
  --o-visualization taxa-bar-plots.qzv

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Unassigned \
  --o-filtered-table table-no-Unassigned.qza

qiime taxa barplot \
  --i-table table-no-Unassigned.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file $SAMPLE_METADATA \
  --o-visualization taxa-bar-plots-no-Unassigned.qzv

qiime taxa collapse \
--i-table table.qza --i-taxonomy taxonomy.qza \
--p-level 6 \
--o-collapsed-table table_collapsed.qza

qiime tools export \
--input-path table_collapsed.qza \
--output-path .

mv feature-table.biom feature-table_absfreq.biom

biom convert \
-i feature-table_absfreq.biom \
-o feature-table_absfreq.tsv \
--to-tsv \
--table-type 'Taxon table'

qiime feature-table relative-frequency \
--i-table table_collapsed.qza \
--o-relative-frequency-table table_collapsed_relfreq.qza

qiime metadata tabulate  \
--m-input-file table_collapsed_relfreq.qza  \
--o-visualization table_collapsed_relfreq.qzv

qiime tools export \
--input-path table_collapsed_relfreq.qza \
--output-path .

mv feature-table.biom feature-table_relfreq.biom

biom convert \
-i feature-table_relfreq.biom \
-o feature-table_relfreq.tsv \
--to-tsv \
--table-type 'Taxon table'

#NCBI database
cat feature-table_absfreq.tsv | grep "#" > header
cat feature-table_absfreq.tsv | grep -v -P "Unassigned|#|BC|NA" | cut -f6 -d';' | tr "\t" "," | sort -t"," -k2,2 -nr | tr "," "\t" > species_counts_noheader.txt

cat header species_counts_noheader.txt > species_counts.txt

rm species_counts_noheader.txt header
