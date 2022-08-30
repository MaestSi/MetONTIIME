#!/bin/bash

#
# Copyright 2021 Simone Maestri. All rights reserved.
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


usage="$(basename "$0") [-w working_dir] [-f metadata_file] [-s sequences_artifact] [-t taxonomy_artifact] [-n num_threads] [-c taxonomic_classifier] [-m max_accepts] [-q min_query_coverage] [-i min_id_thr]"

while :
do
    case "$1" in
      -h | --help)
          echo $usage
          exit 0
          ;;
      -w)
          WORKING_DIR=$(realpath $2)
          shift 2
          ;;
      -f)
           SAMPLE_METADATA=$2
           shift 2
           ;;
      -s)
           DB=$(realpath $2)
           shift 2
           ;;
      -t)
           TAXONOMY=$(realpath $2)
           shift 2
           ;;
      -n)
           THREADS=$2
           shift 2
           ;;
      -c)
           CLASSIFIER=$2
           shift 2
           ;;
      -m)
           MAX_ACCEPTS=$2
           shift 2
           ;;
      -q)
           QUERY_COV=$2
           shift 2
           ;;
      -i)
           ID_THR=$2
           shift 2
           ;;
       --) # End of all options
           shift
           break
           ;;
       -*)
           echo "Error: Unknown option: $1" >&2
           ## or call function display_help
           exit 1
           ;;
        *) # No more options
           break
           ;;
    esac
done

FASTQ_FILES=$(realpath $(find $WORKING_DIR -maxdepth 1 | grep "\.fastq\.gz"))
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
  --p-perc-identity $ID_THR \
  --p-query-cov $QUERY_COV \
  --p-maxaccepts $MAX_ACCEPTS \
  --o-classification taxonomy.qza \
  --o-search-results search_results.qza
elif [ "$CLASSIFIER_UC" == "VSEARCH" ]; then
  qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs.qza  \
    --i-reference-reads $DB \
    --i-reference-taxonomy $TAXONOMY \
    --p-perc-identity $ID_THR \
    --p-query-cov $QUERY_COV \
    --p-maxaccepts 100 \
    --p-maxrejects 100 \
    --p-maxhits $MAX_ACCEPTS \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --p-threads $THREADS \
    --o-classification taxonomy.qza \
    --o-search-results search_results.qza
else
  echo "Classifier $CLASSIFIER is not supported (choose between Blast and Vsearch); running default classifier Blast"
  qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza  \
  --i-reference-reads $DB \
  --i-reference-taxonomy $TAXONOMY \
  --p-perc-identity $ID_THR \
  --p-query-cov $QUERY_COV \
  --p-maxaccepts $MAX_ACCEPTS \
  --o-classification taxonomy.qza \
  --o-search-results search_results.qza
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

for lev in {1..7}; do
  qiime taxa collapse \
  --i-table table.qza --i-taxonomy taxonomy.qza \
  --p-level $lev \
  --o-collapsed-table table_collapsed_absfreq_level$lev.qza

qiime metadata tabulate  \
  --m-input-file table_collapsed_absfreq_level$lev.qza  \
  --o-visualization table_collapsed_absfreq_level$lev.qzv

  qiime tools export \
  --input-path table_collapsed_absfreq_level$lev.qza \
  --output-path .

  mv feature-table.biom feature-table_absfreq_level$lev.biom

  biom convert \
  -i feature-table_absfreq_level$lev.biom \
  -o feature-table_absfreq_level$lev.tsv \
  --to-tsv \
  --table-type 'Taxon table'

  qiime feature-table relative-frequency \
  --i-table table_collapsed_absfreq_level$lev.qza \
  --o-relative-frequency-table table_collapsed_relfreq_level$lev.qza

  qiime metadata tabulate  \
  --m-input-file table_collapsed_relfreq_level$lev.qza  \
  --o-visualization table_collapsed_relfreq_level$lev.qzv

  qiime tools export \
  --input-path table_collapsed_relfreq_level$lev.qza \
  --output-path .

  mv feature-table.biom feature-table_relfreq_level$lev.biom

  biom convert \
  -i feature-table_relfreq_level$lev.biom \
  -o feature-table_relfreq_level$lev.tsv \
  --to-tsv \
  --table-type 'Taxon table'

  cat feature-table_absfreq_level$lev.tsv | grep "#" > header
  cat feature-table_absfreq_level$lev.tsv | grep -v -P "Unassigned|#|BC|NA" | cut -f$lev -d';' | tr "\t" "," | grep -v "^__" | sort -t"," -k2,2 -nr | tr "," "\t" > "feature-table_absfreq_level"$lev"_stringent_noheader.tsv"
  cat header "feature-table_absfreq_level"$lev"_stringent_noheader.tsv" > "feature-table_absfreq_level"$lev"_stringent.tsv"
  rm "feature-table_absfreq_level"$lev"_stringent_noheader.tsv" header
done

mkdir collapsed_feature_tables
mv table_collapsed* feature-table_*freq_level*.biom feature-table_absfreq_level*_stringent.tsv collapsed_feature_tables
