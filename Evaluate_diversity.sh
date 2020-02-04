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

usage="$(basename "$0") [-w WORKING_DIRECTORY] [-m SAMPLE_METADATA] [-d SAMPLING_DEPTH] [-t THREADS] [-c CLUSTERING_THRESHOLD]"

while :
do
    case "$1" in
      -h | --help)
          echo $usage
          exit 0
          ;;
      -w)
          WORKING_DIRECTORY=$(realpath $2)
          shift 2
          echo "Working directory: $WORKING_DIRECTORY"
          ;;
      -m)
           SAMPLE_METADATA=$(realpath $2)
           shift 2
           echo "Sample metadata: $SAMPLE_METADATA"
           ;;
      -d)
           SAMPLING_DEPTH=$2
           shift 2
           echo "Sampling depth: $SAMPLING_DEPTH reads"
           ;;
      -t)
           THREADS=$2
           shift 2
           echo "Number of threads: $THREADS"
           ;;
      -c)
           CLUSTERING_THRESHOLD=$2
           shift 2
           echo "Clustering threshold: $CLUSTERING_THRESHOLD"
           if (( $(bc <<<"$CLUSTERING_THRESHOLD > 1") )); then
             echo "Choose a value for -c parameter in (0, 1]"
             exit 1
           fi
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


source activate MetONTIIME_env

FASTQ_FILES=$(realpath $(find $WORKING_DIRECTORY -maxdepth 1 | grep "\.fastq\.gz" | grep -v "subsampled\.fastq\.gz"))

MANIFEST_SUB=$WORKING_DIRECTORY"/manifest_"$SAMPLING_DEPTH"_subsampled.txt"

if [ ! -f "$MANIFEST_SUB" ]; then
  echo -e sample-id"\t"absolute-filepath > $MANIFEST_SUB
  for f in $FASTQ_FILES; do
    s=$(echo $(basename $f) | sed 's/\.fastq\.gz//g')
    s_sub=$s"_"$SAMPLING_DEPTH"_subsampled"
    f_sub=$WORKING_DIRECTORY"/"$s_sub".fastq.gz"
    seqtk sample $f $SAMPLING_DEPTH | gzip > $f_sub
    echo -e $s"\t"$f_sub >> $MANIFEST_SUB;
  done
fi

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path $MANIFEST_SUB \
--input-format 'SingleEndFastqManifestPhred33V2' \
--output-path "sequences_"$SAMPLING_DEPTH"_subsampled.qza"

qiime vsearch dereplicate-sequences \
--i-sequences "sequences_"$SAMPLING_DEPTH"_subsampled.qza" \
--o-dereplicated-table "table_tmp_"$SAMPLING_DEPTH"_subsampled.qza" \
--o-dereplicated-sequences "rep-seqs_tmp_"$SAMPLING_DEPTH"_subsampled.qza"

qiime vsearch cluster-features-de-novo \
--i-sequences "rep-seqs_tmp_"$SAMPLING_DEPTH"_subsampled.qza" \
--i-table "table_tmp_"$SAMPLING_DEPTH"_subsampled.qza" \
--p-perc-identity $CLUSTERING_THRESHOLD \
--o-clustered-table "table_"$SAMPLING_DEPTH"_subsampled.qza" \
--o-clustered-sequences "rep-seqs_"$SAMPLING_DEPTH"_subsampled.qza" \
--p-threads $THREADS

rm "table_tmp_"$SAMPLING_DEPTH"_subsampled.qza" "rep-seqs_tmp_"$SAMPLING_DEPTH"_subsampled.qza"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $WORKING_DIRECTORY"/rep-seqs_"$SAMPLING_DEPTH"_subsampled.qza" \
  --o-alignment $WORKING_DIRECTORY"/aligned-rep-seqs_"$SAMPLING_DEPTH"_subsampled.qza" \
  --o-masked-alignment $WORKING_DIRECTORY"/masked-aligned-rep-seqs_"$SAMPLING_DEPTH"_subsampled.qza" \
  --o-tree $WORKING_DIRECTORY"/unrooted-tree_"$SAMPLING_DEPTH"_subsampled.qza" \
  --o-rooted-tree $WORKING_DIRECTORY"/rooted-tree_"$SAMPLING_DEPTH"_subsampled.qza" \
  --p-n-threads $THREADS

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $WORKING_DIRECTORY"/rooted-tree_"$SAMPLING_DEPTH"_subsampled.qza" \
  --i-table $WORKING_DIRECTORY"/table_"$SAMPLING_DEPTH"_subsampled.qza" \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file $SAMPLE_METADATA \
  --output-dir $WORKING_DIRECTORY"/core-metrics-results_"$SAMPLING_DEPTH"_subsampled" \

qiime diversity alpha-rarefaction \
  --i-table $WORKING_DIRECTORY"/table_"$SAMPLING_DEPTH"_subsampled.qza" \
  --i-phylogeny $WORKING_DIRECTORY"/rooted-tree_"$SAMPLING_DEPTH"_subsampled.qza" \
  --p-max-depth $SAMPLING_DEPTH \
  --m-metadata-file $SAMPLE_METADATA \
  --o-visualization $WORKING_DIRECTORY"/alpha-rarefaction_"$SAMPLING_DEPTH"_subsampled.qzv"
