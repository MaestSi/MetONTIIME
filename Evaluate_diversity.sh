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

usage="$(basename "$0") [-w working_directory] [-m sample_metadata] [-d sampling_depth] [-t threads]"

while :
do
    case "$1" in
      -h | --help)
          echo $usage
          exit 0
          ;;
      -w)
          working_directory=$(realpath $2)
          shift 2
          echo "Working directory: $working_directory"
          ;;
      -m)
           sample_metadata=$(realpath $2)
           shift 2
           echo "Sample metadata: $sample_metadata"
           ;;
      -d)
           sampling_depth=$2
           shift 2
           echo "Minimum sampling depth: $sampling_depth reads"
           ;;
      -t)
           threads=$2
           shift 2
           echo "Number of threads: $threads"
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

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $working_directory/rep-seqs.qza \
  --o-alignment $working_directory/aligned-rep-seqs.qza \
  --o-masked-alignment $working_directory/masked-aligned-rep-seqs.qza \
  --o-tree $working_directory/unrooted-tree.qza \
  --o-rooted-tree $working_directory/rooted-tree.qza \
  --p-n-threads $threads

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $working_directory/rooted-tree.qza \
  --i-table $working_directory/table.qza \
  --p-sampling-depth $sampling_depth \
  --m-metadata-file $sample_metadata \
  --output-dir $working_directory/core-metrics-results \

qiime diversity alpha-rarefaction \
  --i-table $working_directory/table.qza \
  --i-phylogeny $working_directory/rooted-tree.qza \
  --p-max-depth $sampling_depth \
  --m-metadata-file $sample_metadata \
  --o-visualization $working_directory/alpha-rarefaction.qzv

