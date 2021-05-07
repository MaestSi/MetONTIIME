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

usage="$(basename "$0") [-f FEATURE_TABLE] [-m SAMPLE_METADATA] [-d SAMPLING_DEPTH]"

while :
do
    case "$1" in
      -h | --help)
          echo $usage
          exit 0
          ;;
      -f)
          FEATURE_TABLE=$(realpath $2)
          shift 2
          echo "Collapsed feature table: $FEATURE_TABLE"
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

WORKING_DIRECTORY=$(realpath $(dirname $FEATURE_TABLE))
TN=$(echo $(basename $FEATURE_TABLE) | sed 's/\.qza//')

qiime diversity core-metrics \
  --i-table $FEATURE_TABLE \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file $SAMPLE_METADATA \
  --output-dir $WORKING_DIRECTORY"/core-metrics-results_"$TN"_"$SAMPLING_DEPTH"_subsampled_non_phylogenetic" \

for f in $(find $WORKING_DIRECTORY"/core-metrics-results_"$TN"_"$SAMPLING_DEPTH"_subsampled_non_phylogenetic" | grep "\.qza" | grep -v "distance_matrix"); do
  fn=$(echo $f | sed 's/\.qza//');
  qiime metadata tabulate --m-input-file $f --o-visualization $fn".qzv";
done

qiime diversity alpha-rarefaction \
  --i-table $FEATURE_TABLE \
  --p-max-depth $SAMPLING_DEPTH \
  --m-metadata-file $SAMPLE_METADATA \
  --o-visualization $WORKING_DIRECTORY"/alpha-rarefaction_"$TN"_"$SAMPLING_DEPTH"_subsampled_non_phylogenetic.qzv"
