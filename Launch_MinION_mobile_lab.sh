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

RAW_READS_DIR=$1
RAW_READS_DIR_FULL=$(realpath $RAW_READS_DIR)
source activate MetONTIIME_env
PIPELINE_DIR=$(realpath $( dirname "${BASH_SOURCE[0]}" ))
nohup Rscript $PIPELINE_DIR/MinION_mobile_lab.R $PIPELINE_DIR/config_MinION_mobile_lab.R $RAW_READS_DIR_FULL &
