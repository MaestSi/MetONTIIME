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

DB_FASTA=$1

conda create -n entrez_qiime_env python=2.7 gcc_linux-64
source activate entrez_qiime_env

pip install numpy
pip install cogent

if [ ! -d "taxonomy" ]; then
  mkdir taxonomy
  cd taxonomy
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  gunzip nucl_gb.accession2taxid.gz
  mkdir taxdump
  cd taxdump
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  tar -zxvf taxdump.tar.gz
  cd ../..
fi

sed -i '/^[[:space:]]*$/d' $DB_FASTA

if [ ! -d "entrez_qiime" ]; then
  git clone https://github.com/bakerccm/entrez_qiime.git
fi
python ./entrez_qiime/entrez_qiime.py -i $DB_FASTA -n ./taxonomy/taxdump -a ./taxonomy/nucl_gb.accession2taxid
conda deactivate

source activate MetONTIIME_env

OUTPUT_DIR=$(dirname $(realpath $DB_FASTA))
DB_NAME=$(echo $(basename $DB_FASTA) | sed 's/\.fa.*$//g')
DB=$OUTPUT_DIR"/"$DB_NAME"_sequence.qza"
TAXONOMY_TSV=$OUTPUT_DIR"/"$DB_NAME"_accession_taxonomy.txt"
TAXONOMY=$OUTPUT_DIR"/"$DB_NAME"_taxonomy.qza"

qiime tools import \
--type FeatureData[Sequence] \
--input-format DNAFASTAFormat \
--input-path $DB_FASTA \
--output-path $DB

qiime tools import \
--type FeatureData[Taxonomy] \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $TAXONOMY_TSV \
--output-path $TAXONOMY
