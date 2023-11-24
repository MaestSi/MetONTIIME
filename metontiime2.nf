#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/MetONTIIME2
========================================================================================
 MaestSi/MetONTIIME2 analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/MetONTIIME2
----------------------------------------------------------------------------------------
*/

def helpMessage() {
		log.info"""
	Usage:
	nextflow -c metontiime2.conf run metontiime2.nf --workDir="/path/to/workDir" --resultsDir="/path/to/resultsDir" -profile docker

	Mandatory argument:
	-profile                                                 Configuration profile to use. Available: docker, singularity
	Other mandatory arguments which may be specified in the metontiime2.conf file

	--workDir                                               Path to working directory including fastq.gz files
	--sampleMetadata                                        Path to sample metadata tsv file; if it doesn't exist yet, it is created at runtime
	--dbSequencesFasta                                      Path to database file with sequences in fasta format
	--dbTaxonomyTsv                                         Path to database file with sequence id-to-taxonomy correspondence in tsv format
	--dbSequencesQza                                        Database file name with sequences as QIIME2 artifact (qza)
	--dbTaxonomyQza                                         Database file name with sequence id-to-taxonomy correspondence as QIIME2 artifact (qza)
	--classifier                                            Taxonomy classifier, available: VSEARCH, Blast
	--maxNumReads                                           Maximum number of reads per sample; if one sample has more than maxNumReads, random downsampling is performed
	--minReadLength                                         Minimum length (bp) for a read to be retained
	--maxReadLength                                         Maximum length (bp) for a read to be retained
	--minQual                                               Minimum average PHRED score for a read to be retained
	--extraEndsTrim                                         Number of bases to be trimmed at both ends
	--clusteringIdentity                                    Identity for de novo clustering [0-1]
	--maxAccepts                                            Maximum number of candidate hits for each read, to be used for consensus taxonomy assignment
	--minConsensus                                          Minimum fraction of assignments must match top hit to be accepted as consensus assignment [0.5-1]
	--minQueryCoverage                                      Minimum query coverage for an alignment to be considered a candidate hit [0-1]
	--minIdentity                                           Minimum alignment identity for an alignment to be considered a candidate hit [0-1]
	--taxaLevelDiversity                                    Taxonomy level at which you want to perform non phylogeny-based diversity analyses
	--numReadsDiversity                                     Max num. reads for diversity analyses
	--taxaOfInterest                                        Taxa of interest that you want to retain and to focus the analysis on
	--minNumReadsTaxaOfInterest                             Minimum number of reads assigned to Taxa of interest to retain a sample
	--resultsDir                                            Path to directory containing results
	""".stripIndent()
}

// Show help message
if (params.help) {
	helpMessage()
	exit 0
}

//import db
process importDb {
	input:
	path sequences_db
	path taxonomy_db
	output:
	val 'flag_db'
	script:
    if(params.importDb)
	"""
		mkdir -p ${params.resultsDir}/importDb

		qiime tools import \
		--type \'FeatureData[Sequence]\' \
		--input-path $sequences_db \
		--output-path ${params.resultsDir}/importDb/${params.dbSequencesQza}

		qiime tools import \
		--type \'FeatureData[Taxonomy]\' \
		--input-path $taxonomy_db \
		--input-format HeaderlessTSVTaxonomyFormat \
		--output-path ${params.resultsDir}/importDb/${params.dbTaxonomyQza}
	"""
	else
	"""
		mkdir -p ${params.resultsDir}/importDb
	"""
}

//concatenate fastq files
process concatenateFastq {
	input:
	val workdir
	output:
	val 'flag_concatenate'
	script:
    if(params.concatenateFastq)
	"""
		mkdir -p ${params.resultsDir}
		mkdir -p ${params.resultsDir}/concatenateFastq

		barcodes_dirs=\$(realpath \$(find $workdir -maxdepth 1 -mindepth 1 | grep \"barcode\"));
		echo \$barcodes_dirs
		for b in \$barcodes_dirs; do
			bn=\$(basename \$b);
			f=\$(find \$b | grep \".fastq\");
			zless \$f | gzip > ${params.resultsDir}/concatenateFastq/\$bn.fastq.gz;
		done
	"""
	else
	"""
		mkdir -p ${params.resultsDir}
		mkdir -p ${params.resultsDir}/concatenateFastq
		cp $workdir/*fastq.gz ${params.resultsDir}/concatenateFastq
	"""
}

//filter fastq files by quality and length
process filterFastq {
	input:
	val 'flag_concatenate'
	output:
	val 'flag_filter'
	script:
    if(params.filterFastq)
	"""
		mkdir -p ${params.resultsDir}/filterFastq
		
		fq=\$(find ${params.resultsDir}/concatenateFastq | grep \"\\.fastq\");
		
		for f in \$fq; do 
			sn=\$(echo \$(basename \$f) | sed \'s/\\.fastq.*//\');
			zless \$f | NanoFilt -q ${params.minQual} -l ${params.minReadLength} --maxlength ${params.maxReadLength} | seqtk trimfq -b ${params.extraEndsTrim} -e ${params.extraEndsTrim} - | gzip > ${params.resultsDir}/filterFastq/\$sn.fastq.gz;
			if LC_ALL=C gzip -l ${params.resultsDir}/filterFastq/\$sn.fastq.gz | awk \'NR==2 {exit(\$2!=0)}\'; then rm ${params.resultsDir}/filterFastq/\$sn.fastq.gz; fi
		done
	"""
	else
	"""
		mkdir -p ${params.resultsDir}/filterFastq
		cp ${params.resultsDir}/concatenateFastq/*fastq.gz ${params.resultsDir}/filterFastq
	"""
}

//downsample fastq files
process downsampleFastq {
	input:
	val 'flag_filter'
	output:
	val 'flag_downsample'
	script:
    if(params.downsampleFastq)
	"""
		mkdir -p ${params.resultsDir}/downsampleFastq
		fq=\$(find ${params.resultsDir}/filterFastq/ | grep \"\\.fastq\\.gz\$\");
		for f in \$fq; do 
			sn=\$(basename \$f);
			seqtk sample \$f ${params.maxNumReads} | gzip > ${params.resultsDir}/downsampleFastq/\$sn
		done
	"""
	else
	"""
		mkdir -p ${params.resultsDir}/downsampleFastq
		cp ${params.resultsDir}/filterFastq/*fastq.gz ${params.resultsDir}/downsampleFastq
	"""
}

//import fastq files as QIIME2 artifact
process importFastq {
	input:
	val 'flag_downsample'
	output:
	val 'flag_import'
	script:
    if(params.importFastq)
	"""
		mkdir -p ${params.resultsDir}/importFastq
		
		fq=\$(realpath \$(find ${params.resultsDir}/downsampleFastq/ -maxdepth 1 | grep \"\\.fastq\\.gz\"))
		manifestFile=${params.resultsDir}/importFastq/manifest.txt

		if [ ! -f "\$manifestFile" ]; then
			echo -e sample-id\"\t\"absolute-filepath > \$manifestFile;
			for f in \$fq; do
			s=\$(echo \$(basename \$f) | sed \'s/\\.fastq\\.gz//g\');
			echo -e \$s\"\t\"\$f >> \$manifestFile;
			done
		fi

		if [ ! -f ${params.sampleMetadata} ]; then
			echo -e sample-id\"\t\"sample-name > ${params.sampleMetadata};
			for f in \$fq; do
			s=\$(echo \$(basename \$f) | sed \'s/\\.fastq\\.gz//g\');
			echo -e \$s\"\t\"\$s >> ${params.sampleMetadata};
			done
		fi

		qiime tools import \
		--type 'SampleData[SequencesWithQuality]' \
		--input-path \$manifestFile \
		--input-format 'SingleEndFastqManifestPhred33V2' \
		--output-path ${params.resultsDir}/importFastq/sequences.qza

		ln -s ${params.resultsDir}/importFastq/sequences.qza ./sequences.qza
	"""
	else
	"""
		echo "Skipped"
	"""
}

//perform data QC
process dataQC {
	input:
	val 'flag_import'
	output:
	script:
    if(params.dataQC)
	"""
		mkdir -p ${params.resultsDir}/dataQC

		qiime demux summarize \
		--i-data ${params.resultsDir}/importFastq/sequences.qza \
		--o-visualization ${params.resultsDir}/dataQC/demux_summary.qzv
	"""
	else
	"""
		echo "Skipped"
	"""
}

//dereplicate sequences
process derepSeq {
	input:
	val 'flag_import'
	output:
	val 'flag_derep'
	script:
    if(params.derepSeq)
	"""
		mkdir -p ${params.resultsDir}/derepSeq

		qiime vsearch dereplicate-sequences \
		--i-sequences ${params.resultsDir}/importFastq/sequences.qza \
		--o-dereplicated-table ${params.resultsDir}/derepSeq/table_tmp.qza \
		--o-dereplicated-sequences ${params.resultsDir}/derepSeq/rep-seqs_tmp.qza

		qiime vsearch cluster-features-de-novo \
		--i-sequences ${params.resultsDir}/derepSeq/rep-seqs_tmp.qza \
		--i-table ${params.resultsDir}/derepSeq/table_tmp.qza \
		--p-perc-identity ${params.clusteringIdentity} \
		--o-clustered-table ${params.resultsDir}/derepSeq/table.qza \
		--o-clustered-sequences ${params.resultsDir}/derepSeq/rep-seqs.qza \
		--p-threads ${task.cpus}

		rm ${params.resultsDir}/derepSeq/table_tmp.qza ${params.resultsDir}/derepSeq/rep-seqs_tmp.qza

		qiime feature-table summarize \
		--i-table ${params.resultsDir}/derepSeq/table.qza \
		--o-visualization ${params.resultsDir}/derepSeq/table.qzv \
		--m-sample-metadata-file ${params.sampleMetadata}

		qiime feature-table tabulate-seqs \
		--i-data ${params.resultsDir}/derepSeq/rep-seqs.qza \
		--o-visualization ${params.resultsDir}/derepSeq/rep-seqs.qzv
	"""
	else
	"""
		echo "Skipped"
	"""
}

//assign taxonomy
process assignTaxonomy {
	input:
	val 'flag_derep'
	val 'flag_db'
	output:
	val 'flag_taxa'
	script:
    if(params.assignTaxonomy)
	"""
		mkdir -p ${params.resultsDir}/assignTaxonomy

		classifier_uc=\$(awk '{print toupper(\$0)'} <<< ${params.classifier})

		if [ "\$classifier_uc" == "BLAST" ]; then
			qiime feature-classifier makeblastdb \
			--i-sequences ${params.resultsDir}/importDb/${params.dbSequencesQza} \
			--o-database ${params.resultsDir}/importDb/blastIndexedDb.qza
			
			qiime feature-classifier classify-consensus-blast \
			--i-query ${params.resultsDir}/derepSeq/rep-seqs.qza  \
			--i-blastdb ${params.resultsDir}/importDb/blastIndexedDb.qza \
			--i-reference-taxonomy ${params.resultsDir}/importDb/${params.dbTaxonomyQza} \
			--p-num-threads ${task.cpus} \
			--p-perc-identity ${params.minIdentity} \
			--p-query-cov ${params.minQueryCoverage} \
			--p-maxaccepts ${params.maxAccepts} \
			--p-min-consensus ${params.minConsensus} \
			--o-classification ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
			--o-search-results ${params.resultsDir}/assignTaxonomy/search_results.qza
		elif [ "\$classifier_uc" == "VSEARCH" ]; then
			qiime feature-classifier classify-consensus-vsearch \
			--i-query ${params.resultsDir}/derepSeq/rep-seqs.qza  \
			--i-reference-reads ${params.resultsDir}/importDb/${params.dbSequencesQza} \
			--i-reference-taxonomy ${params.resultsDir}/importDb/${params.dbTaxonomyQza} \
			--p-perc-identity ${params.minIdentity} \
			--p-query-cov ${params.minQueryCoverage} \
			--p-maxaccepts 100 \
			--p-maxrejects 100 \
			--p-maxhits ${params.maxAccepts} \
			--p-min-consensus ${params.minConsensus} \
			--p-strand 'both' \
			--p-unassignable-label 'Unassigned' \
			--p-threads ${task.cpus} \
			--o-classification ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
			--o-search-results ${params.resultsDir}/assignTaxonomy/search_results.qza
		else
			echo "Classifier ${params.classifier} is not supported (choose between Blast and Vsearch)"
		fi

	qiime metadata tabulate \
		--m-input-file ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--o-visualization ${params.resultsDir}/assignTaxonomy/taxonomy.qzv

		qiime taxa filter-table \
		--i-table ${params.resultsDir}/derepSeq/table.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--p-exclude Unassigned \
		--o-filtered-table ${params.resultsDir}/derepSeq/table-no-Unassigned.qza
	"""
	else
	"""
		echo "Skipped"
	"""
}

//collapse feature tables
process collapseTables {
	input:
	val 'flag_taxa'
	output:
	val 'flag_collapse'
	script:
    if(params.collapseTables)
	"""
		mkdir -p ${params.resultsDir}/collapseTables

		num_levels=\$(echo \$(cat ${params.dbTaxonomyTsv} | head -n2 | tail -n1 | cut -f2 | grep -n -o \";\" | wc -l) + 1 | bc)

		for lev in \$( seq 1 \$num_levels); do
			qiime taxa collapse \
			--i-table ${params.resultsDir}/derepSeq/table.qza \
			--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
			--p-level \$lev \
			--o-collapsed-table ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level\$lev.qza

			qiime metadata tabulate  \
			--m-input-file ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level\$lev.qza \
			--o-visualization ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level\$lev.qzv

			qiime tools export \
			--input-path ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level\$lev.qza \
			--output-path ${params.resultsDir}/collapseTables/

			mv ${params.resultsDir}/collapseTables/feature-table.biom ${params.resultsDir}/collapseTables/feature-table-absfreq-level\$lev.biom

			biom convert \
			-i ${params.resultsDir}/collapseTables/feature-table-absfreq-level\$lev.biom \
			-o ${params.resultsDir}/collapseTables/feature-table-absfreq-level\$lev.tsv \
			--to-tsv \
			--table-type 'Taxon table'

			qiime feature-table relative-frequency \
			--i-table ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level\$lev.qza \
			--o-relative-frequency-table ${params.resultsDir}/collapseTables/table-collapsed-relfreq-level\$lev.qza

			qiime metadata tabulate  \
			--m-input-file ${params.resultsDir}/collapseTables/table-collapsed-relfreq-level\$lev.qza \
			--o-visualization ${params.resultsDir}/collapseTables/table-collapsed-relfreq-level\$lev.qzv

			qiime tools export \
			--input-path ${params.resultsDir}/collapseTables/table-collapsed-relfreq-level\$lev.qza \
			--output-path ${params.resultsDir}/collapseTables/

			mv ${params.resultsDir}/collapseTables/feature-table.biom ${params.resultsDir}/collapseTables/feature-table-relfreq-level\$lev.biom

			biom convert \
			-i ${params.resultsDir}/collapseTables/feature-table-relfreq-level\$lev.biom \
			-o ${params.resultsDir}/collapseTables/feature-table-relfreq-level\$lev.tsv \
			--to-tsv \
			--table-type 'Taxon table'
		done
	"""
	else
	"""
		echo "Skipped"
	"""
}

//filter taxa of interest
process filterTaxa {
	input:
	val 'flag_taxa'
	output:
	val 'flag_filter'
	script:
    if(params.filterTaxa)
	"""
		mkdir -p ${params.resultsDir}/filterTaxa

		qiime taxa filter-table \
		--i-table ${params.resultsDir}/derepSeq/table.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--p-include ${params.taxaOfInterest} \
		--o-filtered-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}.qza

		qiime feature-table filter-samples \
		--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}.qza	\
		--p-min-frequency ${params.minNumReadsTaxaOfInterest} \
		--o-filtered-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}.qza

		num_levels=\$(echo \$(cat ${params.dbTaxonomyTsv} | head -n2 | tail -n1 | cut -f2 | grep -n -o \";\" | wc -l) + 1 | bc)

		for lev in \$( seq 1 \$num_levels); do
			qiime taxa collapse \
			--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}.qza \
			--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
			--p-level \$lev \
			--o-collapsed-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level\$lev.qza

			qiime metadata tabulate \
			--m-input-file ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level\$lev.qza \
			--o-visualization ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level\$lev.qzv

			qiime tools export \
			--input-path ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level\$lev.qza \
			--output-path ${params.resultsDir}/filterTaxa/

			mv ${params.resultsDir}/filterTaxa/feature-table.biom ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-absfreq-level\$lev.biom

			biom convert \
			-i ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-absfreq-level\$lev.biom \
			-o ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-absfreq-level\$lev.tsv \
			--to-tsv \
			--table-type 'Taxon table'

			qiime feature-table relative-frequency \
			--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level\$lev.qza \
			--o-relative-frequency-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-relfreq-level\$lev.qza

			qiime metadata tabulate \
			--m-input-file ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-relfreq-level\$lev.qza \
			--o-visualization ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-relfreq-level\$lev.qzv

			qiime tools export \
			--input-path ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-relfreq-level\$lev.qza \
			--output-path ${params.resultsDir}/filterTaxa/

			mv ${params.resultsDir}/filterTaxa/feature-table.biom ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-relfreq-level\$lev.biom

			biom convert \
			-i ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-relfreq-level\$lev.biom \
			-o ${params.resultsDir}/filterTaxa/feature-table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-relfreq-level\$lev.tsv \
			--to-tsv \
			--table-type 'Taxon table'
		done
	"""
	else
	"""
		echo "Skipped"
	"""
}

//produce barplots
process taxonomyVisualization {
	input:
	val 'flag_filter'
	output:
	script:
    if(params.taxonomyVisualization && params.filterTaxa)
	"""
		mkdir -p ${params.resultsDir}/taxonomyVisualization

		qiime taxa barplot \
		--i-table ${params.resultsDir}/derepSeq/table.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/taxonomyVisualization/taxa-bar-plots.qzv

		qiime taxa barplot \
		--i-table ${params.resultsDir}/derepSeq/table-no-Unassigned.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/taxonomyVisualization/taxa-bar-plots-no-Unassigned.qzv

		qiime taxa barplot \
		--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/taxonomyVisualization/taxa-bar-plots-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}.qzv
	"""
	else if(params.taxonomyVisualization && !params.filterTaxa)
	"""	
		mkdir -p ${params.resultsDir}/taxonomyVisualization

		qiime taxa barplot \
		--i-table ${params.resultsDir}/derepSeq/table.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/taxonomyVisualization/taxa-bar-plots.qzv

		qiime taxa barplot \
		--i-table ${params.resultsDir}/derepSeq/table-no-Unassigned.qza \
		--i-taxonomy ${params.resultsDir}/assignTaxonomy/taxonomy.qza \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/taxonomyVisualization/taxa-bar-plots-no-Unassigned.qzv
	"""
	else
	"""
		echo "Skipped"	
	"""

}

//diversity analyses
process diversityAnalyses {
	input:
	val 'flag_taxa'
	val 'flag_collapse'
	output:
	script:
    if(params.diversityAnalyses && params.filterTaxa)
	"""
		mkdir -p ${params.resultsDir}/diversityAnalyses
		mkdir -p ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}

		num_samples=\$(echo \$(cat ${params.sampleMetadata} | wc -l) -1 | bc)

		if [ "\$num_samples" -gt 1 ]; then
			qiime diversity core-metrics \
			--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level${params.taxaLevelDiversity}.qza \
			--p-sampling-depth ${params.numReadsDiversity} \
			--m-metadata-file ${params.sampleMetadata} \
			--o-rarefied-table ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/rarefied-table.qza \
			--o-observed-features-vector ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/observed-features-vector.qza \
			--o-shannon-vector ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/shannon-vector.qza \
			--o-evenness-vector ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/evenness-vector.qza \
			--o-jaccard-distance-matrix ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-distance-matrix.qza \
			--o-bray-curtis-distance-matrix ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-distance-matrix.qza \
			--o-jaccard-pcoa-results ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-pcoa-results.qza \
			--o-bray-curtis-pcoa-results ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-pcoa-results.qza \
			--o-jaccard-emperor ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-emperor.qzv \
			--o-bray-curtis-emperor ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-emperor.qzv;
		fi

		for f in \$(find ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity} | grep \"qza\" | grep -v \"distance-matrix\"); do

			fn=\$(echo \$f | sed \'s/\\.qza//\');
			
			qiime metadata tabulate \
			--m-input-file \$f \
			--o-visualization \$fn".qzv";
		done

		qiime diversity alpha-rarefaction \
		--i-table ${params.resultsDir}/filterTaxa/table-${params.taxaOfInterest}-min-freq-${params.minNumReadsTaxaOfInterest}-collapsed-absfreq-level${params.taxaLevelDiversity}.qza \
		--p-max-depth ${params.numReadsDiversity} \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/diversityAnalyses/taxa-${params.taxaOfInterest}-samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/alpha-rarefaction.qzv
	"""
	else if(params.diversityAnalyses && !params.filterTaxa)
	"""
		mkdir -p ${params.resultsDir}/diversityAnalyses
		mkdir -p ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}

		num_samples=\$(echo \$(cat ${params.sampleMetadata} | wc -l) -1 | bc)

		if [ "\$num_samples" -gt 1 ]; then
			qiime diversity core-metrics \
			--i-table ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level${params.taxaLevelDiversity}.qza \
			--p-sampling-depth ${params.numReadsDiversity} \
			--m-metadata-file ${params.sampleMetadata} \
			--o-rarefied-table ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/rarefied-table.qza \
			--o-observed-features-vector ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/observed-features-vector.qza \
			--o-shannon-vector ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/shannon-vector.qza \
			--o-evenness-vector ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/evenness-vector.qza \
			--o-jaccard-distance-matrix ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-distance-matrix.qza \
			--o-bray-curtis-distance-matrix ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-distance-matrix.qza \
			--o-jaccard-pcoa-results ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-pcoa-results.qza \
			--o-bray-curtis-pcoa-results ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-pcoa-results.qza \
			--o-jaccard-emperor ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/jaccard-emperor.qzv \
			--o-bray-curtis-emperor ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/bray-curtis-emperor.qzv;
		fi

		for f in \$(find ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity} | grep \"\\.qza\" | grep -v \"distance-matrix\"); do
			
			fn=\$(echo \$f | sed \'s/\\.qza//\');
			
			qiime metadata tabulate \
			--m-input-file \$f \
			--o-visualization \$fn".qzv";
		done

		qiime diversity alpha-rarefaction \
		--i-table ${params.resultsDir}/collapseTables/table-collapsed-absfreq-level${params.taxaLevelDiversity}.qza \
		--p-max-depth ${params.numReadsDiversity} \
		--m-metadata-file ${params.sampleMetadata} \
		--o-visualization ${params.resultsDir}/diversityAnalyses/samplingDepth-${params.numReadsDiversity}-level${params.taxaLevelDiversity}/alpha-rarefaction.qzv
	"""
	else
	"""
		echo "Skipped"
	"""
}

workflow {
	sequences_db=Channel.fromPath(params.dbSequencesFasta)
	taxonomy_db=Channel.fromPath(params.dbTaxonomyTsv)
	importDb(sequences_db, taxonomy_db)
	concatenateFastq(params.workDir)
	filterFastq(concatenateFastq.out)
	downsampleFastq(filterFastq.out)
	importFastq(downsampleFastq.out)
	derepSeq(importFastq.out)
	assignTaxonomy(derepSeq.out, importDb.out)
	filterTaxa(assignTaxonomy.out)
	taxonomyVisualization(filterTaxa.out)
	collapseTables(assignTaxonomy.out)
	dataQC(importFastq.out)
	diversityAnalyses(filterTaxa.out, collapseTables.out)
}
