sampleToFastqLocationsSingle = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 2}
process prepareReadsSingle {
  
  maxForks 5

  input:
  tuple val(sample), val(fastq) from sampleToFastqLocationsSingle

  output:
  tuple val(sample), file("reads_kneaddata.fastq") into kneadedReadsSingle

  script:
  """
  ${params.wgetCommand} $fastq -O reads.fastq.gz
  gunzip *gz
  ${params.kneaddataCommand} \
    --input reads.fastq \
    --output .
  """
  afterScript:
  """
  rm -v reads.fastq
  """
}

sampleToFastqLocationsPaired = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 3}
process prepareReadsPaired {
  
  maxForks 5

  input:
  tuple val(sample), val(fastq1), val(fastq2) from sampleToFastqLocationsPaired

  output:
  tuple val(sample), file("reads_kneaddata.fastq") into kneadedReadsPaired

  script:
  """
  ${params.wgetCommand} $fastq1 -O reads.fastq.gz
  ${params.wgetCommand} $fastq2 -O reads_R.fastq.gz
  gunzip *gz
  ${params.kneaddataCommand} \
    --input reads.fastq --input reads_R.fastq --cat-final-output \
    --outputnput name is the same as the channel name, the from part of the input declaration can be omitted. Thus, the above example could be written as  . 
  """

  afterScript:
  """
  rm -v reads.fastq reads_R.fastq.gz
  """
}

kneadedReads = kneadedReadsSingle.mix(kneadedReadsPaired)


process runHumann {
  label 'mem_4c'

  input:
  tuple val(sample), file(kneadedReads) from kneadedReads

  output:
  file("${sample}.bugs.tsv") into outChTaxonAbundances
  file("${sample}.ec_named.tsv") into outChEnzymeAbundances
  file("${sample}.pathway_abundance.tsv") into outChPathwayAbundances
  file("${sample}.pathway_coverage.tsv") into outChPathwayCoverages

  script:
  """
  mv -v $kneadedReads reads.fastq
  ${params.humannCommand} --threads 4 --input reads.fastq --output .

  mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv ${sample}.bugs.tsv
  mv -v reads_humann_temp/reads.log humann.log
  
  humann_renorm_table --input reads_genefamilies.tsv --output reads_genefamilies_cpm.tsv --units cpm --update-snames

  humann_regroup_table --input reads_genefamilies_cpm.tsv --output reads_ec4.tsv --groups uniref50_rxn
  humann_rename_table --input reads_ec4.tsv --output ${sample}.ec_named.tsv --names metacyc-rxn

  mv -v reads_pathabundance.tsv ${sample}.pathway_abundance.tsv
  mv -v reads_pathcoverage.tsv ${sample}.pathway_coverage.tsv
  """

  afterScript:
  """
  rm -rv reads_humann_temp
  """
}

process aggregateTaxonAbundances {
  publishDir params.resultDir

  input:
  file('*.bugs.tsv') from outChTaxonAbundances.collect()


  output:
  file("all_bugs.tsv")

  script:
  """
  grep "" *.bugs.tsv > all_bugs.tsv
  """
}

process aggregateEnzymeAbundances {
  publishDir params.resultDir

  input:
  file('*.ec_named.tsv') from outChEnzymeAbundances.collect()

  output:
  file("all_ecs.tsv")

  script:
  """
  humann_join_tables --input . --output all_ecs.tsv --file_name ec_named
  """
}

process aggregatePathwayAbundances {
  publishDir params.resultDir

  input:
  file('*.pathway_abundance.tsv') from outChPathwayAbundances.collect()

  output:
  file("all_pathway_abundances.tsv")

  script:
  """
  humann_join_tables --input . --output all_pathway_abundances.tsv --file_name pathway_abundance
  """
}

process aggregatePathwayCoverages {
  publishDir params.resultDir

  input:
  file('*.pathway_coverage.tsv') from outChPathwayCoverages.collect()

  output:
  file("all_pathway_coverages.tsv")

  script:
  """
  humann_join_tables --input . --output all_pathway_coverages.tsv --file_name pathway_coverage
  """
}

