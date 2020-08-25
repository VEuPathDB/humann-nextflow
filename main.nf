sampleToFastqLocationsSingle = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 2}
process prepareReadsSingle {
  
  input:
  tuple val(sample), val(fastq) from sampleToFastqLocationsSingle

  output:
  tuple val(sample), file("reads_kneaddata.fastq") into kneadedReadsSingle

  script:
    """
    ${params.wgetCommand} $fastq -O reads.fastq.gz
    gunzip *gz
    kneaddata ${params.kneaddataFlags} \
      --input reads.fastq \
      --output .
    """
}

sampleToFastqLocationsPaired = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 3}
process prepareReadsPaired {
  
  input:
  tuple val(sample), val(fastq1), val(fastq2) from sampleToFastqLocationsPaired

  output:
  tuple val(sample), file("reads_kneaddata.fastq") into kneadedReadsPaired

  script:
    """
    ${params.wgetCommand} $fastq1 -O reads.fastq.gz
    ${params.wgetCommand} $fastq2 -O reads_R.fastq.gz
    gunzip *gz
    kneaddata ${params.kneaddataFlags} \
      --input reads.fastq --input reads_R.fastq --cat-final-output \
      --output . 
    """
}

kneadedReads = kneadedReadsSingle.mix(kneadedReadsPaired)


process runHumann {
  label 'mem_4c'

  input:
  tuple val(sample), file(kneadedReads) from kneadedReads

  output:
  file("${sample}.bugs.tsv") into out_bugs
  file("${sample}.ec_named.tsv") into out_ecs
  file("${sample}.pathway_abundance.tsv") into out_pas
  file("${sample}.pathway_coverage.tsv") into out_pcs

  script:
  """
  mv -v $kneadedReads reads.fastq
  ${params.humannCommand} --threads 4 --input reads.fastq --output .

  mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv ${sample}.bugs.tsv
  
  humann_renorm_table --input reads_genefamilies.tsv --output reads_genefamilies_cpm.tsv --units cpm --update-snames
  humann_regroup_table --input reads_genefamilies_cpm.tsv --output reads_ec4.tsv --groups uniref50_level4ec
  humann_rename_table --input reads_ec4.tsv --output ${sample}.ec_named.tsv --names ec

  mv -v reads_pathabundance.tsv ${sample}.pathway_abundance.tsv
  mv -v reads_pathcoverage.tsv ${sample}.pathway_coverage.tsv

  """
}

process aggregate_bugs {
  publishDir params.resultDir

  input:
  file('*.bugs.tsv') from out_bugs.collect()


  output:
  file("all_bugs.tsv")

  script:
  """
  grep "" *.bugs.tsv > all_bugs.tsv
  """
}

process aggregate_ecs {
  publishDir params.resultDir

  input:
  file('*.ec_named.tsv') from out_ecs.collect()

  output:
  file("all_ecs.tsv")

  script:
  """
  humann_join_tables --input . --output all_ecs.tsv --file_name ec_named
  """
}

process aggregate_pathway_abundances {
  publishDir params.resultDir

  input:
  file('*.pathway_abundance.tsv') from out_pas.collect()

  output:
  file("all_pathway_abundances.tsv")

  script:
  """
  humann_join_tables --input . --output all_pathway_abundances.tsv --file_name pathway_abundance
  """
}

process aggregate_pathway_coverages {
  publishDir params.resultDir

  input:
  file('*.pathway_coverage.tsv') from out_pcs.collect()

  output:
  file("all_pathway_coverages.tsv")

  script:
  """
  humann_join_tables --input . --output all_pathway_coverages.tsv --file_name pathway_coverage
  """
}

