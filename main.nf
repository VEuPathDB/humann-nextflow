sampleToFastqLocationsSingle = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 2}
process prepareReadsSingle {
  errorStrategy { 
    sleep(Math.pow(2, task.attempt) * 500 as long);
    if (task.exitStatus == 8 || task.attempt < 4 ) {
      return 'retry' 
    } else {
      return 'terminate' 
    }
  }
  maxRetries 10
  maxForks 5

  afterScript 'rm -v reads.fastq'

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

}

sampleToFastqLocationsPaired = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 3}
process prepareReadsPaired {
  errorStrategy { 
    sleep(Math.pow(2, task.attempt) * 500 as long);
    if (task.exitStatus == 8 || task.attempt < 4 ) {
      return 'retry' 
    } else {
      return 'terminate' 
    }
  }
  maxRetries 10
  maxForks 5

  afterScript 'rm -v reads.fastq reads_R.fastq.gz'

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
    --output .
  """

}

kneadedReads = kneadedReadsSingle.mix(kneadedReadsPaired)


process runHumann {
  label 'mem_4c'

  afterScript 'mv -v reads_humann_temp/reads.log humann.log; test -f reads_humann_temp/reads_metaphlan_bugs_list.tsv && mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv bugs.tsv ; rm -rv reads_humann_temp'

  input:
  tuple val(sample), file(kneadedReads) from kneadedReads

  output:
  file("${sample}.metaphlan.out") into taxonAbundances
  tuple val(sample), file("${sample}.gene_abundance.tsv") into geneAbundances
  file("${sample}.pathway_abundance.tsv") into pathwayAbundances
  file("${sample}.pathway_coverage.tsv") into pathwayCoverages

  script:
  """
  ln -vs $kneadedReads reads.fastq
  ${params.humannCommand} --threads 4 --input reads.fastq --output .

  mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv ${sample}.metaphlan.out
  
  humann_renorm_table --input reads_genefamilies.tsv --output ${sample}.gene_abundance.tsv --units cpm --update-snames

  mv -v reads_pathabundance.tsv ${sample}.pathway_abundance.tsv
  mv -v reads_pathcoverage.tsv ${sample}.pathway_coverage.tsv
  """

}

functionalUnitNames = [
  eggnog: "eggnog",
  go: "go",
  ko: "kegg-orthology",
  level4ec: "ec",
  pfam: "pfam",
  rxn: "metacyc-rxn",
]

process groupFunctionalUnits {

  input:
  tuple val(sample), file(geneAbundances) from geneAbundances
  each (groupType) from params.functionalUnits

  output:
  file("${sample}.${groupType}.tsv") into functionAbundances

  script:
  group=params.unirefXX + "_" + groupType
  groupName=functionalUnitNames[groupType]
  """
  humann_regroup_table --input $geneAbundances --output ${group}.tsv --groups ${group}
  humann_rename_table --input ${group}.tsv --output ${sample}.${groupType}.tsv --names ${groupName}
  """
}

process aggregateFunctionAbundances {
  publishDir params.resultDir, mode: 'move'

  input:
  file('*') from functionAbundances.collect()
  each (groupType) from params.functionalUnits

  output:
  file("${groupType}s.tsv")

  script:
  """
  joinTablesForGroupType ${groupType}s.tsv $groupType
  """
}

process aggregateTaxonAbundances {
  publishDir params.resultDir, mode: 'move'

  input:
  file('*') from taxonAbundances.collect()


  output:
  file("taxon_abundances.tsv")

  script:
  '''
  merge_metaphlan_tables.py *.metaphlan.out \
   | grep -v '^#' \
   | cut -f 1,3- \
   | perl -pe 'if($.==1){s/.metaphlan//g}' \
   > taxon_abundances.tsv
  '''
}


process aggregatePathwayAbundances {
  publishDir params.resultDir, mode: 'move'

  input:
  file('*') from pathwayAbundances.collect()

  output:
  file("pathway_abundances.tsv")

  script:
  """
  joinTablesForGroupType pathway_abundances.tsv pathway_abundance
  """
}

process aggregatePathwayCoverages {
  publishDir params.resultDir, mode: 'move'

  input:
  file('*') from pathwayCoverages.collect()

  output:
  file("pathway_coverages.tsv")

  script:
  """
  joinTablesForGroupType pathway_coverages.tsv pathway_coverage
  """
}

