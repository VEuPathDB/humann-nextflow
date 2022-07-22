import nextflow.splitter.CsvSplitter

nextflow.enable.dsl=2

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

process prepSra {
  label 'prep'
  input:
  tuple val(sample), path(runAccession)
  output:
  tuple val(sample), path("*_*_kneaddata.fastq")
  script:
  if(params.libraryLayout == 'single')
  """
  gzip -d --force *.fastq.gz
  ${params.kneaddataCommand} \
    --input ${sample}.fastq \
    --output . \
    -db $params.knead_databases
  """
  else if(params.libraryLayout == 'paired')
  """
  gzip -d --force *.fastq.gz
  ${params.kneaddataCommand} \
    --input ${sample}_1.fastq \
    --input ${sample}_2.fastq \
    --cat-final-output \
    --output . \
    -db $params.knead_databases
  """
}

process runHumann {
  afterScript 'mv -v reads_humann_temp/reads.log humann.log; test -f reads_humann_temp/reads_metaphlan_bugs_list.tsv && mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv bugs.tsv ; rm -rv reads_humann_temp'
  input:
  tuple val(sample), path(kneadedReads)
  output:
  file("${sample}.metaphlan.out")
  tuple val(sample), file("${sample}.gene_abundance.tsv")
  file("${sample}.pathway_abundance.tsv")
  file("${sample}.pathway_coverage.tsv")

  script:
  """
  humann_config --update database_folders nucleotide "$params.humann_databases/chocophlan"
  humann_config --update database_folders protein "$params.humann_databases/uniref"
  humann_config --update database_folders utility_mapping "$params.humann_databases/utility_mapping" 
  ln -vs -f $kneadedReads reads.fastq
  ${params.humannCommand} --input reads.fastq --output .
  mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv ${sample}.metaphlan.out
  humann_renorm_table --input reads_genefamilies.tsv --output ${sample}.gene_abundance.tsv --units cpm --update-snames
  mv -v reads_pathabundance.tsv ${sample}.pathway_abundance.tsv
  mv -v reads_pathcoverage.tsv ${sample}.pathway_coverage.tsv
  """
}

process groupFunctionalUnits {
  input:
  tuple val(sample), file(geneAbundances)
  each (groupType)
  val functionalUnitNames
  output:
  file("${sample}.${groupType}.tsv") 
  script:
  group=params.unirefXX + "_" + groupType
  groupName=functionalUnitNames[groupType]
  """
  humann_config --update database_folders utility_mapping "$params.humann_databases/utility_mapping" 
  humann_regroup_table --input $geneAbundances --output ${group}.tsv --groups ${group}
  humann_rename_table --input ${group}.tsv --output ${sample}.${groupType}.tsv --names ${groupName}
  """
}

process aggregateFunctionAbundances {
  publishDir params.resultDir, mode: 'move'

  input:
  file('*') 
  each (groupType) 

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
  file('*') 


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
  file('*') 

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
  file('*') 

  output:
  file("pathway_coverages.tsv")

  script:
  """
  joinTablesForGroupType pathway_coverages.tsv pathway_coverage
  """
}

workflow {

  functionalUnitNames = [
  eggnog: "eggnog",
  go: "go",
  ko: "kegg-orthology",
  level4ec: "ec",
  pfam: "pfam",
  rxn: "metacyc-rxn",
  ]
  
  accessions = fetchRunAccessions(params.inputPath)
  input = Channel.fromSRA(accessions, apiKey: params.apiKey, protocol: "http")
  kneadedReads = prepSra(input)
  
  humannOutput = runHumann(kneadedReads)

  taxonAbundances = humannOutput[0].collect()
  aggregateTaxonAbundances(taxonAbundances)

  functionAbundances = groupFunctionalUnits(humannOutput[1], params.functionalUnits, functionalUnitNames )
  functionAbundancesCollected = functionAbundances.collect()
  aggregateFunctionAbundances(functionAbundancesCollected, params.functionalUnits)

  pathwayAbundances = humannOutput[2].collect()
  aggregatePathwayAbundances(pathwayAbundances)

  pathwayCoverages = humannOutput[3].collect()
  aggregatePathwayCoverages(pathwayCoverages)
}