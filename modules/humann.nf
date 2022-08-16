#!/usr/bin/env nextflow
import nextflow.splitter.CsvSplitter
nextflow.enable.dsl=2


process downloadFiles {
  input:
    val id

  output:
    tuple val(id), path("${id}**.fastq")

  script:
    template 'downloadFiles.bash'
}


process knead {
  label 'download_and_preprocess'

  input:
    tuple val(id), path(readsFastq)

  output:
    tuple val(id), path("*_*_kneaddata.fastq")

  script:
  if(params.libraryLayout == 'single')
    template 'kneadSingle.bash'
  else if(params.libraryLayout == 'paired')
    template 'kneadPaired.bash'
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
    template 'runHumann.bash'
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
    template 'groupFunctionalUnits.bash'
}


process aggregateFunctionAbundances {
  publishDir params.resultDir, mode: 'move'

  input:
    file('*') 
    each (groupType) 

  output:
    file("${groupType}s.tsv")

  script:
    template 'aggregateFunctionAbundances.bash'
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
    template 'aggregatePathwayAbundances.bash'
}

process aggregatePathwayCoverages {
  publishDir params.resultDir, mode: 'move'

  input:
    file('*') 

  output:
    file("pathway_coverages.tsv")

  script:
    template 'aggregatePathwayCoverages.bash'
}


workflow humann {
  take:
    input

  main: 

    functionalUnitNames = [
    eggnog: "eggnog",
     go: "go",
    ko: "kegg-orthology",
    level4ec: "ec",
    pfam: "pfam",
    rxn: "metacyc-rxn",
    ]

    if (params.downloadMethod == 'sra') {
      ids = Channel.fromList(input)
      files = downloadFiles(ids)
    }
  
    kneadedReads = knead(input)
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