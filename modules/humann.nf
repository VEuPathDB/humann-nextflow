#!/usr/bin/env nextflow
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
    if(params.libraryLayout.toLowerCase() == 'single')
      template 'kneadSingle.bash'
    else if(params.libraryLayout.toLowerCase() == 'paired')
      template 'kneadPaired.bash'
}

process runHumann {
  afterScript 'mv -v reads_humann_temp/reads.log humann.log; test -f reads_humann_temp/reads_metaphlan_bugs_list.tsv && mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv bugs.tsv ; rm -rv reads_humann_temp'

  input:
    tuple val(sample), path(kneadedReads)

  output:
    path("${sample}.metaphlan.out"), emit: metaphlan_output
    tuple val(sample), file("${sample}.gene_abundance.tsv"), emit: sample_geneabundance_tuple
    path("${sample}.pathway_abundance.tsv"), emit: pathway_abundance
    path("${sample}.pathway_coverage.tsv"), emit: pathway_coverage

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

    if (params.downloadMethod.toLowerCase() == 'sra') {
      ids = Channel.fromList(input)
      files = downloadFiles(ids)
      kneadedReads = knead(files)
    }
    else {
      kneadedReads = knead(input)
    }
    
    humannOutput = runHumann(kneadedReads)

    taxonAbundances = humannOutput.metaphlan_output.collect()
    aggregateTaxonAbundances(taxonAbundances)

    functionAbundances = groupFunctionalUnits(humannOutput.sample_geneabundance_tuple, params.functionalUnits, functionalUnitNames )
    functionAbundancesCollected = functionAbundances.collect()
    aggregateFunctionAbundances(functionAbundancesCollected, params.functionalUnits)

    pathwayAbundances = humannOutput.pathway_abundance.collect()
    aggregatePathwayAbundances(pathwayAbundances)

    pathwayCoverages = humannOutput.pathway_coverage.collect()
    aggregatePathwayCoverages(pathwayCoverages)
}