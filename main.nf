#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def fetchRunAccessions( tsv ) {
    def splitter = new nextflow.splitter.CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    return splitter.target( reader ).list().collect { row -> row['run_accession'] }
}


//---------------------------------------------------------------
// Includes
//---------------------------------------------------------------

include { humann } from './modules/humann.nf'

//---------------------------------------------------------------
// Main Workflow
//---------------------------------------------------------------


workflow {
  if(params.downloadMethod.toLowerCase() == 'sra') {
    input = fetchRunAccessions(params.inputPath)
  }
  else if (params.downloadMethod.toLowerCase() == 'local') {
    if (params.libraryLayout.toLowerCase() == 'paired') {
      input = channel.fromFilePairs(params.inputPath + "/*_{1,2}.fastq")
    }
    else if (params.libraryLayout.toLowerCase() == 'single') {
      input = channel.fromPath(params.inputPath + "/*.fastq").map { file -> tuple(file.baseName, [file]) }
    }
    else {
      error("Non-valid value for params.libraryLayout")
    }
  }
  else {
    error("Non-valid value for params.downloadMethod")
  }

  if(!params.mateIds_are_equal) {
    if(params.downloadMethod.toLowerCase() == 'sra') {
      params.mateIds_are_equal = 'True'
    }
    else if (params.downloadMethod.toLowerCase() == 'local') {
      params.mateIds_are_equal = 'False'
    }
  }

  // Setting Defaults If Parameters Do Not Exist

  if(!params.query_mate_separator) {
    if(params.downloadMethod.toLowerCase() == 'sra') {
      params.query_mate_separator = "."
    }
    else if (params.downloadMethod.toLowerCase() == 'local') {
      params.query_mate_separator = "/"
    }
  }

  humann(input)

}