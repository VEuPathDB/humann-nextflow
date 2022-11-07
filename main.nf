#!/usr/bin/env nextflow
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


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

if(params.downloadMethod.toLowerCase() == 'sra') {
  input = fetchRunAccessions(params.inputPath)
}
else if (params.downloadMethod.toLowerCase() == 'local') {
  if (params.libraryLayout.toLowerCase() == 'paired') {
    input = Channel.fromFilePairs(params.inputPath + "/*_{1,2}.fastq")
  }
  else if (params.libraryLayout.toLowerCase() == 'single') {
    input = Channel.fromPath(params.inputPath + "/*.fastq").map { file -> tuple(file.baseName, [file]) }
  }
  else {
    throw new Exception("Non-valid value for params.libraryLayout")
  }
}
else {
  throw new Exception("Non-valid value for params.downloadMethod")
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

//---------------------------------------------------------------
// Includes
//---------------------------------------------------------------

include { humann } from './modules/humann.nf'

//---------------------------------------------------------------
// Main Workflow
//---------------------------------------------------------------


workflow {

  humann(input)

}