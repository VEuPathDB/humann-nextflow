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

if(params.downloadMethod == 'sra') {
  input = fetchRunAccessions(params.inputPath)
}
else if (params.downloadMethod == 'local') {
  input = Channel.fromFilePairs(params.localFileLocation + "/*_{1,2}.fastq")
}
else {
  throw new Exception("Non-valid value for params.downloadMethod")
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