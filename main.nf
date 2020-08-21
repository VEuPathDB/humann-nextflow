#!/usr/bin/env nextflow

sampleToFastqChSingle = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter(it.length == 2)

sampleToFastqChPaired = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter(it.length == 3)

process initWorkingDirSingle {
  input:
  tuple val(sampleId), val(fastq) from sampleToFastqChSingle

  output:
  tuple val(sampleId), path(workingDir) into workingDirsChSingle

  script:
  workingDir = "$baseDir/working/$sampleId"
  
  '''
  mkdir -p "$workingDir/input"
  wget $fastq -O "$workingDir/input/$sampleId.fastq"
  '''
}

process initWorkingDirPaired {
  input:
  tuple val(sampleId), val(fastq1), val(fastq2) from sampleToFastqChPaired

  output:
  tuple val(sampleId), path(workingDir) into workingDirsChPaired

  script:
  workingDir = "$baseDir/working/$sampleId"

  '''
  mkdir -p "$workingDir/input"
  wget $fastq1 -O "$workingDir/input/${sampleId}_1.fastq"
  wget $fastq2 -O "$workingDir/input/${sampleId}_2.fastq"
  '''
}

process kneadWorkingDirSingle {
  input:
  tuple val(sampleId), path(workingDir) from workingDirsChSingle

  output:
  tuple val(sampleId), path(workingDir), path(kneadedFastq) into kneadedFastqs

  script:
    db=refs.kneaddataDb
    trimmomatic = refs.trimmomaticLib
    fastq = "input/${sampleId}.fastq"
    kneadedFastq = "trimmed/${sampleId}.fastq"
  '''

  cd $workingDir \
    && kneaddata -reference-db $db --trimmomatic $trimmomatic --cat-final-output --input $fastq --output trimmed  \
    && test -f $kneadedFastq && rm -r input
  '''
}

process kneadWorkingDirPaired {
  input:
  tuple val(sampleId), path(workingDir) from workingDirsChPaired

  output:
  tuple val(sampleId), path(workingDir), path(kneadedFastq) into kneadedFastqs

  script:
    db=refs.kneaddataDb
    trimmomatic = refs.trimmomaticLib
    fastq1 = "input/${sampleId}_1.fastq"
    fastq2 = "input/${sampleId}_2.fastq"
    kneadedFastq = "trimmed/${sampleId}.fastq"
  '''

  cd $workingDir \
    && kneaddata -reference-db $db --trimmomatic $trimmomatic --cat-final-output --input $fastq1 --input $fastq2 --output trimmed  \
    && test -f $kneadedFastq && rm -r input
  '''
}

process runHumannMain {
  input:
  tuple val(sampleId), path(workingDir), path(kneadedFastq) from kneadedFastqs

  output:
  tuple val(sampleId), path(metaphlanOutPath) into metaphlanBugs

  script:
  metaphlanOutPath = "${sampleId}_humann_temp/${sampleId}_metaphlan_bugs_list.tsv"
  genesOutPath = "${sampleId}_genefamilies.tsv
  pathwaysAbundanceOutPath = "${sampleId}_pathabundance.tsv
  pathwaysCoverageOutPath = "${sampleId}_pathcoverage.tsv
  '''
  humann --input $kneadedFastq --output .
  '''
}

process joinAllMetaphlanBugs {

  publishDir params.resultDir

  input:
  tuple val(sampleId), path(metaphlanOutPath) from metaphlanBugs

  output:
  path(resultNameTaxons) into results

  script:
    taxonsOut = params.resultTaxonsPath
  """
  perl -ne 'my (@bars) = /\|/g; print if @bars == 6 ' < $metaphlanOutPath | while read -r line; do
    echo "$sampleId\t$line"
  done >> $taxonsOut
  """
}

process joinAllHumannPathways {

  publishDir params.resultDir

  input:
  tuple val(sampleId), path(metaphlanOutPath) from metaphlanBugs
}
}
