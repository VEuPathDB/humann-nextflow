sampleToFastqLocationsSingle = Channel
  .fromPath(params.sampleToFastqsPath)
  .splitCsv(sep: "\t")
  .filter {it.size() == 2}
process prepareReadsSingle {
  
  input:
  tuple val(sample), val(fastq) from sampleToFastqLocationsSingle

  output:
  file("reads_kneaddata.fastq") into kneadedReadsSingle

  script:
    """
    wget $fastq -O reads.fastq
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
  file("reads_kneaddata.fastq") into kneadedReadsPaired

  script:
    """
    wget $fastq1 -O reads_1.fastq
    wget $fastq2 -O reads_2.fastq
    kneaddata ${params.kneaddataFlags} \
      --input reads_1.fastq --input reads_2.fastq --cat-final-output \
      --output .
    """
}

process runHumann {
  input:
  file(kneadedReads) from kneadedReadsSingle.join(kneadedReadsPaired)

  output:
  file(kneadedReads) into stdout
  script:
  """
  echo $kneadedReads
  """
}
