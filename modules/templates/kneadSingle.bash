#!/usr/bin/env bash

set -euo pipefail
${params.kneaddataCommand} \
  --mateIds_are_equal ${params.mateIds_are_equal} \
  --query_mate_separator ${params.query_mate_separator} \
  --input ${id}.fastq \
  --output . 
mv ${id}_kneaddata.trimmed.fastq ${id}_1_kneaddata.fastq
