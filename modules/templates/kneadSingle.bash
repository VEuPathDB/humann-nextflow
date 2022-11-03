#!/usr/bin/env bash

set -euo pipefail
${params.kneaddataCommand} \
  -i1 ${id}.fastq \
  --output . 
mv ${id}_kneaddata.trimmed.fastq ${id}_1_kneaddata.fastq
