#!/usr/bin/env bash

set -euo pipefail
${params.kneaddataCommand} \
  --input ${id}_1.fastq \
  --output . 
