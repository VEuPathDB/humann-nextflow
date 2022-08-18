#!/usr/bin/env bash

set -euo pipefail
${params.kneaddataCommand} \
  --input ${id}_1.fastq \
  --input ${id}_2.fastq \
  --cat-final-output \
  --output .
