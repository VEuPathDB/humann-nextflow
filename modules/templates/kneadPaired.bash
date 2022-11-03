#!/usr/bin/env bash

set -euo pipefail
${params.kneaddataCommand} \
  -i1 ${id}_1.fastq \
  -i2 ${id}_2.fastq \
  --cat-final-output \
  -o .
