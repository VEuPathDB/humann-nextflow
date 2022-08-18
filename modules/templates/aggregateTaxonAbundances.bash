#!/usr/bin/env bash

set -euo pipefail
merge_metaphlan_tables.py *.metaphlan.out \
   | grep -v '^#' \
   | cut -f 1,3- \
   | perl -pe 'if($.==1){s/.metaphlan//g}' \
   > taxon_abundances.tsv
