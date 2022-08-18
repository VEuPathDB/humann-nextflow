#!/usr/bin/env bash

set -euo pipefail
humann_regroup_table --input $geneAbundances --output ${group}.tsv --groups ${group}
humann_rename_table --input ${group}.tsv --output ${sample}.${groupType}.tsv --names ${groupName}
