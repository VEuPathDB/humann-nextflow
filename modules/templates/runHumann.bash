#!/usr/bin/env bash

set -euo pipefail
ln -vs -f $kneadedReads unfixed.fastq
perl /usr/local/bin/fixKnead.pl -filename unfixed.fastq > reads.fastq
${params.humannCommand} --input reads.fastq --output .
mv -v reads_humann_temp/reads_metaphlan_bugs_list.tsv ${sample}.metaphlan.out
humann_renorm_table --input reads_genefamilies.tsv --output ${sample}.gene_abundance.tsv --units cpm --update-snames
mv -v reads_pathabundance.tsv ${sample}.pathway_abundance.tsv
mv -v reads_pathcoverage.tsv ${sample}.pathway_coverage.tsv
