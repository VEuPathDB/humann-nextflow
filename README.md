# humann-nextflow

Nextflow pipeline that runs HUMAnN functional profiling (plus KneadData preprocessing and MetaPhlAn taxonomic profiling) over metagenomic/metatranscriptomic short-read samples for VEuPathDB microbiome studies.

## Overview

For each sample, humann-nextflow optionally downloads raw reads from SRA, quality-trims and decontaminates them with `kneaddata`, then runs `humann` to jointly profile microbial community taxonomy (via MetaPhlAn) and functional pathway/gene-family abundance. Per-sample results are grouped into configurable functional units (e.g. EC numbers, KEGG orthologs, Pfam, GO, eggNOG) and merged across all samples into single, EDA/MicroBiomeDB-ready abundance and coverage tables. It is the pipeline VEuPathDB uses to turn raw sequencing reads for a microbiome dataset into the taxon and function abundance tables loaded downstream by `eda-nextflow`.

## Requirements

- Nextflow
- Docker (the pipeline runs every process in the `veupathdb/humann` container image)
- Local HUMAnN/KneadData/MetaPhlAn reference databases, bind-mounted into the container (see Reference databases below)
- SRA Toolkit access if downloading reads by accession (`fasterq-dump`, bundled in the container)

### Reference databases

HUMAnN, KneadData, and MetaPhlAn require reference databases to be downloaded separately and mounted into the container:

- UniRef90: `http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz`
- UniRef50: `http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901.tar.gz`
- ChocoPhlAn: `http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v31_CHOCOPhlAn_201901.tar`
- Utility mapping: `http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901.tar.gz`
- KneadData contaminant reference: `http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz`
- MetaPhlAn database: `http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v31_CHOCOPhlAn_201901.tar`

`nextflow.config`'s `docker.runOptions` binds these into the container at `/humann_databases`, `/kneaddata_databases`, and the MetaPhlAn install's `metaphlan_databases` directory; adjust the host paths (or override `runOptions`) to point at your database locations.

## Usage

```
nextflow run VEuPathDB/humann-nextflow -r main \
  --downloadMethod local \
  --inputPath /path/to/fastqs \
  --libraryLayout paired \
  --resultDir /path/to/results \
  -resume -C site.config
```

To pull reads from SRA instead of local FASTQs, set `--downloadMethod sra` and point `--inputPath` at a TSV of run accessions (header `run_accession`, see `data/sample-to-fastqs.tsv`).

The pipeline has a single entry point (the default, unnamed `workflow`); there are no named `-entry` targets.

## Key parameters

- `downloadMethod` — `local` to read FASTQs from `inputPath`, or `sra` to download runs listed in a TSV at `inputPath` via `fasterq-dump`
- `inputPath` — directory of local FASTQ files (`local`) or a TSV of `run_accession`s (`sra`)
- `libraryLayout` — `paired` or `single`; determines FASTQ pairing and the KneadData invocation used
- `resultDir` — output directory for the final merged abundance/coverage tables
- `kneaddataCommand` — full `kneaddata` command line, including reference database paths and trimming options
- `humannCommand` — full `humann` command line, including any DIAMOND search options
- `unirefXX` — which UniRef build (`uniref90` or `uniref50`) HUMAnN's gene family output is grouped against
- `functionalUnits` — list of functional unit types to aggregate per sample (e.g. `level4ec`, `ko`, `pfam`, `go`, `eggnog`, `rxn`)
- `mateIds_are_equal` / `query_mate_separator` — control how KneadData matches paired-end mate IDs; defaults differ automatically between `sra` and `local` download methods

### Cluster submission

The most resource-intensive step is `runHumann`, labelled `mem_4c`. To size and retry it on an LSF cluster:

```
process {
  executor = 'lsf'
  maxForks = 20

  withLabel: 'mem_4c' {
    errorStrategy = { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries = 3
    clusterOptions = { task.attempt == 1 ?
      '-n 4 -M 12000 -R "rusage [mem=12000] span[hosts=1]"'
      : task.attempt == 2 ?
      '-n 4 -M 17000 -R "rusage [mem=17000] span[hosts=1]"'
      : '-n 4 -M 25000 -R "rusage [mem=25000] span[hosts=1]"'
    }
  }
}
```

Memory needs scale with reference database size and input size, so job memory is typically tuned per site.

## Output

Written to `resultDir`, one merged TSV per result type across all samples (equivalent to `humann_join_tables` output, with a header row of sample names):

- `taxon_abundances.tsv` — MetaPhlAn species/clade abundance table
- `<functionalUnit>s.tsv` for each entry in `functionalUnits` (e.g. `level4ecs.tsv`, `kos.tsv`) — grouped functional unit abundances, normalized to counts per million
- `pathway_abundances.tsv` — MetaCyc pathway abundances
- `pathway_coverages.tsv` — MetaCyc pathway coverages
