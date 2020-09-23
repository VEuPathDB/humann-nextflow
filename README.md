# Humann(+wget, kneaddata) in a Nextflow pipeline

Runs humann for a list of samples.

In:
- list of sample + fastq URLs, single or paired end
- config

Out:
- species tables from metaphlan
- abundances of ECs and other functional units
- abundances of pathways
- coverages of pathways

Stuff gets nicely formatted and provided as single TSV per result type.

## Install

pip3, metaphlan, and humann:
```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
~/.local/bin/pip3 install {cython, numpy, metaphlan, humann, kneaddata}
```

diamond:
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.2/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```
minpath:
```
cd ~/lib
wget 'https://omics.informatics.indiana.edu/mg/get.php?justdoit=yes&software=minpath1.4.tar.gz' -O minpath1.4.tar.gz

```

trimmomatic:
```
cd ~/lib
wget 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip'
unzip *zip

```

## Configuration
### Reference databases
HUMAnN and Kneaddata need reference databases:

```
kneaddata_database --download human_genome bowtie2 ~/kneaddata_databases/

humann_databases chocophlan full ~/humann_databases
humann_databases uniref uniref90_diamond ~/humann_databases
humann_databases utility_mapping full ~/humann_databases
```

MetaPHlAn will come preconfigured - its installation comes with its own ChocoPhlAn.

### Choosing the reference protein set
`uniref50_diamond` might be better for you instead of `uniref90_diamond`.

In that case, modify or override `nextflow.config` as follows:
```
params {
  unirefXX = 'uniref50'
}
```

If you use `uniref90_ec_filtered_diamond`, you should probably not aggregate into other functional units than ECs, as the results will be biased:
```
params {
  functionalUnits = ["level4ec"]
}
```

### Arguments to tools
You can provide custom installation paths or extra parameters as needed:
```
params {
  kneaddata --trimmomatic ~/lib/Trimmomatic-0.39
  humannCommand = "humann --memory-use maximum"
}
```

### Cluster submission
The most resource-intensive part of this pipeline is the job that runs HumAnN. It is labelled as `mem_4c`. To run it with 10GB (enough for the EC filtered UniRef90):
```
process {
  executor = 'lsf'
  maxForks = 20
 
  withLabel: 'mem_4c' {
   clusterOptions = '-n 4 -M 10000 -R "rusage [mem=10000] span[hosts=1]"'
  }
}

```
## How to provide input
File with two or three columns in the TSV format:

```
sample	fastq URL [second fastq URL]
```
for single or paired read files. The fastqs will be fetched using `wget`.

The pipeline currently doesn't support providing filesystem paths but it wouldn't be hard to modify it!


## Output format
There's a file for taxa, pathway abundances, pathway coverages, and each functional unit specified in the config.

Each file is a TSV. The header contains the row type, followed by sample names. This is the same as the output of `humann_join_tables` script.

Taxa are in the same format, created by adjusting the usual Metaphlan output by removing the frontmatter and the NCBI id column.
