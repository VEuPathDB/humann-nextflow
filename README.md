# Humann(+wget, kneaddata) in a Nextflow pipeline

Runs humann for a list of samples.

In:
- list of sample + fastq URLs
- configuration

Out:
- species tables from metaphlan
- abundances of ECs and other functional units
- abundances of pathways
- coverages of pathways

## Overview
The fastqs are downloaded locally, trimmed ("kneaded") with `kneaddata`, and - if paired - merged together.
Then `humann` is ran on each input, with specified CPU and memory.
The results are merged and returned as a single file per result type.

![diagram](https://raw.githubusercontent.com/wbazant/humann-nextflow/master/flowchart.svg)



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

### Reference databases
HUMAnN and Kneaddata need reference databases.

```
kneaddata_database --download human_genome bowtie2 ~/kneaddata_databases/

humann_databases chocophlan full ~/humann_databases
humann_databases uniref uniref90_diamond ~/humann_databases
humann_databases utility_mapping full ~/humann_databases
```

MetaPHlAn will come preconfigured - its installation comes with its own ChocoPhlAn.
## Configuration guide

### How to provide input
Prepare a file with two or three columns in the TSV format:

```
sample	fastq URL [second fastq URL]
```
for single or paired read files. The fastqs will be fetched using `wget`.

and place it under the path `data/sample-to-fastqs.tsv` or modify the config.

The pipeline currently doesn't support providing filesystem paths but it wouldn't be hard to modify it!


### Choosing the reference protein set
`uniref50_diamond` might be more economical than `uniref90_diamond`.

In that case, modify or override `nextflow.config` as follows:
```
params {
  unirefXX = 'uniref50'
}
```

If you only want the pathway results, you can use `uniref90_ec_filtered_diamond` or `uniref50_ec_filtered_diamond`.

In that case, you should probably not aggregate into other functional units than ECs, as the results will be biased:
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

Memory use for each job depends on the reference size, and input size. To retry failed jobs with more memory, you could modify the config as follows:

```
  withLabel: 'mem_4c' {
    errorStrategy = { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries = 3
    clusterOptions = { task.attempt == 1 ?
      '-n 4 -M 12000 -R \"rusage [mem=12000] span[hosts=1]\"'
      : task.attempt == 2 ?
      '-n 4 -M 17000 -R \"rusage [mem=17000] span[hosts=1]\"'
      : '-n 4 -M 25000 -R \"rusage [mem=25000] span[hosts=1]\"'
    }
  }

```


## Output format
There's a file for taxa, pathway abundances, pathway coverages, and each functional unit specified in the config.

Each file is a TSV. The header contains the row type, followed by sample names. This is the same as the output of `humann_join_tables` script.

Taxa are in the same format, created by adjusting the usual Metaphlan output by removing the frontmatter and the NCBI id column.
