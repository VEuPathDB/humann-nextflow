# MicrobiomeDB WGS workflow

## Install

pip3, metaphlan, and humann:
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
/home/wbazant/.local/bin/pip3 install {cython, numpy, metaphlan, humann, kneaddata}

diamond:
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.2/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

minpath:
wget 'https://omics.informatics.indiana.edu/mg/get.php?justdoit=yes&software=minpath1.4.tar.gz' -O minpath1.4.tar.gz

trimmomatic:
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

databases:
kneaddata_database --download human_genome bowtie2 ~/kneaddata_databases/

metaphlan_database 2019's chocophlan
humann_database uniref50  # the gene-to-name reference is based on uniref50 so use this one



