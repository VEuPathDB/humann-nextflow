#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($strain,$idList);

# Creating Arguments
&GetOptions("strain|i=s" => \$strain,
            "idList|l=s" => \$idList);

my @idArray = split(/,/, $idList);

foreach my $id (@idArray){

    system("fasterq-dump --split-3 ${id}");
    if(-e "${id}_1.fastq"){
        system("cat $id\_1.fastq >> ${strain}_1.fastq");
    }
    if(-e "${id}.fastq"){
        system("cat $id.fastq >> ${strain}_1.fastq");
    }
    if(-e "${id}_2.fastq"){
	system("cat $id\_2.fastq >> ${strain}_2.fastq");
    }
    
}

