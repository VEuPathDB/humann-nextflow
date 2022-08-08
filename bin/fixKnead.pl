#!/usr/bin/perl

use strict;
use Getopt::Long;
my $filename;
my $numlen;
my $numbers;
&GetOptions("filename=s" => \$filename);

open(FH, '<', $filename) or die $!;

while(<FH>){
    if ($_ =~ /(^\@\S*\.)(\d+)(l.*)/) {
	$numlen = length($2);
	$numbers = substr($2, 0, $numlen/2);
	print "$1$numbers $numbers $3\n";
    }
    else {
        print "$_";
    }
}

close(FH);



