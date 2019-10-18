#!/usr/bin/perl -w

# Copyright (c)   AB_Life 2013
# Writer:         chengchao
# Program Date:   2013.

use strict;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme£¬you must write the detailed time¡¢discriptions¡¢parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

if (scalar(@ARGV)<1) {
	print"Usage:$0 samfile fafile outdir\n";
	exit;
}

my $in = $ARGV[0];
my $fa = $ARGV[1];
my $outdir = $ARGV[2];

my $bam_tmp = $in;
$bam_tmp=~s/\S+\///;
$bam_tmp=~s/sam/bam.tmp/;
my $bam = $in;
$bam=~s/\S+\///;
$bam=~s/\.sam//;
`cd $outdir && samtools view -bT $fa $in -o $bam_tmp`;
`cd $outdir && samtools sort $bam_tmp $bam`;
`cd $outdir && samtools index $bam\.bam && rm -rf $bam_tmp`;
print "Done~\n";

