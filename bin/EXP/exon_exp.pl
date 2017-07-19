#!/usr/bin/perl -w
# 
# Copyright (c)   AB_Life 2011
# Writer:         xuxiong <xuxiong19880610@163.com>
# Program Date:   2012.03.16
# Modifier:       xuxiong <xuxiong19880610@163.com>
# Last Modified:  2011.03.16
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $exp = $ARGV[0];
my $count=0;

open IN ,"$exp";
while(<IN>){
	chomp;
	next if($_!~/^\d+/);
	next if(/^BaseNumber/);
	my @line=split;
	my $total = sum(@line);
	$count++ if $total>0;
}

$exp=~s/_EXP_exon_base$//;
print $exp,"\t",$count,"\n";