#!/usr/bin/perl -w

my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#分别获取所有已知的和新的intronR事件。

my %opts;
GetOptions( \%opts, "gff=s", "i=s", "o=s", "ok=s", "gm=s" );

if (   !defined( $opts{o} )
	|| !defined( $opts{ok} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{i} )
	|| !defined( $opts{gm} ) )
{
	print <<"		Usage End.";

			Version:$ver

		Usage:perl $0

			-gff          GFF3 file        must be given 

			-gm           gene model       must be given

			-i            intronR file     must be given
			
			-o            novel_intronR          must be given
			
			-ok           known_intronR           must be given

		Usage End.

	exit;
}

my $gff_file     = $opts{gff};
my $intronR_file = $opts{i};
my $out_file     = $opts{o};
my $ok_file      = $opts{ok};
my $mod_file     = $opts{gm};

my %gene_model = ();
open MOD, $mod_file || die;
while (<MOD>) {
	chomp;
	next if ( $_ =~ /chrom/ );
	my @line = split;
	m/(\S+)\.(\d+)$/;
	$gene_model{$1} = "$1.$2";
}
close MOD;

my %gene = ();
my $forward_start;
my $forward_end;
my $forward_gene    = "";
my $forward_isoform = "";
my $forward_chr     = "";
my $flag            = 0;
my $count           = 0;
my $gene;
my $isoform;
my @intron  = ();
my @exon    = ();
my %intronR = ();
open GFF, $gff_file || die;

while (<GFF>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\s+/;
	next if ( $feature =~ /chromosome/ || $_ =~ /^\#/ || $_ =~ /gene/ );
	if ( $feature eq "mRNA" ) {
		$info =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/;
		my $gene = $2;
		$gene{$gene} = $start;
	}
}
close GFF;

open GFF, $gff_file || die;
while (<GFF>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\t/;
	next if ( $feature =~ /chromosome/ || $_ =~ /^\#/ || $_ =~ /gene/ );
	$info =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/;
	my $gene    = $2;
	my $isoform = "$2.$3";
	if (   ( defined( $gene{$gene} ) )
		&& ( defined( $gene{$forward_gene} ) )
		&& ( $gene{$gene} == $gene{$forward_gene} ) )
	{
		if ( defined( $gene_model{$gene} ) ) {
			if ( $gene_model{$gene} eq $isoform ) {
				if ( $isoform eq $forward_isoform && $feature eq "exon" ) {
					if ( $flag > 0 ) {
						if ( $strand eq "+" ) {
							push @intron, "$forward_end:$start";
						}
						else {
							push @intron, "$end:$forward_start";
						}
					}
					$flag++;
					$forward_start = $start;
					$forward_end   = $end;
				}
			}
			else {
				if ( $isoform eq $forward_isoform && $feature eq "exon" ) {
					push @exon, "$start:$end";
				}
			}
		}
		$forward_isoform = $isoform;
	}
	else {
		$flag = 0;
		foreach my $i (@intron) {
			my ( $istart, $iend ) = split /\:/, $i;
			foreach my $j (@exon) {
				my ( $jstart, $jend ) = split /\:/, $j;
				if ( $jstart <= $istart && $jend >= $iend ) {
					$intronR{"$forward_chr:$istart:$iend"} = 1;
					$count++;
					last;
				}
			}
		}
		undef(@exon);
		undef(@intron);
	}
	$forward_gene = $gene;
	$forward_chr  = $chr;
}
print "$count\n";
open OUT, ">$out_file"  || die;
open OK,  ">$ok_file"   || die;
open IN,  $intronR_file || die;
while (<IN>) {
	chomp;
	my ( $chr, $start, $end, $strand, $type, $other, $jr, $er, $gene ) = split;
	next if ( /Chrosome/ || $type =~ /\\/ );
	if ( !defined( $intronR{"$chr:$start:$end"} ) ) {
		print OUT "$chr\t$start\t$end\t$strand\t$type\t$other\t$jr\t$er\t$gene\n";
	}
	else {
		print OK "$chr\t$start\t$end\t$strand\t$type\t$other\t$jr\t$er\t$gene\n";
	}
}
close OUT;
close OK;
close IN;
close GFF;
