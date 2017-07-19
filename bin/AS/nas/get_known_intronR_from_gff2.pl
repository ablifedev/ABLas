#!/usr/bin/perl -w

my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my %opts;
GetOptions( \%opts, "gff=s", "o=s","gm=s" );

if (   !defined( $opts{o} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{gm} ) )
{
	print <<"		Usage End.";

			Version:$ver

		Usage:perl $0

			-gff          GFF3 file        must be given 

			-gm           gene model       must be given
			
			-o            outfile          must be given			

		Usage End.

	exit;
}

my $gff_file     = $opts{gff};
my $out_file     = $opts{o};
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
	next if ( $feature =~ /chromosome/ || $_ =~ /^\#/ || $feature =~ /gene/ );
	if ( $feature eq "mRNA" ) {
		$info =~ m/ID=([\w\.]+);Parent=(\w+);/;
		my $gene = $2;
		$gene{$gene} = $start;
	}
}
close GFF;

open GFF, $gff_file || die;
open OUT, ">$out_file"  || die;
while (<GFF>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\s+/;
	next if ( $_=~/^#/ || $feature =~ /chromosome/ || $feature =~ /gene/ || $feature =~ /protein/ );
	if($info =~ m/Parent=(\w+)\.(\d+)/){
		$gene="$1";
		$isoform = "$1.$2";
	}elsif($info =~ m/ID=([\w\.]+);Parent=(\w+);/){
		$gene="$2";
		$isoform = "$1";
	}
	$flag =0 if $feature =~ /mRNA/;
	if (   ( defined( $gene{$gene} ) )
		&& ( defined( $gene{$forward_gene} ) )
		&& ( $gene{$gene} == $gene{$forward_gene} ) )
	{
		if ( defined( $gene_model{$gene} ) ) {
			if ( $gene_model{$gene} eq $isoform ) {
				if ( $isoform eq $forward_isoform && $feature eq "exon" ) {
					push @exon, "$start:$end";
				}
			}
			else {
				if ( $isoform eq $forward_isoform && $feature eq "exon" ) {
					if ( $flag > 0 ) {
						if ( $strand eq "+" ) {
							print "$strand\t$forward_end\t$start\n";
							push @intron, "$forward_end:$start";
						}
						else {
							print "$strand\t$end\t$forward_start\n";
							push @intron, "$end:$forward_start";
						}
					}
					$flag++;
					$forward_start = $start;
					$forward_end   = $end;
				}
			}
		}
		$forward_isoform = $isoform;
	}
	else {
		$flag = 0;
		
		foreach my $j (@exon) {
			my ( $jstart, $jend ) = split /\:/, $j;
			foreach my $i (@intron) {
				my ( $istart, $iend ) = split /\:/, $i;
				if ( $jstart <= $istart && $jend >= $iend ) {
					$intronR{"$forward_chr:$istart:$iend"} = 1;
					print OUT "$forward_chr:$istart:$iend:$jstart:$jend\n";
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


close OUT;
close GFF;
