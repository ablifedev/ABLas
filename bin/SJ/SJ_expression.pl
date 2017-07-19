#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions( \%opts, "bed=s", "gff=s", "mod=s", "o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{mod} )
	|| !defined( $opts{bed} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-bed	bed file	must be given

			-gff	gff file	must be given

			-mod	gene model	must be given

			-o		outfile     must be given

	Usage End.

	exit;
}

my $bed_file = $opts{bed};
my $gff_file = $opts{gff};
my $mod_file = $opts{mod};
my $out_file = $opts{o};

my %gene_model = ();
open MOD, $mod_file || die;
while (<MOD>) {
	chomp;
	m/(\S+)\.(\d+)$/;
	$gene_model{$1} = "$1.$2";
}
close MOD;

my %junction        = ();
my %num             = ();
my $gene            = "";
my $forward_gene    = "";
my $isoform         = "";
my $forward_isoform = "";
my $forward_start   = 0;
my $forward_end     = 0;
my $left            = 0;
my $right           = 0;
my $flag            = 0;
my $i               = 0;

my $count = 0;
open GFF, $gff_file || die;

while (<GFF>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\s+/;
	next
	  if ( $feature =~ /chromosome/
		|| $_ =~ /^\#/
		|| $_ =~ /gene/
		|| $_ =~ /ATC/
		|| $_ =~ /ATM/ );
	$info =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/;
	$gene    = $2;
	$isoform = "$2.$3";
	if (   ( defined( $gene_model{$gene} ) )
		&& ( $gene_model{$gene} eq $isoform ) )
	{

		if ( $isoform eq $forward_isoform ) {
			if ( $feature eq "exon" ) {
				if ( $flag > 0 ) {
					if ( $strand eq "+" ) {
						$left                          = $forward_end - 50;
						$right                         = $start + 49;
						$junction{"$chr:$left:$right"} = $forward_gene;
					}
					else {
						$left                          = $end - 50;
						$right                         = $forward_start + 49;
						$junction{"$chr:$left:$right"} = $forward_gene;
					}
				}
				$flag++;
				$forward_start = $start;
				$forward_end   = $end;
			}
		}
		else {
			$flag = 0;
		}
		$forward_isoform = $isoform;
		$forward_gene    = $gene;
	}
}
close GFF;

print $count;

open BED, $bed_file    || die;
open OUT, ">$out_file" || die;
while (<BED>) {
	chomp;
	my @line = split;
	$line[3] =~ /\((\d+)\)/;
	my $read = $1;
	$num{"$line[0]:$line[1]:$line[2]"} = $read;
}
close BED;

foreach my $keys ( sort keys %junction ) {
	my @a = split /\:/, $keys;
	if ( defined $num{$keys} ) {
		print OUT "$junction{$keys}:$a[1]:$a[2]", "\t$num{$keys}\n";
	}
	else {
		print OUT "$junction{$keys}:$a[1]:$a[2]", "\t0\n";
	}
}
close OUT;
