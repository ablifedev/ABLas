#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions( \%opts, "i=s", "o=s", "h" );

if ( !defined( $opts{i} ) || !defined( $opts{o} ) ) {
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i           gfffile          must be given

			-o           bedfile          must be given

	Usage End.

	exit;
}

my $infile  = $opts{i};
my $outfile = $opts{o};

my $forward_chr     = "";
my $isoform         = "";
my $forward_isoform = "";
my $forward_start   = 0;
my $forward_end     = 0;
my $left            = 0;
my $right           = 0;
my $flag            = 0;

my $key               = "";
my %negative_junction = ();
my %positive_junction = ();

open OUT, ">$outfile" || die;
open IN,  "$infile"   || die;

while (<IN>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\s+/;
	next if ( $feature =~ /chromosome/ || $feature =~ /gene/ || $feature =~ /protein/ || $_ =~ /^\#/);
	if($feature =~ /^mRNA|transcript$/){
		$info =~ m/ID=([\w\:\.\-]+)/;
		$isoform="$1";
	}
	if($feature =~ /^exon$/){
		$info =~ m/Parent=([\w\:\.\-]+)/;
		$isoform="$1";
	}

	if ( $isoform eq $forward_isoform ) {
		if ( $feature eq "exon" ) {
			if ( $flag > 0 ) {
				if ( $strand eq "+" ) {
					$left  = $forward_end - 50;
					$right = $start + 49;
					$key   = $chr . ":" . $left . ":" . $right;
					if ( !defined( $positive_junction{$key} ) ) {
						$positive_junction{$key} = 1;
						print OUT
"$forward_chr\t$left\t$right\tJUNC\t50\t$strand\t$left\t$right\t255,0,0\t2\t50,50\t0,0\n";
					}
				}
				else {
					$left  = $end - 50;
					$right = $forward_start + 49;
					$key   = $chr . ":" . $left . ":" . $right;
					if ( !defined( $negative_junction{$key} ) ) {
						$negative_junction{$key} = 1;
						print OUT
"$forward_chr\t$left\t$right\tJUNC\t50\t$strand\t$left\t$right\t255,0,0\t2\t50,50\t0,0\n";
					}
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
	$forward_chr     = $chr;
}
close IN;
close OUT;
