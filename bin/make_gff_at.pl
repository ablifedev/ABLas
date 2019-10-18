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
GetOptions( \%opts, "i=s", "o=s", "h" );

if ( !defined( $opts{i} ) || !defined( $opts{o} ) ) {
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i			infile          must be given

			-o			outfile         must be given
			

	Usage End.

	exit;
}

my $infile  = $opts{i};
my $outfile = $opts{o};

my $gene    = "";
my $isoform = "";
my %region  = ();
open OUT, ">$outfile" || die;
open IN,  $infile     || die;
while (<IN>) {
	chomp;
	next if(/^#/);
	my @line = split(/\t/);

	#	if ($line[2] =~ /chromosome/){
	#		print OU "$line[0]\t$line[4]\n";
	#	}
	# if ( $line[8] =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/ ) {
	# 	$gene    = $2;
	# 	$isoform = "$2.$3";
	# }
	
		my $key =
	    $line[0] . ":"
	  . $line[2] . ":"
	  . $line[6] . ":"
	  . $line[3] . ":"
	  . $line[4] ;
	next if defined($region{$key});
	$region{$key}=1;
	
	if ( $line[2] =~ /^mRNA$/ ) {
		print OUT "$line[0]\tmRNA\t$line[6]\t$line[3]\t$line[4]\n";
	}
	if ( $line[2] =~ /^exon$/ ) {
		print OUT "$line[0]\texon\t$line[6]\t$line[3]\t$line[4]\n";
	}
	if ( $line[2] =~ /intron/ ) {
		print OUT "$line[0]\tintron\t$line[6]\t$line[3]\t$line[4]\n";
	}
	if ( $line[2] =~ /five/ ) {
		print OUT "$line[0]\tfive\t$line[6]\t$line[3]\t$line[4]\n";
	}
	if ( $line[2] =~ /three/ ) {
		print OUT "$line[0]\tthree\t$line[6]\t$line[3]\t$line[4]\n";
	}
	if ( $line[2] =~ /intergenic/ ) {
		print OUT "$line[0]\tintergenic\t\\\t$line[3]\t$line[4]\t\\\n";
	}
}
close IN;
close OUT;
