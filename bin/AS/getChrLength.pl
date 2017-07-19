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
		
			-i			samFile               must be given

			-o			chrLengthFile         must be given
			

	Usage End.

	exit;
}

my $infile  = $opts{i};
my $outfile = $opts{o};

my @line    = ();
my $chr = "";
my $chrLength = "";
open OUT, ">$outfile" || die;
open IN,  $infile     || die;
while (<IN>) {
	chomp;
	@line = split;

	if ( $line[0] =~ /^\@SQ$/ ) {
		$line[1] =~ /^SN:(\w+)$/;
		$chr = $1;
		$line[2] =~ /^LN:(\w+)$/;
		$chrLength = $1;
		print OUT "$chr\t$chrLength\n";
	}
	
	if ( $line[0] !~ /^\@/ ) {
	last;
	}

}
close IN;
close OUT;
