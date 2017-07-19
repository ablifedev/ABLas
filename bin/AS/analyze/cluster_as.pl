#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#将初步cluster的AS事件进行分型。

my %opts;
GetOptions( \%opts, "i=s", "o=s", "h" );

if (   !defined( $opts{i} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i		infile	must be given
			
			-o		outfile	must be given

	Usage End.

	exit;
}

#############Time_start#############
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $infile   = $opts{i};
my $outfile   = $opts{o};



open GFFAS, $gffas_file || die;
while (<GFFAS>) {    #把所有的GFFjunction添加进hash
	#chr1	1240861	1246239	-	ES	0.33	2	2	CPSF3L.1.5
	chomp;
	my @line = split;
	$key = $line[0] . ":" . $line[1] . ":" . $line[2]. ":" . $line[4];
	if ( $line[3] eq "+" ) {
		$positive_junction{$key} = 1;
	}
	if ( $line[3] eq "-" ) {
		$negative_junction{$key} = 1;
	}
}
close GFFAS;

open ASOUT, $asout_file      || die;
open OUTK,  ">$knownas_file" || die;
open OUTN,  ">$newas_file"   || die;
while (<ASOUT>) {

	#chr1	1225148	1225401	-	ES	0.25	2	6	ACAP3.1.2
	chomp;
	# next if ( $_ !~ /^c/i );
	my @line = split;
	$key = $line[0] . ":" . $line[1] . ":" . $line[2]. ":" . $line[4];
	if ( $line[3] eq "+" ) {
		if ( defined( $positive_junction{$key} ) ) {
			print OUTK $_ , "\n";
		}
		else {
			print OUTN $_, "\n";
		}
	}
	else {
		if ( defined( $negative_junction{$key} ) ) {
			print OUTK $_, "\n";
		}
		else {
			print OUTN $_, "\n";
		}
	}
}
close ASOUT;
close OUTK;
close OUTN;



############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";
#######Sub_format_datetime#######
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}

