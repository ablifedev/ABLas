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
GetOptions( \%opts, "i=s", "gff=s", "ok=s", "on=s", "h" );

if (   !defined( $opts{i} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{ok} )
	|| !defined( $opts{on} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i		asout file	must be given

			-gff	gffas file	must be given

			-ok		knownas file	must be given
			
			-on		newas file	must be given

	Usage End.

	exit;
}

#############Time_start#############
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $asout_file   = $opts{i};
my $gff_file     = $opts{gff};
my $knownas_file = $opts{ok};
my $newas_file   = $opts{on};

my $key               = "";
my %negative_junction = ();
my %positive_junction = ();

open GFF, $gff_file || die;
while (<GFF>) {    #把所有的GFFjunction添加进hash
	 #chr1	hg18_knownGene	intron	114202162	114242012	0.000000	+	.	ID=intron:AK123199.1.1:4;Parent=AK123199.1.1
	chomp;
	my @line = split;
	if ( $line[2] =~ /intron/ ) {
		$key = $line[0] . ":" . ( $line[3] - 1 ) . ":" . ( $line[4] + 1 );
		if ( $line[6] eq "+" ) {
			$positive_junction{$key} = 1;
		}
		if ( $line[6] eq "-" ) {
			$negative_junction{$key} = 1;
		}
	}
}
close GFF;

open ASOUT, $asout_file      || die;
open OUTK,  ">$knownas_file" || die;
open OUTN,  ">$newas_file"   || die;
while (<ASOUT>) {

	#chr1	1017346	1041303	-	IntronR	0.75	3	1	C1orf159.1
	chomp;
	next if ( $_ !~ /^c/i );
	my @line = split;
	if ( $line[4] =~ /IntronR/ ) {
		$key = $line[0] . ":" . $line[1] . ":" . $line[2];
		if ( $line[3] eq "+" ) {
			if ( defined( $positive_junction{$key} ) ) {
				print OUTK $_, "\n";
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

