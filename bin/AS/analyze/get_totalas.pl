#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#将几个样本发现的available_as_events的条目cat到一起 ,包含intronR,通过本程序挑出
#其中的并集作为最后的TOTAL AS,输出到TOTAL。

my %opts;
GetOptions( \%opts, "i=s",  "o=s",
	"h" );

if (   !defined( $opts{i} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i	    all_cat_as_file	must be given

			-o		outfile     must be given

	Usage End.

	exit;
}

#############Time_start#############
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $fuhetiaojian_file = $opts{i};
my $out_file  = $opts{o};

my $key               = "";
my %negative_junction = ();
my %positive_junction = ();

open IN, $fuhetiaojian_file   || die;
open OUT,  ">$out_file" || die;
while (<IN>) {    #把所有的AS添加进hash

	#chr1	23617089	23623645	-	ES	0.18	9	39	TCEA3.1.1
	chomp;
	# next if ( $_ !~ /^c/i );
	my @line = split;
	$key = $line[0] . ":" . $line[1] . ":" . $line[2]. ":" . $line[4];
	if ( $line[3] eq "+" ) {
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT $_, "\n";
		}
	}
	else {
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT $_, "\n";
		}
	}
}
close IN;

close OUT;

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

