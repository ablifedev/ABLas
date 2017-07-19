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
		
			-i	    all_cat_sj_file	must be given

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
my %junction = ();

open IN, $fuhetiaojian_file   || die;
open OUT,  ">$out_file" || die;
while (<IN>) {    #把所有的AS添加进hash

	#chr20   348067  348258  JUNC00000021    3       +       348067  348258  255,0,0 2       45,57   0,134
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	# print $blockSizes,"\n";
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;
	my $key = "$chr:$start:$end:$strand";
	if ( !defined( $junction{$key} ) ) {
  	print OUT $_, "\n";
  	$junction{$key}=1;
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

