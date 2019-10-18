#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;

my $buf = "-";
my %opts;
GetOptions( \%opts, "gffsj=s","sj=s", "ok=s", "on=s","h" );

if (   !defined( $opts{gffsj} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{ok} )
	|| !defined( $opts{on} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-sj           junctions.bed file        must be given;

		-gffsj        annotation sj file        must be given;

		-ok           konwn sj outfile          must be given;

		-on           new sj outfile            must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################
my $sj_file = $opts{sj};
my $gffsj_file = $opts{gffsj};
my $ok_file = $opts{ok};
my $on_file = $opts{on};
open OK, ">$ok_file" || die;
open ON, ">$on_file" || die;

my %sj=();
&load_gffsj($gffsj_file);


&load_sj($sj_file);


close OK;
close ON;


###子过程-----load_sj($sj_file)
sub load_sj {
	my ($sj_file) = @_;
	open SJ, $sj_file || die;
	while (<SJ>) {              #把所有的junction添加进hash
		#默认为读取通过tophat获取的junctions.bed文件
		#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
		#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
		chomp;
		next if ( $_ =~ /^track/i );
		my $start          = 0;
		my $end            = 0;
		my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
			$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
		my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
		if ( $strand eq "+" ) {
			$start = $left + $leftsize;
			$end   = $right - $rightsize + 1;
			if(defined($sj{$chr}->{$strand}->{"$start:$end"})){
				print OK $_,"\n";
			}else{
				print ON $_,"\n";
			}
		}
		else {
			$start = $right - $rightsize + 1;
			$end   = $left + $leftsize;
			if(defined($sj{$chr}->{$strand}->{"$start:$end"})){
				print OK $_,"\n";
			}else{
				print ON $_,"\n";
			}
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

###子过程-----load_gffsj($sj_file)
sub load_gffsj {
	my ($sj_file) = @_;
	open SJ, $sj_file || die;
	while (<SJ>) {              #把所有的junction添加进hash
		#gff的sj
		#chr20   199843  204678  -
		chomp;
		next if ( $_ =~ /^track/i );
		my ( $chr, $start, $end, $strand) = split /\t/ ;
		$sj{$chr}->{$strand}->{"$start:$end"} = 1;
	}
	close SJ;
	#调试信息
	print "done reading GFFSJ ...........", "\n\n";
}

############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

my $time_used = time() - $start_time;
my $h = $time_used/3600;
my $m = $time_used%3600/60;
my $s = $time_used%3600%60;
printf("\nAll Time used : %d hours\, %d minutes\, %d seconds\n\n",$h,$m,$s);


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