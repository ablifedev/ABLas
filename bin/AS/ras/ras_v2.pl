#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#发现不同转录组中发生显著变化的剪接事件，将两组AS结果按照AS事件junction所占的比例进行比较

my %opts;
GetOptions( \%opts, "i1=s", "i2=s", "j1=s", "j2=s","dif=s","o=s", "h" );

if (   !defined( $opts{i1} )
	|| !defined( $opts{i2} )
	|| !defined( $opts{j1} )
	|| !defined( $opts{j2} )
	|| !defined( $opts{dif} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ras(2 vs 1)
		
				Version: $ver
		
		Usage:perl $0
		
			-i1		asout file 1	must be given
			
			-i2		asout file 2	must be given

			-j1		junction file 1	must be given

			-j2 	junction file 2	must be given

			-dif	dif of ratio	must be given

			-o		outfile     must be given

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $file1    = $opts{i1};
my $file2    = $opts{i2};
my $file3    = $opts{j1};
my $file4    = $opts{j2};
my $dif    = $opts{dif};
my $out_file = $opts{o};

my $dif_fu = 0 - $dif;
my $key       = "";
my %junction  = ();
my %junction_p  = ();
my %junction2 = ();

my %readsnum_1 =();
my %readsnum_2 =();

my $mytemp = "";
my $mychange = 0;

open IN1, $file1       || die;
open IN2, $file2       || die;
open IN3, $file3       || die;
open IN4, $file4       || die;
open OUT, ">$out_file" || die;
my $file = $out_file."_up";
open OUTUP, ">$file" || die;
$file = $out_file."_down";
open OUTDOWN, ">$file" || die;

print OUT "#chr\tstart\tend\tstrand\tAS_TYPE\ti1|percent(reads/majortypereads)\ti2|percent(reads/majortypereads)\ti2-i1\tgenesymbol\n";

my %junc1=();
my %junc2=();

my $start          = 0;
my $end            = 0;

while (<IN3>) {              #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	$start = $left + $leftsize;
	$end   = $right - $rightsize + 1;

	$junc1{$start}=$readsNum;
	$junc1{$end}=$readsNum;
}
close IN3;

while (<IN4>) {              #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	$start = $left + $leftsize;
	$end   = $right - $rightsize + 1;

	$junc2{$start}=$readsNum;
	$junc2{$end}=$readsNum;
}
close IN4;



while (<IN1>) {    #把所有的junction添加进hash

	#chr14   23771834        23772265        +       5pMXE   0.18    2       9
	chomp;
	next if ( $_ =~ /^#/ );
	my @line = split;
	$key =
	    $line[0] . ":"
	  . $line[1] . ":"
	  . $line[2] . ":"
	  . $line[3] . ":"
	  . $line[4];
	$junction{$key} = $line[5]."(".$line[6]."/".$line[7].")";
	$junction_p{$key} = $line[5];
	$readsnum_1{$key} = $line[6];
}
close IN1;

while (<IN2>) {    #把所有的junction添加进hash

	#chr1    933522  939227  +       5pMXE   0.22    2       7       ISG15.1.1
	chomp;
	next if ( $_ =~ /^#/ );
	my @line = split;
	$key =
	    $line[0] . ":"
	  . $line[1] . ":"
	  . $line[2] . ":"
	  . $line[3] . ":"
	  . $line[4];
	$junction2{$key} = $line[5];
	$readsnum_2{$key} = $line[6];
	$mytemp = $line[5]."(".$line[6]."/".$line[7].")";
	if ( defined( $junction{$key} ) ) {	
		next if $readsnum_1{$key}+$readsnum_2{$key}<10;
		$mychange = $line[5]-$junction_p{$key};
		print OUT
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$mytemp\t$mychange\t$line[8]\n";
		print OUTUP
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$mytemp\t$mychange\t$line[8]\n" if $mychange>=$dif;
		print OUTDOWN
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$mytemp\t$mychange\t$line[8]\n" if $mychange<=$dif_fu;
	}
	else {
		next if $readsnum_2{$key}<5;
		if(not defined($junc1{$line[1]}) && not defined($junc1{$line[2]})){
			next;
		}
		$junc1{$line[1]}=0 if not defined($junc1{$line[1]});
		$junc1{$line[2]}=0 if not defined($junc1{$line[2]});
		my $reads_num = $junc1{$line[1]} + $junc1{$line[2]};
		print OUT
		  "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t0\t$mytemp\t$line[5]\t$line[8]\n";
		print OUTUP
		  "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$reads_num\t$mytemp\t$line[5]\t$line[8]\n" if $line[5]>=$dif;
	}
}
close IN2;

open IN1, $file1 || die;
while (<IN1>) {    #把所有的junction添加进hash

	#chr1    933522  939227  +       5pMXE   0.22    2       7       ISG15.1.1
	chomp;
	next if ( $_ =~ /^#/ );
	my @line = split;
	$key =
	    $line[0] . ":"
	  . $line[1] . ":"
	  . $line[2] . ":"
	  . $line[3] . ":"
	  . $line[4];
	if ( !defined( $junction2{$key} ) ) {
		next if $readsnum_1{$key}<5;
		if(not defined($junc2{$line[1]}) && not defined($junc2{$line[2]})){
			next;
		}
		$junc2{$line[1]}=0 if not defined($junc2{$line[1]});
		$junc2{$line[2]}=0 if not defined($junc2{$line[2]});
		my $reads_num = $junc2{$line[1]} + $junc2{$line[2]};
		$mychange = 0-$junction_p{$key};
		print OUT
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t0\t$mychange\t$line[8]\n";
		print OUTDOWN
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$reads_num\t$mychange\t$line[8]\n" if $mychange<=$dif_fu;
	}
}
close IN1;

close OUT;

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