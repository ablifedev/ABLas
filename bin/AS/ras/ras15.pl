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
GetOptions( \%opts, "i1=s", "i2=s", "o=s", "h" );

if (   !defined( $opts{i1} )
	|| !defined( $opts{i2} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i1		asout file	must be given
			
			-i2		asout file	must be given

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
my $out_file = $opts{o};

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
open OUT, ">$out_file" || die;
my $file = $out_file."_up";
open OUTUP, ">$file" || die;
$file = $out_file."_down";
open OUTDOWN, ">$file" || die;

print OUT "#chr\tstart\tend\tstrand\tAS_TYPE\ti1|percent(reads/majortypereads)\ti2|percent(reads/majortypereads)\ti2-i1\tgenesymbol\n";

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
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$mytemp\t$mychange\t$line[8]\n" if $mychange>=0.15;
		print OUTDOWN
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t$mytemp\t$mychange\t$line[8]\n" if $mychange<=-0.15;
	}
	else {
		next if $readsnum_2{$key}<5;
		print OUT
		  "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t0\t$mytemp\t$line[5]\t$line[8]\n";
		print OUTUP
		  "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t0\t$mytemp\t$line[5]\t$line[8]\n" if $line[5]>=0.15;
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
		$mychange = 0-$junction_p{$key};
		print OUT
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t0\t$mychange\t$line[8]\n";
		print OUTDOWN
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$junction{$key}\t0\t$mychange\t$line[8]\n" if $mychange<=-0.15;
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