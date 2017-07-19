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
GetOptions( \%opts, "i=s", "o=s", "h" );

if (   !defined( $opts{i} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-i		exonExpressionFile	must be given

			-o		outfile     		must be given

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $file    = $opts{i};
my $out_file = $opts{o};

my %exon =();
my $key       = "";

my $num = 0;

open IN, $file       || die;
open OUT, ">$out_file" || die;

print OUT "#chr\ttype\tstrand\tstart\tend\tgene\tbaseNum\n";

while (<IN>) {    #提取uniq_exon

	#chr1    1       +       869151  869824  56      1848
	chomp;
	next if ( $_ =~ /^#/i );
	my @line = split;
	$key =
	    $line[0] . ":"
	  . $line[1] . ":"
	  . $line[2] . ":"
	  . $line[3] . ":"
	  . $line[4];
	if (!defined($exon{$key})){
		$exon{$key} = 1;
		print OUT "$_\n";
		$num++;
	}
	
}
close IN;

close OUT;

print "Total exon number: ",$num,"\n\n";
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