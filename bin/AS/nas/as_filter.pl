#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $buf = "-";
my %opts;
GetOptions( \%opts, "as=s", "num_alt_avg=s", "num_model_avg=s", "num_alt_min=s", "num_model_min=s", "num_alt_ratio=s", "num_model_ratio=s","num_sum_alt_model=s","o=s", "h" );

if (   !defined( $opts{as} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-as                   as result file        must be given;

		-num_alt_avg          Require at least N-many alt average sj.                          optional,default is 0;

		-num_model_avg        Require at least N-many model average sj.                        optional,default is 0;

		-num_alt_min          Require at least N-many alt min sj.                              optional,default is 0;

		-num_model_min        Require at least N-many model min sj.                            optional,default is 0;

		-num_sum_alt_model    Require the sum of alt and model average sjs to be at least N.   optional,default is 0;

		-num_alt_ratio        Require at least N alt ratio                                     optional,default is 0;

		-num_model_ratio      Require at least N model ratio                                   optional,default is 0;

		-o                    out file                                                          must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $as_file = $opts{as};
my $num_alt_avg = $opts{num_alt_avg} || 0;
my $num_model_avg = $opts{num_model_avg} || 0;
my $num_alt_min = $opts{num_alt_min} || 0;
my $num_model_min = $opts{num_model_min} || 0;
my $num_sum_alt_model = $opts{num_sum_alt_model} || 0;
my $num_alt_ratio = $opts{num_alt_ratio} || 0;
my $num_model_ratio = $opts{num_model_ratio} || 0;
my $out_file = $opts{o};
open OUT, ">$out_file" || die;

&load_as($as_file);

close OUT;


# Chr01   724121  724320  +       724121  724210  A3SS    0.00    0       16      0       16      PAC:26822888
#子过程-----load_as($as_file)
sub load_as {
	my ($as_file) = @_;
	my $count = 0;
	open AS, $as_file    || die;
	while (<AS>) {
		chomp;
		print OUT $_,"\n" if ( $_ =~ /^\#/ );
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $alt_avg=$line[8];
		my $model_avg=$line[9];
		my $alt_min=&minsj($line[10]);
		my $model_min=&minsj($line[11]);
		my $sum_alt_model=&sumsj($line[10]) + &sumsj($line[11]);
		my $alt_ratio = $line[7];
		my $model_ratio = 1 - $line[7];
		next if ($alt_avg<$num_alt_avg || $model_avg<$num_model_avg || $alt_min<$num_alt_min || $model_min<$num_model_min || $sum_alt_model<$num_sum_alt_model || $alt_ratio<$num_alt_ratio || $model_ratio<$num_model_ratio );
		print OUT join("\t",@line),"\n";
	}
	close AS;
	#调试信息
	print "done reading as ...........", "\n\n";
}

sub sumsj {
	my ($sjlist) = @_;
	my @num = split(/;/,$sjlist);
	return sum(@num);

}

sub minsj {
	my ($sjlist) = @_;
	my @num = split(/;/,$sjlist);
	return min(@num);

}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/;/,$sjlist);
	my $avg = (sum(@num)) / (scalar(@num)) ;
	return $avg;
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