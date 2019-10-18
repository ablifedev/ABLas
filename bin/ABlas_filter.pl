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
GetOptions( \%opts, "gff=s", "num_alt_avg=s", "num_model_avg=s", "num_alt_min=s", "num_model_min=s", "num_alt_ratio=s", "num_model_ratio=s","num_sum_alt_model=s","o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff              as result file        must be given;

		-num_alt_avg          Require at least N-many alt average sj.           must be given;

		-num_model_avg          Require at least N-many model average sj.           must be given;

		-num_alt_min          Require at least N-many alt min sj.           must be given;

		-num_model_min          Require at least N-many model min sj.           must be given;

		-num_sum_alt_model  Require the sum of alt and model average sjs to be at least N.   must be given;

		-num_alt_ratio        Require at least N alt ratio                      must be given;

		-num_model_ratio        Require at least N model ratio                      must be given;

		-o                out file               must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $gff_file = $opts{gff};
my $num_alt_avg = $opts{num_alt_avg} || 0;
my $num_model_avg = $opts{num_model_avg} || 0;
my $num_alt_min = $opts{num_alt_min} || 0;
my $num_model_min = $opts{num_model_min} || 0;
my $num_sum_alt_model = $opts{num_sum_alt_model} || 0;
my $num_alt_ratio = $opts{num_alt_ratio} || 0;
my $num_model_ratio = $opts{num_model_ratio} || 0;
my $out_file = $opts{o};
open OUT, ">$out_file" || die;

&load_gff($gff_file);

close OUT;


# chr1	uc009viu.1	SE	8776	8232	.	-	1	flanks "8776^,8232-"; structure "0,1-2^"; splice_chain ",8417-8364^"; skipped_exon_no "1"; flank_exon_start "9622,8232"; flank_exon_end "8776,8131"; skipped_exon_start "8417"; skipped_exon_end "8364"; splice_position "8417,8364";
#子过程-----load_gff($gff_file)
sub load_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		print OUT $_,"\n" if ( $_ =~ /^\#/ );
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $alt_avg=&avgsj($line[5]);
		my $model_avg=&avgsj($line[7]);
		my $alt_min=&minsj($line[5]);
		my $model_min=&minsj($line[7]);
		my $sum_alt_model=&sumsj($line[5]) + &sumsj($line[7]);
		my $alt_ratio = $line[9];
		my $model_ratio = 1 - $line[9];
		next if ($alt_avg<$num_alt_avg || $model_avg<$num_model_avg || $alt_min<$num_alt_min || $model_min<$num_model_min || $sum_alt_model<$num_sum_alt_model || $alt_ratio<$num_alt_ratio || $model_ratio<$num_model_ratio );
		print OUT join("\t",@line),"\n";
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

sub sumsj {
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	return sum(@num);

}

sub minsj {
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	return min(@num);

}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
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