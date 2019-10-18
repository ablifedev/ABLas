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
GetOptions( \%opts, "gffdepth=s", "num_avg_depth=s", "num_sj=s", "num_alt_ratio=s", "num_model_ratio=s","num_mid_depth=s", "o=s", "h" );

if (   !defined( $opts{gffdepth} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gffdepth              as result file        must be given;

		-num_sj          Require at least N-many  sj.           must be given;

		-num_avg_depth          Require at least N-many intron average depth.           must be given;

		-num_mid_depth          Require at least N-many alt min sj.           must be given;

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

my $gff_file = $opts{gffdepth};
my $num_sj = $opts{num_sj} || 0;
my $num_avg_depth = $opts{num_avg_depth} || 0;
my $num_mid_depth = $opts{num_mid_depth} || 0;
my $num_alt_ratio = $opts{num_alt_ratio} || 0;
my $num_model_ratio = $opts{num_model_ratio} || 0;
my $out_file = $opts{o};
open OUT, ">$out_file" || die;

&load_gff($gff_file);

close OUT;


# NW_001867363.1	rna14936,rna14939	IR	1409754	1411273	5335,3.97,3	+	1	flanks "1409754-,1411273^"; structure "1^2-,0"; splice_chain "1409806^1411152-,"; splice_position "1409806,1411152"; ID "IR_1"; ID "IR_1";	0.798792756539235
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
		my ($total_base,$avg_base,$mid_base)=split(/,/,$line[5]);
		my $alt_ratio = $line[9];
		my $model_ratio = 1 - $line[9];
		next if ($line[7]<$num_sj || $avg_base<$num_avg_depth || $mid_base<$num_mid_depth || $alt_ratio<$num_alt_ratio || $model_ratio<$num_model_ratio );
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