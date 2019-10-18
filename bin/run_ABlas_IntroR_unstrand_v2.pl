#!/usr/bin/perl -w
my $ver = "1.0.0";

##v2
#1. 将读取sam转换为读取bam，处理基因同时读取bam增加速度减少内存占用

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use Bio::DB::Sam;

my %opts;
GetOptions( \%opts, "bam=s", "fa=s","chrlen=s", "gff=s", "o=s", "anno=s", "od=s","h" );

if (   !defined( $opts{bam} )
	|| !defined( $opts{chrlen} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| !defined( $opts{fa} )
	|| !defined( $opts{anno} )
	|| !defined( $opts{od} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-bam         bam file        must be given;
		
		-chrlen      chromsome lenght       must be given;

		-gff         gff file               must be given;

		-fa         fa file               must be given;

		-o           out file               must be given;

		-anno        anno file              must be given;

		-od          out directory          must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $bam_file    = $opts{bam};
my $fa_file    = $opts{fa};
my $chrlen_file = $opts{chrlen};
my $gff_file    = $opts{gff};
my $out_file    = $opts{o};
my $anno_file    = $opts{anno};
my $outdir    = $opts{od};

my $forward_chr="";

open OUT, ">run_intronR.sh";

open CHRLEN, $chrlen_file || die;
while (<CHRLEN>) {
	chomp;
	next if(/\#/);
	my @line=split(/\t/);
	if($line[0] ne $forward_chr && $forward_chr ne ""){
		print OUT "perl5.14 $Bin/ABlas_IntroR_unstrand_v2.pl -bam $bam_file -fa $fa_file -chrlen $chrlen_file -gff $gff_file -o $out_file"."_$forward_chr"." -anno $anno_file -od $outdir -chr $forward_chr && \n";
	}
	$forward_chr = $line[0];
}
close CHRLEN;
print OUT "perl5.14 $Bin/ABlas_IntroR_unstrand_v2.pl -bam $bam_file -fa $fa_file -chrlen $chrlen_file -gff $gff_file -o $out_file"."_$forward_chr"." -anno $anno_file -od $outdir -chr $forward_chr &&\n";
`perl /public/bin/qsub-sge.pl --queue all.q --resource vf=10.0G --maxproc 60 run_intronR.sh`;
`cd $outdir && cat $out_file\* > $out_file `;
`rm -rf $out_file\_\*`;



############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

my $time_used = time() - $start_time;
my $h         = $time_used / 3600;
my $m         = $time_used % 3600 / 60;
my $s         = $time_used % 3600 % 60;
printf( "\nAll Time used : %d hours\, %d minutes\, %d seconds\n\n", $h, $m,
	$s );

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
