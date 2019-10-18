#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use Bio::DB::Sam;

my $buf = "-";
my %opts;
GetOptions( \%opts, "gff=s", "fa=s", "bam=s","readlen=s","oh=s","o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{fa} )
	|| !defined( $opts{bam} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff        annotation file        must be given;

		-fa         fa                     must be given;

		-bam        bam                    must be given;

		-readlen      readlen              default is 75

		-oh        overhang                   default is 8

		-o          out file               must be given;

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
my $fa_file = $opts{fa};
my $bam_file = $opts{bam};
my $out_file = $opts{o};
my $readlen = 75;
$readlen = $opts{readlen} if defined($opts{readlen});
my $oh = 8;
$oh = $opts{oh} if defined($opts{oh});
open OUT, ">$out_file" || die;

my $sam = Bio::DB::Sam->new(-fasta=>"$fa_file",
                            -bam  =>"$bam_file");

open GFF, "$gff_file" || die;
while(<GFF>){
	chomp;
	# Chr1	AT1G01490.1	IR	182358	182061	.	-	3	flanks "182358-,182061^"; structure "1^2-,0"; splice_chain "182253^182135-,"; splice_position "182253,182135"; ID "IR_11";
	# Chr1	AT1G01540.1	IR	197775	198183	.	+	147	flanks "197775-,198183^"; structure "1^2-,0"; splice_chain "197904^197974-,"; splice_position "197904,197974"; ID "IR_12";
	my @line = split(/\t/);
	$line[-1]=~/splice_position\s\"(\d+)\,(\d+)\"\;/;
	my $left = 0;
	my $right = 0;
	if($line[6] eq "+"){
		$left = $1;
		$right = $2;
	}else{
		$left = $2;
		$right = $1;
	}
	my $start = $left - $readlen + 4;
	my $end = $right - 4;
	my $reads = 0;
	my @alignments = $sam->get_features_by_location(-seq_id => $line[0],
                                                        -start  => $start,
                                                        -end    => $end);
	for my $a (@alignments) {
		my $cigar     = $a->cigar_str;
		my $thisstart  = $a->start;
		next if $cigar =~ /^\d+M\d+N\d+M$/;
		$reads++ if($thisstart>=$start && $thisstart<=$end);
	}

	my $ir_density = $reads/($end - $start + 1);
	my $sj_density = $line[7]/($readlen+1-2*$oh);
	my $tmp = $line[7];
	$line[5]=$ir_density;
	$line[7]=$sj_density;
	print OUT join("\t",@line),"\t",$tmp,"\n";
}

close GFF;
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