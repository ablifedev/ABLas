#!/usr/bin/perl -w

#程序功能：将gff的ucscid转为gene symbol。

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions( \%opts, "g=s", "o=s", "h" );

if ( !defined( $opts{g} ) || !defined( $opts{o} ) ) {
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-g           gfffile         	   must be given

			-o           outfile               must be given

	Usage End.

	exit;
}

my $gfffile = $opts{g};
my $outfile = $opts{o};

my $chr     = "";
my $source  = "";
my $feature = "";
my $start   = "";
my $end     = "";
my $other1  = "";
my $strand  = "";
my $other2  = "";
my $info    = "";

my $forwardmrna_start = 0;
my $forwardmrna_end   = 0;

my $forwardmrna_strand = "";

my $forward_exon_start = 0;
my $forward_exon_end   = 0;

my $forward_exon_info = "";

my $flag = 1;

my $intron = "intron";

my $num = 0;

my $forward_feature = "";
my $forward_line	= "";

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
###################################

open GFF, "$gfffile"  || die;
open OUT, ">$outfile" || die;
while (<GFF>) {
	chomp;
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	) = split /\t/;
	if ( $strand eq "+" ) {
		if ( $feature =~ /mRNA|transcript/ ) {
			if (( $forward_exon_end < $forwardmrna_end )&&( $forwardmrna_strand eq "+" ) ){
				my $intron_start = $forward_exon_end + 1;
				my $intron_end   = $forwardmrna_end;
				$forward_exon_info =~ s/ID=exon/ID=intron/;
				print OUT "$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$forward_exon_info\n";
				$num++;
			}
			if (( $forward_exon_start > $forwardmrna_start )&&( $forwardmrna_strand eq "-" ) ) {
				my $intron_start = $forwardmrna_start;
				my $intron_end   = $forward_exon_start-1;
				$forward_exon_info =~ s/ID=exon/ID=intron/;
				print OUT "$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$forward_exon_info\n";
				$num++;
			}
			$forwardmrna_start = $start;
			$forwardmrna_end   = $end;
			$forwardmrna_strand = $strand;
			$flag              = 1;
			if ( $forward_feature eq "gene" ) {
				print OUT $forward_line, "\n";
			}
		}

		#正向
		if ( $feature =~ /exon/ ) {
			if ( $flag == 1 ) {
				if ( $start > $forwardmrna_start ) {
					my $intron_start = $forwardmrna_start;
					my $intron_end   = $start - 1;
					$info =~ s/ID=exon/ID=intron/;
					print OUT
"$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$info\n";
				}
				$flag = 0;
			}
			else {
				if ( $start > $forward_exon_end + 1 ) {
					my $intron_start = $forward_exon_end + 1;
					my $intron_end   = $start - 1;
					$info =~ s/ID=exon/ID=intron/;
					print OUT
"$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$info\n";
				}
			}
			$forward_exon_start = $start;
			$forward_exon_end   = $end;
			$forward_exon_info  = $info;

		}
		$forward_feature = $feature;
		$forward_line    = $_;
		if ( $feature ne "gene" ) {
			print OUT $_, "\n";
		}
	}
		
		
		
		if ( $strand eq "-" ) {
		if ( $feature =~ /mRNA|transcript/ ) {
			if (( $forward_exon_end < $forwardmrna_end )&&( $forwardmrna_strand eq "+" ) ){
				my $intron_start = $forward_exon_end + 1;
				my $intron_end   = $forwardmrna_end;
				$forward_exon_info =~ s/ID=exon/ID=intron/;
				print OUT "$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$forward_exon_info\n";
				$num++;
			}
			if (( $forward_exon_start > $forwardmrna_start )&&( $forwardmrna_strand eq "-" ) ) {
				my $intron_start = $forwardmrna_start;
				my $intron_end   = $forward_exon_start-1;
				$forward_exon_info =~ s/ID=exon/ID=intron/;
				print OUT "$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$forward_exon_info\n";
				$num++;
			}
			$forwardmrna_start = $start;
			$forwardmrna_end   = $end;
			$forwardmrna_strand = $strand;
			$flag              = 1;
			if ( $forward_feature eq "gene" ) {
				print OUT $forward_line, "\n";
			}
		}

		#负向
		if ( $feature =~ /exon/ ) {
			if ( $flag == 1 ) {
				if ( $end < $forwardmrna_end ) {
					my $intron_start = $end+1;
					my $intron_end   = $forwardmrna_end;
					$info =~ s/ID=exon/ID=intron/;
					print OUT
"$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$info\n";
				}
				$flag = 0;
			}
			else {
				if ( $end + 1 < $forward_exon_start ) {
					my $intron_start = $end + 1;
					my $intron_end   = $forward_exon_start - 1;
					$info =~ s/ID=exon/ID=intron/;
					print OUT
"$chr\t$source\t$intron\t$intron_start\t$intron_end\t$other1\t$strand\t$other2\t$info\n";
				}
			}
			$forward_exon_start = $start;
			$forward_exon_end   = $end;
			$forward_exon_info  = $info;

		}
		$forward_feature = $feature;
		$forward_line    = $_;
		if ( $feature ne "gene" ) {
			print OUT $_, "\n";
		}
	}
}

close GFF;
close OUT;
print $num;

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}
