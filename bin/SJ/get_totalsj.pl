#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#获取几个样本中总共找到的sj。
#输出结果格式：chr1	870846	871694	-，分别为：染色体	junction_start	junction_end	strand

my %opts;
GetOptions( \%opts, "bed1=s", "bed2=s", "bed3=s", "bed4=s", "bed5=s","bed6=s", "o=s",
	"h" );

if (   !defined( $opts{bed1} )
#	|| !defined( $opts{bed2} )
#	|| !defined( $opts{bed3} )
#	|| !defined( $opts{bed4} )
#	|| !defined( $opts{bed5} )
#	|| !defined( $opts{bed6} )
	|| !defined( $opts{o} ) )
{
	print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-bed1	bed file	must be given

			-bed2	bed file	must be given

			-bed3	bed file	must be given
			
			-bed4	bed file	must be given
			
			-bed5	bed file	must be given
			
			-bed6	bed file	must be given

			-o		outfile     must be given

	Usage End.

	exit;
}

#############Time_start#############
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $bed_file1 = $opts{bed1};
my $bed_file2 = $opts{bed2};
my $bed_file3 = $opts{bed3};
my $bed_file4 = $opts{bed4};
my $bed_file5 = $opts{bed5};
my $bed_file6 = $opts{bed6};
my $out_file  = $opts{o};

my $key               = "";
my %negative_junction = ();
my %positive_junction = ();

open BED1, $bed_file1   || die;
open OUT,  ">$out_file" || die;
while (<BED1>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED1;

open BED2, $bed_file2   || die;
while (<BED2>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED2;

open BED3, $bed_file3   || die;
while (<BED3>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED3;

open BED4, $bed_file4   || die;
while (<BED4>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED4;

open BED5, $bed_file5   || die;
while (<BED5>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED5;

open BED6, $bed_file6   || die;
while (<BED6>) {    #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;

	if ( $strand eq "+" ) {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $positive_junction{$key} ) ) {
			$positive_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
	else {
		$key = $chr . ":" . $start . ":" . $end;
		if ( !defined( $negative_junction{$key} ) ) {
			$negative_junction{$key} = 1;
			print OUT "$chr\t$start\t$end\t$strand\n";
		}
	}
}
close BED6;

close OUT;

############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";
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

