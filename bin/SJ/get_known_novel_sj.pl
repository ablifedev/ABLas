#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#从一个样本的junction.bed文件中提取已知和新的junctions。

my %opts;
GetOptions( \%opts, "asj=s", "allsj=s", "o=s", "o2=s", "h" );

if (!defined( $opts{asj} )
 || !defined( $opts{allsj} )
 || !defined( $opts{o2} )
 || !defined( $opts{o} ) )
{
 print <<"	Usage End.";
	
		Description:This programme is used for ~
		
				Version: $ver
		
		Usage:perl $0
		
			-asj           gff_sj_file              must be given

			-allsj		   tophat_sj               must be given

			-o             novelsjfile          must be given
			
			-o2            knownsjfile          must be given

	Usage End.

 exit;
}

my $anno_sj  = $opts{asj};
my $all_sj   = $opts{allsj};
my $outfile  = $opts{o};
my $outfile2 = $opts{o2};

my %normal_sj = ();

#将已注释的SJ添加进哈希
open ASJ, $anno_sj || die;
while (<ASJ>) {
	chomp;
	next if ( $_ =~ /^track/i );
	my $start          = 0;
	my $end            = 0;
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	$start = $left + $leftsize;
	$end   = $right - $rightsize + 1;
	my $key = "$chr:$start:$end:$strand";
	$normal_sj{$key} = 1;
}
close ASJ;

#将SpliceMap跑出的的所有SJ与%normal_sj对比
open ALLSJ, $all_sj      || die;
open OUT,   ">$outfile"  || die;
open OUT2,  ">$outfile2" || die;
while (<ALLSJ>) {

 	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	# print $blockSizes,"\n";
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	my $start = $left + $leftsize;
	my $end   = $right - $rightsize + 1;
	my $key = "$chr:$start:$end:$strand";
 	if ( !defined( $normal_sj{$key} ) ) {
 		#新的sj，格式chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
  		print OUT $_, "\n";
 	}
 	else {
 		#已知的sj，格式chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
 		print OUT2 $_, "\n";
 	}
}
close ALLSJ;
close OUT;
close OUT2
