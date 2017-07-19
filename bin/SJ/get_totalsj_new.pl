#!/usr/bin/perl -w

my $ver = "1.0.0";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's
#explanation,Meanwhile,annotation your programme in English if possible.

#从多个样本总共找到的sj文件中(由gettotalsj.pl获取)提取新的sj和已知的sj。
#输出结果格式：chr1	870846	871694	-，分别为：染色体	junction_start	junction_end	strand

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
		
			-asj           gffsjfile              must be given

			-allsj		   all_sj               must be given

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

 #Supercontig_1.1 39330   39566   (2)[15_2](2/0)  50      -
 chomp;
 my ( $chr, $left, $right, $feature, $length, $strand ) =
   ( split /\s+/ )[ 0 .. 5 ];
 my $start = $left + 50;
 my $end   = $right - 49;
 my $key = "$chr:$start:$end:$strand";
 $normal_sj{$key} = 1;
}
close ASJ;

#将SpliceMap跑出的的所有SJ与%normal_sj对比
open ALLSJ, $all_sj      || die;
open OUT,   ">$outfile"  || die;
open OUT2,  ">$outfile2" || die;
while (<ALLSJ>) {

 	#chr1	870846	871694	-
	#chr1	870846	871694	+
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $start, $end, $strand) = split /\s+/ ;
	my $key = "$chr:$start:$end:$strand";
 	if ( !defined( $normal_sj{$key} ) ) {
  		print OUT $_, "\n";
 	}
 	else {
 		print OUT2 $_, "\n";
 	}
}
close ALLSJ;
close OUT;
close OUT2
