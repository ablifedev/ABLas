#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;

my %opts;
GetOptions( \%opts, "gff=s","mod=s","br=s","sj=s","o=s","h" );

if (   !defined( $opts{sj} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{mod} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-sj          sj file                must be given;

		-gff         gff file               must be given;

		-br          boundary reads         must be given;

		-mod         gene model file        must be given;

		-o           out file               must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $sj_file    = $opts{sj};
my $gff_file    = $opts{gff};
my $br_file    = $opts{br};
my $mod_file    = $opts{mod};
my $out_file    = $opts{o};

my %gene_model = ();

#将gene model保存进哈希以便查询
open MOD, $mod_file || die;
while (<MOD>) {
	chomp;
	next if ( $_ =~ /chrom/ );
	my @line = split;
	$gene_model{$line[1]} = $line[0];
}
close MOD;

#读取sj文件
my %positive_sj=();
my %negative_sj=();
&load_all_sj($sj_file);

#读取boudary reads file
my %reads=();
open BR, $br_file || die;
while (<BR>) {
	chomp;
	next if ( $_ =~ /chrom/ );
	my @line = split(/\t/);
	$reads{$line[0]}{$line[1]}{"+"}=$line[2];
	$reads{$line[0]}{$line[1]}{"-"}=$line[3];
}
close BR;

###遍历gff，按基因检测AS事件
my (
	$gene, $isoform, $forward_chr, $forward_gene, $forward_isoform,
	$gene_count, $forward_start, $forward_end,  $exon_region_value, 
	$intron_region_value
  )
  = ( "", "", "", "", "", 1, 0, 0, 0, 0 );
my $forward_strand = "";
my %exon           = ();
my @exon_hash      = ();	#基因的exon区域
my @intron_hash    = ();	#基因的intron区域
my %sj_start	   = ();	#基因中所有起点的sj信息
my %sj_end  	   = ();	#基因中所有终点的sj信息
my %gene_region    = ();
my %as_count = (
	"IntronR"  => 0
);

#chr1    hg18_knownGene  mRNA    1310954 1323585 .       -       .       ID=CCNL2.1.5;Parent=CCNL2.1
#chr1    hg18_knownGene  exon    1323476 1323585 0.000000        -       .       ID=exon:CCNL2.1.5:1;Parent=CCNL2.1.5
#chr1    hg18_knownGene  intron  1320781 1323475 0.000000        -       .       ID=intron:CCNL2.1.5:2;Parent=CCNL2.1.5
#chr1    hg18_knownGene  exon    1320637 1320780 0.000000        -       .       ID=exon:CCNL2.1.5:2;Parent=CCNL2.1.5
open OUT, ">$out_file" || die;
print OUT "#chr\tAss1\tAss2\tstrand\tIntronR\tratio\tintron_ave_boudaryreads\tintron_sj\tisoform\n";

open GFF, $gff_file    || die;
while (<GFF>) {
	chomp;
	next if ( $_ =~ /^\#/ );
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	  )
	  = split /\t/;
	next if ( $feature =~ /chromosome/ || $feature =~ /gene/ || $feature =~ /protein/ );
	if($feature =~ /^mRNA|transcript$/){
		$info =~ m/ID=([\w\:\.\-]+);.*Parent=([\w\:\.\-]+)/;
		$isoform="$1";
		$gene="$2";
	}
	if($feature =~ /^exon$/){
		$info =~ m/Parent=([\w\:\.\-]+)/;
		$isoform="$1";
	}
	next if ( ( not defined( $gene_model{$gene} ) ) || ( $isoform ne $gene_model{$gene} ) );
	if ( $gene ne $forward_gene ) {
		##处理基因方法，检测AS事件
		&handling_gene( \%exon, \@exon_hash, \@intron_hash, \%gene_region, $forward_strand, $forward_chr );
		%exon        = ();
		%gene_region = ();
		%sj_start	 = ();
		%sj_end		 = ();
		@exon_hash      = ();
		@intron_hash    = ();
		$gene_count++;
		$exon_region_value   = 1;
		$intron_region_value = 0;
	}
	if ( $feature =~ /^exon$/ )
	{
		push @exon_hash,"$start:$end";
		if ( $exon_region_value > 1 ) {    #给intron区域赋值
			my $temp_start = 0;
			my $temp_end   = 0;
			if ( $strand eq "+" ) {
				$temp_start = $forward_end;
				$temp_end   = $start;
				$sj_start{$temp_start}=$temp_end;
				$sj_end{$temp_end}=$temp_start;
			}
			if ( $strand eq "-" ) {
				$temp_start = $end;
				$temp_end   = $forward_start;
				$sj_start{$temp_end}=$temp_start;
				$sj_end{$temp_start}=$temp_end;
			}
			push @intron_hash,"$temp_start:$temp_end";

			$gene_region{$intron_region_value} = "$temp_start:$temp_end";
		}
		$gene_region{$exon_region_value} = "$start:$end";

		$exon_region_value += 1;
		$intron_region_value -= 1;
		$forward_end = $end;
		$forward_start = $start;
		$forward_isoform = $isoform;
		$forward_strand  = $strand;
		$exon{$start}    = $end;	#将exon起始终止位点记录进入hash %exon
	}
	$forward_gene = $gene;
	$forward_chr = $chr;
}
close GFF;

close OUT;


# ###
# print "Done detect IntronR...\n";
# ###
# ###END detect IntronR

sub handling_gene{
	my ( $exon, $exon_hash, $intron_hash, $gene_region, $strand, $chr ) = @_;
	my $intron_num=0;
	foreach my $index (sort {$a<=>$b} keys %{$gene_region}) {
		next if $index>0;
		my ($intron_start,$intron_end) = split /\:/, $gene_region->{$index};
		my $sj_num=0;
		my $avg_boudary_reads = 0;
		if($strand eq "+"){
			$sj_num=$positive_sj{$chr}->{"$intron_start:$intron_end"} if defined($positive_sj{$chr}->{"$intron_start:$intron_end"});
			$avg_boudary_reads = ($reads{$chr}{$intron_start}{"+"}+$reads{$chr}{$intron_end}{"+"})/2;
			next if $reads{$chr}{$intron_start}{"+"}<2 or $reads{$chr}{$intron_end}{"+"}<2;
			next if $avg_boudary_reads < 3;
		}else{
			$sj_num=$negative_sj{$chr}->{"$intron_end:$intron_start"} if defined($negative_sj{$chr}->{"$intron_end:$intron_start"});
			$avg_boudary_reads = ($reads{$chr}{$intron_start}{"-"}+$reads{$chr}{$intron_end}{"-"})/2;
			next if $reads{$chr}{$intron_start}{"-"}<2 or $reads{$chr}{$intron_end}{"-"}<2;
			next if $avg_boudary_reads < 3;
		}
		next if ($sj_num+$avg_boudary_reads)<=0;
		my $ratio = $avg_boudary_reads/($sj_num+$avg_boudary_reads) ;
		# next if $ratio < 0.15;
		print OUT "$chr\t$intron_start\t$intron_end\t$strand\tIntronR\t$ratio\t$avg_boudary_reads\t$sj_num\t$forward_isoform\n";
	}

		
}


###子过程-----load_all_sj($sj_file)
sub load_all_sj {
	my ($sjFile) = @_;
	open SJ, $sjFile || die;
	while (<SJ>) {              #把所有的junction添加进hash
		#默认为读取通过tophat获取的junctions.bed文件
		#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
		#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
		chomp;
		next if ( $_ =~ /^track/i );
		my $start          = 0;
		my $end            = 0;
		my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
			$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
		my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
		if ( $strand eq "+" ) {
			$start = $left + $leftsize;
			$end   = $right - $rightsize + 1;
			$positive_sj{$chr}->{"$start:$end"} = $readsNum;
		}
		else {
			$start = $right - $rightsize + 1;
			$end   = $left + $leftsize;
			$negative_sj{$chr}->{"$start:$end"} = $readsNum;
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}


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
