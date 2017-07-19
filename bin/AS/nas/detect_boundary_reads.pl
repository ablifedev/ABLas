#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $version="1.0.2";


# 程序功能：找到所有gene model的exon boundary reads，统计其数目
# 实现思路：记录所有gene model exon的 边界位点，
# 并按位置分block存储，读取sam时对应到相应block找到是否落在边界且满足boundaryreads的条件，做统计。
# boundary reads的条件：intron一端至少有4个碱基。

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$fGFF,$model,$edge,$BLOCK_SIZE,$SEARCH_DEPTH);
GetOptions(
				"help|?" =>\&USAGE,
				"g:s"=>\$fGFF,
				"i:s"=>\$Index,
				"m:s"=>\$model,
				"edge:i"=>\$edge,
				"BLOCK_SIZE:i"=>\$BLOCK_SIZE,
				"SEARCH_DEPTH:i"=>\$SEARCH_DEPTH,
				) or &USAGE;
&USAGE unless ($Index and $fGFF) ;

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################
my $BEGIN_TIME=time();
my $current_dir = `pwd`;chomp($current_dir);

$fGFF=&AbsolutePath("file",$fGFF);
#$Index=&AbsolutePath("file",$Index);
#my $TMR=&GetTMR("$Index.statMappedInfo");
my $TMR= 0;
#print $TMR,"\n";
#######################################################################################

$edge||=4 ;
$BLOCK_SIZE=1000;
$SEARCH_DEPTH||=30;

my %gene_model = ();

#将gene model保存进哈希以便查询
open MOD, $model || die;
while (<MOD>) {
	chomp;
	next if ( $_ =~ /^chrom/ );
	my @line = split;
	$gene_model{$line[0]} = 1;
}
close MOD;

my %Gff;
&load_gff($fGFF,\%Gff);

my %ClusterGff;
&clusterGene(\%Gff,$BLOCK_SIZE,\%ClusterGff);

my %reads;
&getsam($Index,$BLOCK_SIZE,\%ClusterGff,\%reads);
print "\nDone. getCover\n",time()-$BEGIN_TIME,"s\n";

&Output_arrReads(\%Gff,\%reads,$Index);

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

# ------------------------------------------------------------------
# Get gene info form GFF file.
# ------------------------------------------------------------------

sub load_gff{
	my ($infile,$hash_ref)=@_;
	my $transcript_count=0;
	my $gene_id="";
	my $transcript="";
	open (IN ,"<$infile") || die "Can't open $!";
	$/="\n";
	while (<IN>) { ###open gff file
		chomp;
		next if ($_!~/^[\w\.]+/) ;
		my @line=split(/\s+/,$_);
		if($line[2]=~/^exon$/){
			$line[-1] =~ m/Parent=([\w\:\.\-]+)/;
			my $isoform="$1";
			next if not defined( $gene_model{$isoform} );
			$hash_ref->{$line[0]}{$line[3]}="s";
			$hash_ref->{$line[0]}{$line[4]}="e";
		}
	}
	close IN;
}

# ------------------------------------------------------------------
# Cluster Genes in order to save time later.
# ------------------------------------------------------------------
#sub clusterGene{
#	#We cluster the genes based by the start position of the genes;
#	my ($refGff,$SIZE,$refCluster)=@_;
#	foreach my $chro (keys %{$refGff}) {
#		foreach my $gene (keys %{$refGff->{$chro}}) {
#			my $blockNum=int($refGff->{$chro}{$gene}->[0]->[2]/$SIZE);
#			$refCluster->{$chro}{$blockNum}{$gene}=$refGff->{$chro}{$gene};
#		}
#	}
#}

sub clusterGene{
	#We cluster the genes based by the start position of the genes;
	my ($refGff,$SIZE,$refCluster)=@_;
	foreach my $chro (keys %{$refGff}) {
		foreach my $pos (keys %{$refGff->{$chro}}) {
			my $blockNum=int($pos/$SIZE);
			$refCluster->{$chro}{$blockNum}{$pos}=$refGff->{$chro}{$pos};
		}
	}
}

# ------------------------------------------------------------------
# Get the reads and calculate the coverage of genes.
# ------------------------------------------------------------------
#SBS_0005:6:26:19464:16954#0     0       Chr1    3623    255     50M     *       0       0       CTCCCTCCAAATTATTAGATATACCAAACCAGAGAAAACAAATACATAAT      hhhhhhhhhhhhhhhhghhhhhhhgh
#SBS_0005:6:86:18478:14172#0     0       Chr1    3631    255     50M     *       0       0       AAATTATTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAA      `Z\_`Z_]W_Z_\_W_caY\ccccR^
#SBS_0005:6:93:6946:18247#0      0       Chr1    3643    255     50M     *       0       0       ATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGATTACAG      hhhhhhhhhhhhhhhhhhhhhhhhhh

sub getsam{
	my ($Key,$SIZE,$refCluster,$reads_ref)=@_;
	open (IN,"<","$Key") or die $!;
	while (<IN>) {
		next if (/^@/) ;
		chomp;
		my @line=split(/\t+/,$_);
		next if ( $line[5] !~ /^\d+M$/ );
		$TMR++;
		my ($End,$Strand);
		my $d=$line[1] &(0x0010);
		if ($d ==0) {
			$Strand="+";
		}
		elsif($d==16){
			$Strand="-";
		}
		&belongGen($line[3],$line[5],$Strand,$line[2],$SIZE,$refCluster,$reads_ref);
	}
	print "Total mapped reads count:$TMR\n";
	close (IN);
}

# ------------------------------------------------------------------
# Find which gene the reads belong to.
# ------------------------------------------------------------------
sub belongGen{
	my ($start,$match,$strand,$chro,$SIZE,$refCluster,$reads)=@_;
	my $blockNum_start=int($start/$SIZE);
	my $end;
#	if ($chro!~/^[Cc]hr/) {
#		$chro = "chr".$chro;
#	}
	if ($match=~/^(\d+)M$/){
		$end=$start+$1-1;
	}
	#This reads may belong to this block or the prior N block
	my $blockNum_end=int($end/$SIZE);
	if ($blockNum_start==$blockNum_end) {
		if (exists $refCluster->{$chro}{$blockNum_start}) {
			foreach my $position (sort keys %{$refCluster->{$chro}{$blockNum_start}}) {
				if ($position>=$start and $position<=$end){
					if($refCluster->{$chro}{$blockNum_start}{$position} eq "s"){
						if($position-$start >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							# print $start,"\t",$end,"\t",$position,"\t","s","\n";
							return 0;
						}
					}else{
						if($end-$position >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							# print $start,"\t",$end,"\t",$position,"\t","e","\n";
							return 0;
						}
					}
				}
			}
		}
	}elsif(exists $refCluster->{$chro}{$blockNum_end} && ! exists $refCluster->{$chro}{$blockNum_start}){
		foreach my $position (sort keys %{$refCluster->{$chro}{$blockNum_end}}) {
				if ($position>=$start and $position<=$end){
					if($refCluster->{$chro}{$blockNum_end}{$position} eq "s"){
						if($position-$start >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}else{
						if($end-$position >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}
				}
			}
	}
	elsif(exists $refCluster->{$chro}{$blockNum_start} && !exists $refCluster->{$chro}{$blockNum_end}){
		foreach my $position (sort keys %{$refCluster->{$chro}{$blockNum_start}}) {
				if ($position>=$start and $position<=$end){
					if($refCluster->{$chro}{$blockNum_start}{$position} eq "s"){
						if($position-$start >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}else{
						if($end-$position >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}
				}
			}
	}elsif(exists $refCluster->{$chro}{$blockNum_start} && exists $refCluster->{$chro}{$blockNum_end}){
		foreach my $position (sort keys %{$refCluster->{$chro}{$blockNum_start}}) {
				if ($position>=$start and $position<=$end){
					if($refCluster->{$chro}{$blockNum_start}{$position} eq "s"){
						if($position-$start >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}else{
						if($end-$position >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}
				}
			}
		foreach my $position (sort keys %{$refCluster->{$chro}{$blockNum_end}}) {
				if ($position>=$start and $position<=$end){
					if($refCluster->{$chro}{$blockNum_end}{$position} eq "s"){
						if($position-$start >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}else{
						if($end-$position >= $edge){
							$reads->{$chro}{$position}{$strand}++;
							return 0;
						}
					}
				}
			}
	}else{
		return 0;
	}
	
}

# ------------------------------------------------------------------
# Output the result
# ------------------------------------------------------------------
sub Output_arrReads{
	my ($refCluster,$reads_ref,$Key)=@_;
	open (OUT,">","$Key.boundary_reads_stat") or die $!;
#	print OUT  "Chro\tPos\t+ reads\t- reads\n";
	foreach my $chro (sort {$a cmp $b} keys %{$refCluster}) {
		foreach my $pos (sort {$a cmp $b} keys %{$refCluster->{$chro}}) {
			$reads_ref->{$chro}{$pos}{"+"}=0 if not defined($reads_ref->{$chro}{$pos}{"+"});
			$reads_ref->{$chro}{$pos}{"-"}=0 if not defined($reads_ref->{$chro}{$pos}{"-"});
			my $positive = $reads_ref->{$chro}{$pos}{"+"};
			my $negative = $reads_ref->{$chro}{$pos}{"-"};
			print OUT $chro,"\t",$pos,"\t",$positive,"\t",$negative,"\n";
		}
		
	}
	close OUT;
}

sub cat
#function:quit redundance
#input:($para,@array), para is the merge length
#output:(@array),
#for example (0,1,3,4,7,5,8)->(1,3,4,8) (1,1,3,4,7,5,8)->(1,8)
{
	my($merge,@input) = @_;
	my $i = 0;
	my @output = ();
	my %hash = ();
	my $each = 0;
	my $begin = "";
	my $end = 0;
	for ($i=0;$i<@input;$i+=2){
		my $Qb = $input[$i];
		my $Qe = $input[$i+1];
		if($Qb > $Qe) { next; }
		if(defined($hash{$Qb}))	{
			if($hash{$Qb} < $Qe) {
				$hash{$Qb} = $Qe;
			}
		}
		else {
			$hash{$Qb} = $Qe;
		}
		$Qb = 0;
	}

	foreach $each (sort {$a <=> $b} keys %hash){
		if($begin eq ""){
			$begin = $each;
			$end = $hash{$each};
		}
		else{
			if($hash{$each} > $end){
				if($each > $end + $merge){
					push(@output,$begin);
					push(@output,$end);
					$begin = $each;
					$end = $hash{$each};
				}
				else {
					$end = $hash{$each};
				}
			}
		}
	}
	if(keys %hash > 0){
		push(@output,$begin);
		push(@output,$end);
	}
	%hash = ();
	return(@output);
}

sub AbsolutePath {		#获取指定目录或文件的决定路径
	my ($type,$input) = @_;
	my $return;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	}
	elsif($type eq 'file'){
		my $pwd = `pwd`;
		chomp $pwd;
		my $dir=dirname($input);
		my $file=basename($input);
		chdir($dir);
		$return = `pwd`;
		chomp $return;
		$return .="\/".$file;
		chdir($pwd);
	}
	return $return;
}

sub overlap_indentify {
	my ($start1,$end1,$start2,$end2)=@_;
	my %hash=();
	for (my $i=$start1;$i<=$end1;$i++) {
		$hash{$i}++;
	}
	for (my $j=$start2;$j<=$end2;$j++) {
		$hash{$j}++;
	}
	my @overlap=grep {$hash{$_}==2} keys(%hash);
	undef %hash;
	return scalar(@overlap);
}


#sub intersection_union{
#	my ($a_ref,$b_ref,$inter,$union)=@_;
#	my %a = map{$_ => 1} @{$a_ref};
#	my %b = map{$_ => 1} @{$b_ref};
#	my @{$inter} = grep {$a{$_}} @{$b_ref};
#	my %merge = map {$_ => 1} @{$a_ref},@{$b_ref};
#	my @{$union} = keys (%merge);
#}

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: XiongXu<xiongxu\@ablife.cc> <xuxiong19880610\@163.com>

Example:perl ~/work/perl/RNA-seq/filter_out_gene_reads.pl -g ~/test/Homo_sapiens/hg18_knowgene.gff -i accepted_hits.sam

Description:
	This program is used to filter out the reads mapped to the gene region.

Usage:
  -g               <file>   Gff File, forced;
  -i               <str>    Index of inFiles and outFiles. Index.rmap, Index.geneExpression.xls, Index.randCheck
  -m               <str>    gene model
  -BLOCK_SIZE      <int>    Block Size, used to save time, default 1000;
  -SEARCH_DEPTH    <int>    Search Depth, used to save time with BLOCK_SIZE, default 30;


USAGE
	print $usage;
	exit;
}
