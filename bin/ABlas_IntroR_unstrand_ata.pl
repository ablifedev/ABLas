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
GetOptions( \%opts, "bam=s", "fa=s","chrlen=s", "gff=s", "o=s", "anno=s", "od=s","chr=s","h" );

if (   !defined( $opts{bam} )
	|| !defined( $opts{chrlen} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| !defined( $opts{fa} )
	|| !defined( $opts{anno} )
	|| !defined( $opts{od} )
	|| !defined( $opts{chr} )
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

		-chr        chr              must be given;

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
my $thischr    = $opts{chr};

my $sam = Bio::DB::Sam->new(-fasta=>"$fa_file",
                                 -bam  =>"$bam_file");


my %gffexon=();
my $forward_chr="";

open GFF, $gff_file || die;
while (<GFF>) {
	chomp;
	next if(/\#/);
	my @line=split(/\t/,$_,6);
	next if $line[0] ne $thischr;
	# if($line[0] ne $forward_chr && $forward_chr ne ""){print "done $forward_chr\n";last;}
	if($line[2]=~/^exon$/){
		# print $_,"\n";
		for(my $i=$line[3];$i<=$line[4];$i++){
			$gffexon{$line[0]}->{$i} = 1;
		}
	}
	$forward_chr = $line[0];
}
close GFF;

###
print "Done record GFF exon...\n";
###

###detect IntronR

###遍历gff，按基因检测AS事件
my (
	$gene, $this_chr, $this_strand,
	$gene_count, $forward_start, $forward_end,  $exon_region_value, 
	$intron_region_value
  ) = ( "", "", "", 1, 0, 0, 0, 0 );
my %exon           = ();
my @exon_hash      = ();	#基因的exon区域
my @intron_hash    = ();	#基因的intron区域
my %sj_start	   = ();	#基因中所有起点的sj信息
my %sj_end  	   = ();	#基因中所有终点的sj信息
my %gene_region    = ();
my %as_count = (
	"IntronR"  => 0
);
my $mrna_start = 0;
my $mrna_end = 0;

#chr1    hg18_knownGene  mRNA    1310954 1323585 .       -       .       ID=CCNL2.1.5;Parent=CCNL2.1
#chr1    hg18_knownGene  exon    1323476 1323585 0.000000        -       .       ID=exon:CCNL2.1.5:1;Parent=CCNL2.1.5
#chr1    hg18_knownGene  intron  1320781 1323475 0.000000        -       .       ID=intron:CCNL2.1.5:2;Parent=CCNL2.1.5
#chr1    hg18_knownGene  exon    1320637 1320780 0.000000        -       .       ID=exon:CCNL2.1.5:2;Parent=CCNL2.1.5
open OUT, ">$out_file" || die;
open GFF, $gff_file    || die;
while (<GFF>) {
	chomp;
	next if ( $_ =~ /^\#/ );
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	  )
	  = split /\t+/;
	next if $chr ne $thischr;
	next if ( $feature =~ /^chromosome$/ || $_ =~ /^transposable_element_gene$/);
	if($feature =~ /^gene$/){$gene="";next;}
	if ( $feature =~ /^(mRNA|miRNA|mRNA_TE_gene|ncRNA|rRNA|snoRNA|snRNA|tRNA|transcript)$/ ) {
		##处理基因方法，检测AS事件
		# print "$this_chr\t$gene\t$this_strand\n";
		&handling_gene( \%exon, \@exon_hash, \@intron_hash, \%gene_region, $this_strand, $this_chr, $gene, $mrna_start, $mrna_end) if($gene ne "");
		$info =~ m/^ID=([\w|\.|\%|\:|\_|\-|\s]+);/;
		$gene = $1;
		$this_chr = $chr;
		$this_strand = $strand;
		%exon        = ();
		%gene_region = ();
		%sj_start	 = ();
		%sj_end		 = ();
		@exon_hash      = ();
		@intron_hash    = ();
		$gene_count++;
		$exon_region_value   = 1;
		$intron_region_value = 0;
		$mrna_start = $start;
		$mrna_end = $end;
	}
	if ( $feature =~ /exon/ )
	{
		next if $gene eq "";
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
		$exon{$start}    = $end;	#将exon起始终止位点记录进入hash %exon
	}
}
close GFF;

close OUT;

###
print "Done detect IntronR...\n";
###
###END detect IntronR

sub handling_gene{
	my ( $exon, $exon_hash, $intron_hash, $gene_region, $strand, $chr, $gene, $genestart, $geneend) = @_;
	my $sum_intron_len = 0;
	my $sum_intron_base = 0;
	my $ave_base_all = 0;

	my %base = ();
	print $genestart, "\t",$geneend, "\n";
	my @alignments = $sam->get_features_by_location(-seq_id => $chr,
                                                            -start  => $genestart,
                                                            -end    => $geneend);
	for my $a (@alignments) {
		my $start  = $a->start;
		my $end    = $a->end;
		my $cigar     = $a->cigar_str;
		next if $cigar =~ /^\d+M\d+N\d+M$/;
		# print "$chr\t$start\t$end\t$strand\n";
		for(my $i=$start;$i<=$end;$i++){
			$base{$chr}->{$i} = 0 if not defined($base{$chr}->{$i});
			$base{$chr}->{$i} ++; 
		}
	}


	# foreach my $index (sort {$a<=>$b} keys %{$gene_region}) {
	# 	next if $index>0;
	# 	my ($intron_start,$intron_end) = split /\:/, $gene_region->{$index};
	# 	$sum_intron_len += $intron_end - $intron_start + 1;
	# 	$sum_intron_base += &sum($chr, $intron_start, $intron_end);
	# }
	# $ave_base_all = $sum_intron_base/$sum_intron_len if $sum_intron_len!=0;
	# print "sum_intron_base:$sum_intron_base\tsum_intron_len:$sum_intron_len\tave_base_all:$ave_base_all\n";

	LABEL1:foreach my $index (sort {$a<=>$b} keys %{$gene_region}) {
		next if $index>0;
		my $ends_exon_len = 0;
		my $ends_exon_base = 0;
		my $ends_exon_ave_base = 0; 

		my ($intron_start,$intron_end) = split /\:/, $gene_region->{$index};
		print "intron_start:$intron_start\tintron_end:$intron_end\n";

		for(my $j=$intron_start+1;$j<$intron_end;$j++){
			if(defined($gffexon{$chr}->{$j}) && $gffexon{$chr}->{$j}==1){
				# print "next\n";
				next LABEL1;
			}
		}

		my $intron_len = $intron_end - $intron_start - 1;
		my $intron_base = &sum($chr, $intron_start+1, $intron_end-1,\%base);
		my $intron_max_base = &max_depth($chr, $intron_start+1, $intron_end-1,\%base);
		my $intron_mid_base = &mid_depth($chr, $intron_start+1, $intron_end-1,\%base);
		next LABEL1 if $intron_len==0;
		my $intron_ave_base = $intron_base/$intron_len;
		print "intron_max_base:$intron_max_base\tintron_ave_base:$intron_ave_base\tintron_base:$intron_base\tintron_len:$intron_len\tintron_mid_depth:$intron_mid_base\n";
		next LABEL1 if $intron_base < 100;
		next LABEL1 if $intron_mid_base < 5;
		# print $intron_base,"\t",$intron_max_base,"\n";
		# next LABEL1 if $intron_ave_base < $ave_base_all * 2;		##筛选条件1：intron区base平均depth要大于总intron区base平均depth的2倍.

		##筛选条件2：左右两边border reads至少有一边大于1。
		my $dd = 0;
		for(my $i=$intron_start;$i<=$intron_start+4;$i++){
			# next LABEL1 if &sum($chr, $i, $i)<1;
			if(&sum($chr, $i, $i,\%base)<1){
				$dd++;
				last;
			}
		}
		for(my $i=$intron_end-4;$i<=$intron_end;$i++){
			# next LABEL1 if &sum($chr, $i, $i)<1;
			if(&sum($chr, $i, $i,\%base)<1){
				$dd++;
				last;
			}
		}
		next LABEL1 if $dd==2;
		print "have border reads\n";

		##筛选条件3：intron区base平均depth要大于左右两边exon区base平均depth的15%.
		my ($left_exon_start,$left_exon_end) = split /\:/, $gene_region->{0-$index};
		$ends_exon_len += $left_exon_end - $left_exon_start + 1 ;
		$ends_exon_base += &sum($chr, $left_exon_start, $left_exon_end,\%base);
		my ($right_exon_start,$right_exon_end) = split /\:/, $gene_region->{1-$index};
		$ends_exon_len += $right_exon_end - $right_exon_start + 1 ;
		$ends_exon_base += &sum($chr, $right_exon_start, $right_exon_end,\%base);
		$ends_exon_ave_base = $ends_exon_base/$ends_exon_len;
		print "ends_exon_ave_base:$ends_exon_ave_base\tintron_ave_base:$intron_ave_base\n";
		next LABEL1 if($intron_mid_base < $ends_exon_ave_base * 0.15 && $intron_ave_base < 15);
		print "find a intronR\n";

		my $idnum = 0;
		$as_count{"IntronR"}++;
		$idnum = $as_count{"IntronR"};
		my $asid = "IntronR\_$idnum";

		if($strand eq "+"){
			my $structure = "1\^2\-\,0";
			my $splice_chain = "$intron_start\^$intron_end\-\,";
			my $splice_position = "$intron_start,$intron_end";
			my $info = "flanks \"$left_exon_start\-\,$right_exon_end\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
			print OUT "$chr\t$gene\tIR\t$left_exon_start\t$right_exon_end\t\.\t$strand\t\.\t$info\n";
		}
		if($strand eq "-"){
			my $structure = "1\^2\-\,0";
			my $splice_chain = "$intron_end\^$intron_start\-\,";
			my $splice_position = "$intron_end,$intron_start";
			my $info = "flanks \"$left_exon_end\-\,$right_exon_start\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
			print OUT "$chr\t$gene\tIR\t$left_exon_end\t$right_exon_start\t\.\t$strand\t\.\t$info\n";
		}
	}
	%base = ();
}

#####方法：取最大值
sub max_depth {
	my ( $chr, $start, $end, $base) = @_;
	my @basenum = ();
	for ( my $j = $start ; $j <= $end ; $j++ ) {
		$base->{$chr}->{$j} = 0 if not defined($base->{$chr}->{$j});
		push @basenum, $base->{$chr}->{$j};
	}
	return 0 unless (@basenum);
	my $max_value = shift @basenum;
	foreach (@basenum) {
		if ( $_ > $max_value ) {
			$max_value = $_;
		}
	}
	return $max_value ;
}

#####方法：取base总和
sub sum {
	my ( $chr, $start, $end, $base) = @_;
	my @basenum = ();
	for ( my $j = $start ; $j <= $end ; $j++ ) {
		$base->{$chr}->{$j} = 0 if not defined($base->{$chr}->{$j});
		push @basenum, $base->{$chr}->{$j};
	}

	return 0 unless (@basenum);
	my $sum_value = shift @basenum;
	foreach (@basenum) {
			$sum_value += $_;
	}
	return $sum_value;
}

#####方法：取中位数（参数为起始结束位点）
sub mid_depth {
	my ( $chr, $start, $end , $base) = @_;
	my @basenum = ();
	for ( my $j = $start ; $j <= $end ; $j++ ) {
		$base->{$chr}->{$j} = 0 if not defined($base->{$chr}->{$j});
		push @basenum, $base->{$chr}->{$j};
	}
	return -1 unless (@basenum);
	my @list = sort { $a <=> $b } @basenum;
	# print join(":",@basenum),"\n";
	# print join(":",@list),"\n";
	my $count = @list;
	if ( $count == 0 ) {
		return -1;
	}
	if ( ( $count % 2 ) == 1 ) {
		return $list[ int( ( $count - 1 ) / 2 ) ];
	}
	elsif ( ( $count % 2 ) == 0 ) {
		return (
			$list[ int( ( $count - 1 ) / 2 ) ] + $list[ int( ($count) / 2 ) ] )
		  / 2;
	}
}

#####方法：取中位数（参数为起始结束位点）
sub mid25_depth {
	my ( $chr, $start, $end , $base) = @_;
	my @basenum = ();
	for ( my $j = $start ; $j <= $end ; $j++ ) {
		$base->{$chr}->{$j} = 0 if not defined($base->{$chr}->{$j});
		push @basenum, $base->{$chr}->{$j};
	}
	return -1 unless (@basenum);
	my @list = sort { $a <=> $b } @basenum;
	# print join(":",@basenum),"\n";
	# print join(":",@list),"\n";
	my $count = @list;
	if ( $count == 0 ) {
		return -1;
	}
	if ( ( $count % 4 ) == 0 ) {
		return $list[ int( $count / 4 ) ];
	}
	elsif ( ( $count % 4 ) == 1 ) {
		return $list[ int( ( $count - 1 ) / 4 ) ];
	}elsif ( ( $count % 4 ) == 2 ) {
		return $list[ int( ( $count - 2 ) / 4 ) ];
	}elsif ( ( $count % 4 ) == 3 ) {
		return $list[ int( ( $count - 3 ) / 4 ) ];
	}
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
