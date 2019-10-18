#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;

my %opts;
GetOptions( \%opts, "sam=s", "chrlen=s", "gff=s", "o=s", "anno=s", "od=s","h" );

if (   !defined( $opts{sam} )
	|| !defined( $opts{chrlen} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| !defined( $opts{anno} )
	|| !defined( $opts{od} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-sam         splice junction        must be given;
		
		-chrlen      chromsome lenght       must be given;

		-gff         gff file               must be given;

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

my $sam_file    = $opts{sam};
my $chrlen_file = $opts{chrlen};
my $gff_file    = $opts{gff};
my $out_file    = $opts{o};
my $anno_file    = $opts{anno};
my $outdir    = $opts{od};

my $group       = 0;
my $count       = 0;
my $len         = 0;

my %base_p = ();    #positive
my %base_n = ();    #negative

open CHRLEN, $chrlen_file || die;
while (<CHRLEN>) {

	#chr1	201139347
	chomp;
	my ( $chr, $len ) = split /\s+/;
	$group = int( $len / 100000000 );
	$count = int( $len % 100000000 );

	#	print $group, "\t", $count, "\n";
	if ( $group == 0 ) {
		$base_p{$chr} = "";
		vec( $base_p{$chr}, $len, 16 ) = 0;
		$base_n{$chr} = "";
		vec( $base_n{$chr}, $len, 16 ) = 0;

	}
	if ( $group > 0 ) {
		$base_p{$chr} = "";
		vec( $base_p{$chr}, 100000000, 16 ) = 0;
		$base_n{$chr} = "";
		vec( $base_n{$chr}, 100000000, 16 ) = 0;

		for ( my $j = 1 ; $j < $group ; $j++ ) {
			my $newchr = $chr . "_" . $j;
			$base_p{$newchr} = "";
			vec( $base_p{$newchr},100000000, 16 ) = 0;
			$base_n{$newchr} = "";
			vec( $base_n{$newchr},100000000, 16 ) = 0;

		}
		my $newchr = $chr . "_" . $group;
		$base_p{$newchr} = "";
		vec( $base_p{$newchr},$count, 16 ) = 0;
		$base_n{$newchr} = "";
		vec( $base_n{$newchr},$count, 16 ) = 0;
	}
}
close CHRLEN;

###
print "Done creat vec...\n";
###

open SAM, $sam_file    || die;

my $num = 0;
#####将合适的reads添加进hash
while (<SAM>) {
	$num++;
	# last if $num > 1000000;

#SRR533613.8839386	0	chr1	3041921	255	50M	*	0	0	TCTCCCTTTCCACTCCCTGTCTGTTACTACATTTCTAAGAAGAAGGCAGG	hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhZhhhhhhhhhhhWhhhhhhh	NM:i:0	NH:i:1
	chomp;
	next if (/^\@/);
	my ( $strand, $chr, $start, $other, $feature ) = ( split(/\t/,$_,7) )[ 1 .. 5 ];
	next if ( $feature !~ /^\d+M$/ && $feature !~ /^\d+M\d+N\d+M$/ );    #非普通mapped reads或junctions reads则跳过
	my $read_strand = $strand & (0x0010) ? "-" : "+";
	my $left_block  = 0;
	my $middle_block= 0;
	my $right_block = 0;
	##记录base数
	if($read_strand eq "+"){    #正向
		if($feature =~ m/^(\d+)M$/){
			$len = $1;
			##记录base数
			for ( my $i = $start ; $i < $start + $len ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_p{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_p{$newchr}, $c, 16 ) += 1;
				}
			}
		}elsif($feature =~ /^(\d+)M(\d+)N(\d+)M$/){
			$left_block = $1;
			$middle_block = $2;
			$right_block = $3;
			for ( my $i = $start ; $i < $start + $left_block ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_p{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_p{$newchr}, $c, 16 ) += 1;
				}
			}
			for ( my $i = $start + $left_block + $middle_block ; $i < $start + $left_block + $middle_block + $right_block ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_p{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_p{$newchr}, $c, 16 ) += 1;
				}
			}

		}
	}elsif($read_strand eq "-"){    #负向
		if($feature =~ m/^(\d+)M$/){
			$len = $1;
			##记录base数
			for ( my $i = $start ; $i < $start + $len ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_n{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_n{$newchr}, $c, 16 ) += 1;
				}
			}
		}elsif($feature =~ /^(\d+)M(\d+)N(\d+)M$/){
			$left_block = $1;
			$middle_block = $2;
			$right_block = $3;
			for ( my $i = $start ; $i < $start + $left_block ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_n{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_n{$newchr}, $c, 16 ) += 1;
				}
			}
			for ( my $i = $start + $left_block + $middle_block ; $i < $start + $left_block + $middle_block + $right_block ; $i++ ) {
				my $g = int( $i / 100000000 );
				my $c = int( $i % 100000000 );
				if ( $g == 0 ) {
					vec( $base_n{$chr}, $i, 16 ) += 1;
				}
				elsif ( $g > 0 ) {
					my $newchr = $chr . "_" . $g;
					vec( $base_n{$newchr}, $c, 16 ) += 1;
				}
			}
		}
	}

}

close SAM;

###
print "Done reading SAM...\n";
###

###do EXP

my $outfile = $outdir."/"."EXP_exon";
open EXON, ">$outfile" || die;
print EXON "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
$outfile = $outdir."/"."EXP_mRNA";
open MRNA, ">$outfile" || die;
print MRNA "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
$outfile = $outdir."/"."EXP_five";
open FIVE, ">$outfile" || die;
print FIVE "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
$outfile = $outdir."/"."EXP_three";
open THREE, ">$outfile" || die;
print THREE "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
$outfile = $outdir."/"."EXP_intron";
open INTRON, ">$outfile" || die;
print INTRON "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";


#my $reads_num = 0;
my $base_num  = 0;

open ANNO, $anno_file    || die;
#NC_008467.1     mRNA       +       36663   37331   POPTRDRAFT_750402

while (<ANNO>) {
	chomp;
	my ($chr,$feature,$strand,$start,$end) = split;
#	$reads_num = &sumofreads($chr,$start,$end);
	$base_num  = &sum($chr,$start,$end,$strand);
	
	print EXON "$chr\t$feature\t$strand\t$start\t$end\t$base_num\n" if $feature eq "exon";
	print MRNA "$chr\t$feature\t$strand\t$start\t$end\t$base_num\n" if $feature eq "mRNA";
	print INTRON "$chr\t$feature\t$strand\t$start\t$end\t$base_num\n" if $feature eq "intron";
	print THREE "$chr\t$feature\t$strand\t$start\t$end\t$base_num\n" if $feature eq "three";
	print FIVE "$chr\t$feature\t$strand\t$start\t$end\t$base_num\n" if $feature eq "five";

}


close ANNO;
close EXON;
close FIVE;
close THREE;
close INTRON;

###
print "Done do EXP...\n";
###
###END do EXP 


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
my $flag = 0;
my %as_count = (
	"IntronR"  => 0
);

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
	next if ( $feature =~ /chromosome/ || $_ =~ /match/ || $_ =~ /region/);
	if($feature =~ /^gene$/){
		$flag = 0;
		next;
	}
	if ( $feature =~ /^(mRNA|miRNA|mRNA_TE_gene|ncRNA|rRNA|snoRNA|snRNA|tRNA|transcript|pseudogenic_transcript)$/ ) {
		##处理基因方法，检测AS事件
		&handling_gene( \%exon, \@exon_hash, \@intron_hash, \%gene_region, $this_strand, $this_chr, $gene ) if($gene ne "");
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
		$flag = 1;
	}
	if ( $feature =~ /exon/ )
	{
		next if $flag == 0;
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
	my ( $exon, $exon_hash, $intron_hash, $gene_region, $strand, $chr, $gene) = @_;
	my $sum_intron_len = 0;
	my $sum_intron_base = 0;
	my $ave_base_all = 0;
	# my $intron_num = 0;
	# foreach my $index (sort {$a<=>$b} keys %{$gene_region}) {
	# 	next if $index>0;
	# 	my ($intron_start,$intron_end) = split /\:/, $gene_region->{$index};
	# 	$sum_intron_len += $intron_end - $intron_start - 1;
	# 	$sum_intron_base += &sum($chr, $intron_start+1, $intron_end-1, $strand);
	# 	$intron_num++;
	# }
	# $ave_base_all = $sum_intron_base/$sum_intron_len if $sum_intron_len!=0;
	print $gene, "\n";

	LABEL1:foreach my $index (sort {$a<=>$b} keys %{$gene_region}) {
		next if $index>0;
		my $ends_exon_len = 0;
		my $ends_exon_base = 0;
		my $ends_exon_ave_base = 0; 

		my ($intron_start,$intron_end) = split /\:/, $gene_region->{$index};
		print "intron_start:$intron_start\tintron_end:$intron_end\n";

		my $intron_len = $intron_end - $intron_start - 1;
		next if $intron_len<=0;
		my $intron_base = &sum($chr, $intron_start+1, $intron_end-1, $strand);
		
		my $intron_mid_base = &mid_depth($chr, $intron_start+1, $intron_end-1,$strand);
		my $intron_ave_base = $intron_base/$intron_len;
		print "intron_ave_base:$intron_ave_base\tintron_base:$intron_base\tintron_len:$intron_len\tintron_mid_depth:$intron_mid_base\n";
		next LABEL1 if $intron_base < 100;
		next LABEL1 if $intron_mid_base < 5;
		
		# next LABEL1 if($intron_ave_base < $ave_base_all * 2 && $intron_num>=2);		##筛选条件1：intron区base平均depth要大于总intron区base平均depth的2倍.
		# print $ave_base_all,"\t",$intron_ave_base,"\t",$intron_base,"\t",$intron_len,"\n";

		##筛选条件2：左右两边border reads都要大于1。
		my $dd = 0;
		for(my $i=$intron_start;$i<=$intron_start+4;$i++){
			# next LABEL1 if &sum($chr, $i, $i)<1;
			if(&sum($chr, $i, $i,$strand)<1){
				$dd++;
				last;
			}
		}
		for(my $i=$intron_end-4;$i<=$intron_end;$i++){
			# next LABEL1 if &sum($chr, $i, $i)<1;
			if(&sum($chr, $i, $i,$strand)<1){
				$dd++;
				last;
			}
		}
		next LABEL1 if $dd==2;
		print "have border reads\n";

		##筛选条件3：intron区base平均depth要大于左右两边exon区base平均depth的15%.
		my ($left_exon_start,$left_exon_end) = split /\:/, $gene_region->{0-$index};
		$ends_exon_len += $left_exon_end - $left_exon_start + 1 ;
		$ends_exon_base += &sum($chr, $left_exon_start, $left_exon_end, $strand);
		my ($right_exon_start,$right_exon_end) = split /\:/, $gene_region->{1-$index};
		$ends_exon_len += $right_exon_end - $right_exon_start + 1 ;
		$ends_exon_base += &sum($chr, $right_exon_start, $right_exon_end, $strand);
		$ends_exon_ave_base = $ends_exon_base/$ends_exon_len;
		print "ends_exon_ave_base:$ends_exon_ave_base\tintron_ave_base:$intron_ave_base\n";
		next LABEL1 if $intron_ave_base < $ends_exon_ave_base * 0.15;
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
}

#####方法：取base总和
sub sum {
	my ( $chr, $start, $end, $strand ) = @_;
	my @basenum = ();
	if($strand eq "+"){
		for ( my $j = $start ; $j <= $end ; $j++ ) {
			my $g = int( $j / 100000000 );
			my $c = int( $j % 100000000 );
			if ( $g == 0 ) {
				push @basenum, int( vec( $base_p{$chr}, $j, 16 ) );
			}
			elsif ( $g > 0 ) {
				my $newchr = $chr . "_" . $g;
				push @basenum, int( vec( $base_p{$newchr}, $c, 16 ) );
			}
		}
	}
	if($strand eq "-"){
		for ( my $j = $start ; $j <= $end ; $j++ ) {
			my $g = int( $j / 100000000 );
			my $c = int( $j % 100000000 );
			if ( $g == 0 ) {
				push @basenum, int( vec( $base_n{$chr}, $j, 16 ) );
			}
			elsif ( $g > 0 ) {
				my $newchr = $chr . "_" . $g;
				push @basenum, int( vec( $base_n{$newchr}, $c, 16 ) );
			}
		}
	}

	return 0 unless (@basenum);
	my $sum_value = shift @basenum;
	foreach (@basenum) {
			$sum_value += $_;
	}
	return $sum_value;
}

#####方法：取base中位数
sub mid_depth {
	my ( $chr, $start, $end, $strand ) = @_;
	my @basenum = ();
	if($strand eq "+"){
		for ( my $j = $start ; $j <= $end ; $j++ ) {
			my $g = int( $j / 100000000 );
			my $c = int( $j % 100000000 );
			if ( $g == 0 ) {
				push @basenum, int( vec( $base_p{$chr}, $j, 16 ) );
			}
			elsif ( $g > 0 ) {
				my $newchr = $chr . "_" . $g;
				push @basenum, int( vec( $base_p{$newchr}, $c, 16 ) );
			}
		}
	}
	if($strand eq "-"){
		for ( my $j = $start ; $j <= $end ; $j++ ) {
			my $g = int( $j / 100000000 );
			my $c = int( $j % 100000000 );
			if ( $g == 0 ) {
				push @basenum, int( vec( $base_n{$chr}, $j, 16 ) );
			}
			elsif ( $g > 0 ) {
				my $newchr = $chr . "_" . $g;
				push @basenum, int( vec( $base_n{$newchr}, $c, 16 ) );
			}
		}
	}

	return 0 unless (@basenum);
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
