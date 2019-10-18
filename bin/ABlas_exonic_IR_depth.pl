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
GetOptions( \%opts, "bam=s", "fa=s", "gff=s", "o=s", "h" );

if (   !defined( $opts{bam} )
	|| !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| !defined( $opts{fa} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-bam         bam file        must be given;
		
		-gff         gff file               must be given;

		-fa         fa file               must be given;

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

my $bam_file    = $opts{bam};
my $fa_file    = $opts{fa};
my $gff_file    = $opts{gff};
my $out_file    = $opts{o};

my $sam = Bio::DB::Sam->new(-fasta=>"$fa_file",
                                 -bam  =>"$bam_file");


open OUT,">$out_file";

open GFF, $gff_file || die;
while (<GFF>) {
	chomp;
	my @line=split(/\t/);
	$line[5]="intron_base,intron_ave_base,intron_mid_base" if(/\#/);
	print OUT join("\t",@line[0..8]),"\talt_ratio\n" if(/\#/);
	next if(/\#/);
	
	# NW_001867363.1	rna14997,rna14998	IR	5773801	5776872	.	+	20	flanks "5773801-,5776872^"; structure "1^2-,0"; splice_chain "5773874^5774094-,"; splice_position "5773874,5774094"; ID "IR_3";
	$line[8]=~/splice_position\s\"(\d+),(\d+)\"\;/;
	my $start = $1;
	my $end = $2;
	if($start>$end){my $temp=$start;$start=$end;$end=$temp;}
	# print $start,"\t",$end,"\n";
	$line[5] = &handling_sj($line[0],$start,$end,$line[6]);
	# print OUT $_;
	# print OUT "\t",$depth,"\n";
	my ($total_base,$avg_base,$mid_base)=split(/,/,$line[5]);
	next if ($avg_base+$line[7]<=0);
	my $alt_ratio = $avg_base/($avg_base+$line[7]);
	print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";

}
close GFF;


close OUT;

###
print "Done detect IntronR...\n";
###
###END detect IntronR

sub handling_sj{
	my ( $Chr, $Start, $End, $Strand) = @_;
	my %base = ();
	# print $genestart, "\t",$geneend, "\n";
	my @alignments = $sam->get_features_by_location(-seq_id => $Chr,
                                                            -start  => $Start,
                                                            -end    => $End);
	for my $a (@alignments) {
		my $start  = $a->start;
		my $end    = $a->end;
		my $reads_strand = $a->strand;
		if($reads_strand==1){
			$reads_strand="+";
		}else{
			$reads_strand="-";
		}
		my $cigar     = $a->cigar_str;
		# print "$chr\t$start\t$end\t$reads_strand\t$cigar\n";
		next if $reads_strand ne $Strand;
		next if (($cigar !~ /^\d+M\d+N\d+M$/) && ($cigar !~ /^\d+M$/));

		if($cigar =~ m/^(\d+)M$/){
			my $len = $1;
			##记录base数
			for ( my $i = $start ; $i < $start + $len ; $i++ ) {
				$base{$Chr}->{$i} = 0 if not defined($base{$Chr}->{$i});
				$base{$Chr}->{$i} ++;
			}
		}elsif($cigar =~ /^(\d+)M(\d+)N(\d+)M$/){
			my $left_block = $1;
			my $middle_block = $2;
			my $right_block = $3;
			for ( my $i = $start ; $i < $start + $left_block ; $i++ ) {
				$base{$Chr}->{$i} = 0 if not defined($base{$Chr}->{$i});
				$base{$Chr}->{$i} ++;
			}
			for ( my $i = $start + $left_block + $middle_block ; $i < $start + $left_block + $middle_block + $right_block ; $i++ ) {
				$base{$Chr}->{$i} = 0 if not defined($base{$Chr}->{$i});
				$base{$Chr}->{$i} ++;
			}

		}
	}

	my $intron_len = $End - $Start - 1;
	return 0 if $intron_len<=0;
	my $intron_base = &sum($Chr, $Start+1, $End-1,\%base);
	my $intron_mid_base = &mid_depth($Chr, $Start+1, $End-1,\%base);
	my $intron_ave_base = sprintf("%.2f",$intron_base/$intron_len);
	undef(%base);

	return "$intron_base\,$intron_ave_base\,$intron_mid_base";
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

#####方法：取base中位数
sub mid_depth {
	my ( $chr, $start, $end, $base ) = @_;
	my @basenum = ();
	for ( my $j = $start ; $j <= $end ; $j++ ) {
		$base->{$chr}->{$j} = 0 if not defined($base->{$chr}->{$j});
		push @basenum, $base->{$chr}->{$j};
	}


	return 0 unless (@basenum);
	my @list = sort { $a <=> $b } @basenum;
	# print join(":",@basenum),"\n";
	# print join(":",@list),"\n";
	my $count = @list;
	if ( $count == 0 ) {
		return 0;
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
