#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;

my $buf = "-";
my %opts;
GetOptions( \%opts, "gff=s", "sj=s", "o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff        annotation file        must be given;

		-sj         splice junction        must be given;

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
my $sj_file  = $opts{sj};
my $out_file = $opts{o};


#调试信息
print "1_END;", "\n";
my %negative_start = ();    #负向junction起始位点信息记录hash
my %positive_start = ();    #正向同上
my %negative_end   = ();    #负向junction结束位点信息记录hash
my %positive_end   = ();    #正向同上
my $start          = 0;
my $end            = 0;
open SJ, $sj_file || die;

while (<SJ>) {              #把所有的junction添加进hash

	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
	#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
	chomp;
	next if ( $_ =~ /^track/i );
	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
		$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;
	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
	if ( $strand eq "+" ) {
		$start = $left + $leftsize;
		$end   = $right - $rightsize + 1;
		if ( !defined( $positive_start{$chr}->{$start} ) ) {
			$positive_start{$chr}->{$start} = "$end\:$readsNum";
		}
		else {
			$positive_start{$chr}->{$start} .= "\t$end\:$readsNum";

		}
		if ( !defined( $positive_end{$chr}->{$end} ) ) {
			$positive_end{$chr}->{$end} = "$start\:$readsNum";
		}
		else {
			$positive_end{$chr}->{$end} .= "\t$start\:$readsNum";

			#调试信息
			#print $positive_end{$chr}->{$end},"\n";
		}
	}
	else {
		$start = $right - $rightsize + 1;
		$end   = $left + $leftsize;
		if ( !defined( $negative_start{$chr}->{$start} ) ) {
			$negative_start{$chr}->{$start} = "$end\:$readsNum";
		}
		else {
			$negative_start{$chr}->{$start} .= "\t$end\:$readsNum";
		}
		if ( !defined( $negative_end{$chr}->{$end} ) ) {
			$negative_end{$chr}->{$end} = "$start\:$readsNum";
		}
		else {
			$negative_end{$chr}->{$end} .= "\t$start\:$readsNum";
		}
	}
}
close SJ;

#调试信息
print "1. done reading SJ ...........", "\n\n";

#my $debug_count = 0;	### debug

my (
	$gene, $isoform, $forward_gene, $forward_isoform,
	$gene_count, $forward_end,  $exon_region_value, $intron_region_value
  )
  = ( "", "", "", "", 1, 0, 0, 0 );
my $forward_strand = "";
my %exon           = ();
my %exon_hash      = ();
my %gene_region    = ();
my %alt_splice     = ();
my %intron         = ();
my %normal_sj      = ();
open OUT, ">$out_file" || die;
open GFF, $gff_file    || die;
# chr4_gl000193_random    hg19_knownGene  transcript      49163   88375   .       -       .       ID=uc003izx.3;Parent=uc003izx
while (<GFF>) {
	chomp;
	next if ( $_ =~ /^\#/ );
	my (
		$chr,    $source, $feature, $start, $end,
		$other1, $strand, $other2,  $info
	  )
	  = split /\t/;
	next if ( $feature =~ /chromosome/ || $_ =~ /gene/ );
	$info =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/;
	$gene    = $2;
	$isoform = "$2.$3";
	if ( $gene ne $forward_gene ) {
		&handling_gene( \%exon, \%gene_region, $forward_strand, $chr );
		%exon        = ();
		%gene_region = ();
		undef(%exon);
		undef(%gene_region);
		$gene_count++;
		$exon_region_value   = 1;
		$intron_region_value = -1;
	}
	if (   ( $feature =~ /exon/ )
		&& ( defined( $gene_model{$gene} ) )
		&& ( $isoform eq $gene_model{$gene} ) )
	{
		$exon_hash{$gene} = "" if ( !defined( $exon_hash{$gene} ) );
		$exon_hash{$gene} .= "$start:$end\n";
		if ( $exon_region_value > 1 ) {    #给intron区域赋值
			my $temp_start = $forward_end;
			my $temp_end   = $start;
			if ( $strand eq "-" ) {
				$temp_start = $end;
				$temp_end   = $forward_end;
			}
			$intron{$gene} = "" if ( !defined( $intron{$gene} ) );
			$intron{$gene} .= "$temp_start:$temp_end\n";

			#调试代码
			#print $intron{$gene},"\n";
			for ( my $i = $temp_start + 1 ; $i < $temp_end ; $i++ ) {
				$gene_region{$i} = $intron_region_value;
			}
		}
		for ( my $i = $start ; $i <= $end ; $i++ ) {    #给exon区域赋值
			$gene_region{$i} = $exon_region_value;
		}
		$exon_region_value += 1;
		$intron_region_value -= 1;
		if ( $strand eq "+" ) {
			$forward_end = $end;
		}
		else {
			$forward_end = $start;
		}
		$forward_isoform = $isoform;
		$exon{$start}    = $end;
		$forward_strand  = $strand;
	}
	$forward_gene = $gene;

}
close GFF;

#调试信息
print "3_END;", "\n";
my @a3ss_array = ();
my @a5ss_array = ();
foreach my $Chr ( sort keys %alt_splice ) {
	foreach my $as ( sort { $a <=> $b } keys %{ $alt_splice{$Chr} } ) {
		foreach my $line ( split /\n/, $alt_splice{$Chr}->{$as} ) {
			if ( $line =~ /a3ss/i ) {
				push @a3ss_array, $line . "\n";
			}
			elsif ( $line =~ /a5ss/i ) {
				push @a5ss_array, $line . "\n";
			}
			else {
				print OUT $line . "\n";
			}
		}
	}
}
undef(%alt_splice);

#调试信息
print "4_END;", "\n";
###############
my $compare = 0;
foreach my $a3ss (@a3ss_array) {
	foreach my $a5ss (@a5ss_array) {
		next if ( ( split /\s+/, $a3ss )[8] ne ( split /\s+/, $a5ss )[8] );
		my $flag = &find_MXE( $a3ss, $a5ss, ( split /\s+/, $a3ss )[8] );
		my ( $chr, $start, $end, $strand, $type, $percent, $curr, $major,
			$gene ) = ( split /\s+/, $a3ss )[ 0 .. 8 ];
		if ( $flag > 1 ) {
			my $tmp_major = ( split /\s+/, $a5ss )[7] + $major;
			my $tmp_curr  = ( split /\s+/, $a5ss )[6] + $curr;
			print OUT "$chr\t$start\t$end\t$strand\tMXE\|",
			  ( split /\s+/, $a5ss )[1], "\|";
			$normal_sj{$chr}->{"$start:$end"} = 1;
			print OUT ( split /\s+/, $a5ss )[2], "\t",
			  ( int( $tmp_curr / ( $tmp_curr + $tmp_major ) * 100 ) ) / 100;
			print OUT "\t$tmp_curr\t$tmp_major\t$gene\n";
			$a3ss .= "\t#";
			$a5ss .= "\t#";
		}
		elsif ( $flag == 1 ) {
			my ( $start1, $end1 ) = ( split /\s+/, $a5ss )[ 1 .. 2 ];
			if ( $start > $start1 ) {
				my $tmp = $start;
				$start = $end1;
				$end   = $tmp;
			}
			else {
				$start = $end;
				$end   = $start1;
			}
			my $tmp_major = ( split /\s+/, $a5ss )[7] + $major;
			my $tmp_curr  = ( split /\s+/, $a5ss )[6] + $curr;
			print OUT "$chr\t$start\t$end\t$strand\tcassetteExon", "\t",
			  ( int( $tmp_curr / ( $tmp_curr + $tmp_major ) * 100 ) ) / 100;
			print OUT "\t$tmp_curr\t$tmp_major\t$gene\n";
			$a3ss .= "\t#";
			$a5ss .= "\t#";
		}
	}
}

LABEL1: foreach (@a3ss_array) {
	next if (/\#$/);
	my ( $chr, $start, $end, $strand, $type, $percent, $curr, $major, $isoform )
	  = ( split /\s+/ )[ 0 .. 8 ];
	$isoform =~ s/\.\d+$//;
  LABEL2: foreach my $tmp ( split /\n/, $exon_hash{$isoform} ) {
		my ( $start1, $end1 ) = ( split /\:/, $tmp );
		if ( $start < $start1 && $end > $end1 ) {
			my $line = $_;
			$line =~ s/(A3SS)/$1\&ES\|$start1:$end1/;
			print OUT $line;
			next LABEL1;
		}
	}
	print OUT $_;

}
LABEL1: foreach (@a5ss_array) {
	next if (/\#$/);
	my ( $chr, $start, $end, $strand, $type, $percent, $curr, $major, $isoform )
	  = ( split /\s+/ )[ 0 .. 8 ];
	$isoform =~ s/\.\d+$//;
  LABEL2: foreach my $tmp ( split /\n/, $exon_hash{$isoform} ) {
		my ( $start1, $end1 ) = ( split /\:/, $tmp );
		if ( $start < $start1 && $end > $end1 ) {
			my $line = $_;
			$line =~ s/(A5SS)/$1\&ES\|$start1:$end1/;
			print OUT $line;
			next LABEL1;
		}
	}
	print OUT $_;
}

#调试信息
print "5_END;", "\n";

#将已注释的SJ添加进哈希
#open SJ, $anno_sj || die;
#while (<SJ>) {
#
#	#chr1	870846	871694	(4)[36_4](4/0)	50	-	870846	871694	255,0,0	2	50,50	0,798
#	chomp;
#	# next if ( $_ !~ /^c/i );
#	my ( $chr, $left, $right, $feature, $length, $strand ) =
#	  ( split /\s+/ )[ 0 .. 5 ];
#	$start                            = $left + 50;
#	$end                              = $right - 49;
#	$normal_sj{$chr}->{"$start:$end"} = 1;
#}
#close SJ;
#
##找出未能定义的SJ
#my %other_sj = ();
#open SJ, $sj_file || die;
#while (<SJ>) {
#
#	#chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
#	#chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
#	chomp;
#	next if ( $_ =~ /^track/i );
#	my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd, 
#		$itemRgb, $blockCount, $blockSizes, $blockStarts) = ( split /\s+/ );
#	my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];
#	$start = $left + $leftsize;
#	$end   = $right - $rightsize + 1;
#	if ( !defined( $normal_sj{$chr}->{"$start:$end"} ) ) {
#		$other_sj{$chr}->{$start} = "$end:$strand:$readsNum";
#	}
#}
#close SJ;
#
##调试信息
#print "6_END;", "\n";
#my $t_end     = 0;
#my $t_strand  = 0;
#my $reads_num = 0;
#
##找出基因里的other SJ
#open GFF, $gff_file || die;
#while (<GFF>) {
#	chomp;
#	my (
#		$chr,    $source, $feature, $start, $end,
#		$other1, $strand, $other2,  $info
#	  )
#	  = split /\s+/;
#	next if ( $feature ne "gene" );
#	$info =~ m/^ID=(\w+:)?(\S+)\.(\d+)(:\w+)?;.+$/;
#	$gene = $2;
#	for ( my $i = $start ; $i <= $end ; $i++ ) {
#		if ( defined( $other_sj{$chr}->{$i} ) ) {
#			( $t_end, $t_strand, $reads_num ) = split /\:/,
#			  $other_sj{$chr}->{$i};
#			if ( $strand eq $t_strand ) {
#				print OUT
#				  "$chr\t$i\t$t_end\t$strand\tother\t1\t$reads_num\t0\t$gene\n";
#			}
#		}
#	}
#}

close OUT;

#调试信息
print "7_END;", "\n";

sub find_MXE {
	my ( $forward_line, $line, $isoform ) = @_;
	my ( $start1, $end1, $strand1 ) = ( split /\s+/, $forward_line )[ 1 .. 3 ];
	my ( $start2, $end2, $strand2 ) = ( split /\s+/, $line )[ 1 .. 3 ];
	my $flag = 0;
	$isoform =~ s/\.\d+$//;
	if (   ( $start1 > $end2 && $strand1 eq "-" )
		|| ( $start2 > $end1 && $strand1 eq "+" ) )
	{
		my $f;
		my $r;
		if ( $start1 > $end2 ) {
			$r = $start1;
			$f = $end2;
		}
		else {
			$r = $start2;
			$f = $end1;
		}
		foreach ( split /\n/, $intron{$isoform} ) {
			my ( $start, $end ) = split /\:/;
			$flag = 1 if ( $start < $f && $end > $r );
		}
		if ($flag) {

			#print "start1:$start1\tend1:$end1\tstart2:$start2\tend2:$end2\n";
			foreach ( split /\n/, $exon_hash{$isoform} ) {
				my ( $start, $end ) = split /\:/;
				$flag++
				  if ( ( $start >= $start1 - 1 && $end <= $end1 + 1 )
					|| ( $start >= $start2 - 1 && $end <= $end2 + 1 ) );
			}
		}
	}
	return $flag;
}

sub handling_gene {
	my ( $exon, $gene_region, $strand, $chr ) = @_;
	my $num = scalar( keys %{$exon} );
	my ( $count, $start, $end, $forward, $next ) = ( 0, 0, 0, 0, 0 );
	my @sort_keys      = ();
	my @start_junction = ();
	my @end_junction   = ();
	if ( $strand eq "+" ) {
		@sort_keys = sort { $a <=> $b } keys %{$exon};
	}
	else {
		@sort_keys = sort { $b <=> $a } keys %{$exon};
	}
	foreach my $point (@sort_keys) {
		if ( $strand eq "-" ) {
			$end   = $point;
			$start = $exon->{$end};
		}
		else {
			$start = $point;
			$end   = $exon->{$start};
		}
		if ( $strand eq "+" ) {
			@start_junction = split /\s+/, $positive_end{$chr}->{$start}
			  if ( defined( $positive_end{$chr}->{$start} ) )
			  ;    #当前exon起始位点确定的所有junction
			@end_junction = split /\s+/, $positive_start{$chr}->{$end}
			  if ( defined( $positive_start{$chr}->{$end} ) )
			  ;    #结束位点确定的所有junction
		}
		else {
			@start_junction = split /\s+/, $negative_end{$chr}->{$start}
			  if ( defined( $negative_end{$chr}->{$start} ) );
			@end_junction = split /\s+/, $negative_start{$chr}->{$end}
			  if ( defined( $negative_start{$chr}->{$end} ) );
		}
		if( @start_junction < 1 && @end_junction < 1 ){
			$count++;
			next;
		}
		if ( $count == 1 ) {    #检测互斥5'UTR，通过第二个exon
			&handling_utr( $chr, $start, \@start_junction, $gene_region,
				$strand, 5 )
			  if ( @start_junction > 0 );
			&handling_normal( $chr, $end, \@end_junction, $gene_region, $strand,
				3 )
			  if ( @end_junction > 0 );
		}
		elsif ( $count == $num - 2 )
		{                       #通过倒数第二个exon检测互斥3'UTR
			&handling_utr( $chr, $end, \@end_junction, $gene_region, $strand,
				3 )
			  if ( @end_junction > 0 );
			&handling_normal( $chr, $start, \@start_junction, $gene_region,
				$strand, 5 )
			  if ( @start_junction > 0 );
		}
		else {                  #检测其他的exon
			&handling_normal( $chr, $start, \@start_junction, $gene_region,
				$strand, 5 )
			  if ( @start_junction > 0 );
			&handling_normal( $chr, $end, \@end_junction, $gene_region, $strand,
				3 )
			  if ( @end_junction > 0 );
		}
		$count++;
		undef(@start_junction);
		undef(@end_junction);
	}
}

sub handling_utr {
	my ( $chr, $current, $junction, $gene_region, $strand, $utr_type ) = @_;
	my $sum     = 0;
	my $as_type = 3;
	foreach my $i ( @{$junction} ) {
		my ( $position, $reads ) = split /\:/, $i;
			$sum += $reads;
	}
	$as_type = 5 if ( $utr_type != $as_type );
	my @as_tmp      = ();
	my $major_reads = 0;
	foreach my $i ( @{$junction} ) {    #处理相同3'端的junction
		my ( $position, $reads ) = split /\:/, $i;
		my $junc_num = $reads;
		if (   !defined( $gene_region->{ $position + 1 } )
			&& !defined( $gene_region->{ $position - 1 } ) )
		{                               #含有互斥的UTR
			if ( $current > $position ) {
				push @as_tmp,
				  "$chr\t$position\t$current\t$strand\t$utr_type" . "pMXE\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$position:$current"} = 1;
			}
			else {
				push @as_tmp,
				  "$chr\t$current\t$position\t$strand\t$utr_type" . "pMXE\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$current:$position"} = 1;
			}
		}
		elsif (
			(
				   $gene_region->{ $position - 1 } > 0
				&& $gene_region->{ $position + 1 } < 0
			)
			|| (   $gene_region->{ $position - 1 } < 0
				&& $gene_region->{ $position + 1 } > 0 )
		  )
		{    #主要的剪切形式，同genemodel
			$major_reads = $junc_num;
		}
		else {    #否则为可变的3'/5'剪切位点
			if ( $current > $position ) {
				push @as_tmp,
				  "$chr\t$position\t$current\t$strand\tA$as_type" . "SS\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$position:$current"} = 1;
			}
			else {
				push @as_tmp,
				  "$chr\t$current\t$position\t$strand\tA$as_type" . "SS\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$current:$position"} = 1;
			}
		}
	}

	foreach my $i (@as_tmp) {
		my ( $f, $r ) = ( split /\s+/, $i )[ 1 .. 2 ];
		if ( !defined( $alt_splice{$chr}{$f} ) ) {
			$alt_splice{$chr}{$f} = "$i\t$major_reads\t$forward_isoform\n";
		}
		else {
			my $temp = "";
			my $flag = 0;
			foreach my $k ( split /\n/, $alt_splice{$chr}{$f} ) {
				if ( $flag == 0 && $r < ( split /\s+/, $k )[2] ) {
					$temp .= "$i\t$major_reads\t$forward_isoform\n";
					$flag = 1;
				}
				$flag = 1 if ( $r == ( split /\s+/, $k )[2] );
				$temp .= $k . "\n";
			}
			$temp .= "$i\t$major_reads\t$forward_isoform\n" if ( $flag == 0 );
			$alt_splice{$chr}{$f} = $temp;
			undef($temp);
			undef($flag);
		}
	}
	undef(@as_tmp);
}

sub handling_normal {
	my ( $chr, $current, $junction, $gene_region, $strand, $direction ) = @_;
	my $sum     = 0;
	my $as_type = 3;
	foreach my $i ( @{$junction} ) {
		my ( $position, $reads ) = split /\:/, $i;
		next
		  if ( !defined( $gene_region->{ $position + 1 } )
			|| !defined( $gene_region->{ $position - 1 } ) );
		$sum += $reads;
	}
	$as_type = 5 if ( $direction != $as_type );
	my @as_tmp      = ();
	my $major_reads = 0;
	foreach my $i ( @{$junction} ) {    #处理相同3'端的junction
		my ( $position, $reads ) = split /\:/, $i;
		my $junc_num = $reads;
		if (   !defined( $gene_region->{ $position + 1 } )
			|| !defined( $gene_region->{ $position - 1 } ) )
		{                               #跳过未知junction
			next;
		}
		elsif (
			(
				   $gene_region->{ $position - 1 } > 0
				&& $gene_region->{ $position + 1 } < 0
			)
			|| (   $gene_region->{ $position - 1 } < 0
				&& $gene_region->{ $position + 1 } > 0 )
		  )
		{                               #主要的剪切形式，同genemodel
			if (
				abs( $gene_region->{$current} - $gene_region->{$position} ) ==
				1 )
			{
				$major_reads = $junc_num;
			}
			else {
				if ( $current > $position ) {
					push @as_tmp,
					  "$chr\t$position\t$current\t$strand\tES\t"
					  . ( int( $junc_num / $sum * 100 ) / 100 )
					  . "\t$junc_num";
					$normal_sj{$chr}->{"$position:$current"} = 1;
				}
				else {
					push @as_tmp,
					  "$chr\t$current\t$position\t$strand\tES\t"
					  . ( int( $junc_num / $sum * 100 ) / 100 )
					  . "\t$junc_num";
					$normal_sj{$chr}->{"$current:$position"} = 1;
				}
			}
		}
		else {    #否则为可变的3'/5'剪切位点
			if ( $current > $position ) {
				push @as_tmp,
				  "$chr\t$position\t$current\t$strand\tA$as_type" . "SS\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$position:$current"} = 1;
			}
			else {
				push @as_tmp,
				  "$chr\t$current\t$position\t$strand\tA$as_type" . "SS\t"
				  . ( int( $junc_num / $sum * 100 ) / 100 )
				  . "\t$junc_num";
				$normal_sj{$chr}->{"$current:$position"} = 1;
			}
		}
	}
	foreach my $i (@as_tmp) {
		my ( $f, $r ) = ( split /\s+/, $i )[ 1 .. 2 ];
		if ( !defined( $alt_splice{$chr}{$f} ) ) {
			$alt_splice{$chr}{$f} = "$i\t$major_reads\t$forward_isoform\n";
		}
		else {
			my $temp = "";
			my $flag = 0;
			foreach my $k ( split /\n/, $alt_splice{$chr}{$f} ) {
				if ( $flag == 0 && $r < ( split /\s+/, $k )[2] ) {
					$temp .= "$i\t$major_reads\t$forward_isoform\n";
					$flag = 1;
				}
				$flag = 1 if ( $r == ( split /\s+/, $k )[2] );
				$temp .= $k . "\n";
			}
			$temp .= "$i\t$major_reads\t$forward_isoform\n" if ( $flag == 0 );
			$alt_splice{$chr}{$f} = $temp;
			undef($temp);
			undef($flag);
		}
	}
	undef(@as_tmp);
}

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