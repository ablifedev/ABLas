#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $buf = "-";
my %opts;
GetOptions( \%opts, "gff=s", "t=s", "sj=s","o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{t} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff        annotation file        must be given;

		-t          as type [SE]           must be given;

		-sj         sj                     must be given;

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
my $as_type = $opts{t};
my $sj_file = $opts{sj};
my $out_file = $opts{o};
open OUT, ">$out_file" || die;
print OUT "#Chromosome\tIsoform\tType\tStart\tEnd\tAlt_sj_reads\tStrand\tModel_sj_reads\tFeature\tAlt_ratio\n";

my %positive_sj    = ();    #正向sj的reads数
my %negative_sj    = ();    #负向同上
&load_all_sj($sj_file);

if($as_type eq "SE"){
	&load_SE_gff($gff_file);
}
if($as_type eq "CE"){
	&load_CE_gff($gff_file);
}
if($as_type eq "MXE"){
	&load_MXE_gff($gff_file);
}
if($as_type eq "AFE"){
	&load_AFE_gff($gff_file);
}
if($as_type eq "ALE"){
	&load_ALE_gff($gff_file);
}
if($as_type eq "II"){
	&load_II_gff($gff_file);
}
if($as_type eq "EE"){
	&load_EE_gff($gff_file);
}
if($as_type eq "A3SS"){
	&load_A3SS_gff($gff_file);
}
if($as_type eq "A5SS"){
	&load_A5SS_gff($gff_file);
}
if($as_type eq "IR"){
	&load_IR_gff($gff_file);
}
close OUT;


# chr1	uc009viu.1	SE	8776	8232	.	-	1	flanks "8776^,8232-"; structure "0,1-2^"; splice_chain ",8417-8364^"; skipped_exon_no "1"; flank_exon_start "9622,8232"; flank_exon_end "8776,8131"; skipped_exon_start "8417"; skipped_exon_end "8364"; splice_position "8417,8364";
#子过程-----load_SE_gff($gff_file)
sub load_SE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$flank_start:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$flank_end"});
			my $alt_sj_reads = $positive_sj{$chr}->{"$flank_start:$flank_end"};
			my $model_sj_reads = "";
			for(my $i=0;$i<=$#pos;$i++){
				if($i==0){
					$positive_sj{$chr}->{"$flank_start:$pos[$i]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[$i]"});
					print "$chr\t+\t$flank_start\t$pos[$i]\t",$positive_sj{$chr}->{"$flank_start:$pos[$i]"},"\n" if $positive_sj{$chr}->{"$flank_start:$pos[$i]"} != 0 ;
					$model_sj_reads = $positive_sj{$chr}->{"$flank_start:$pos[$i]"};
				}elsif($i==$#pos){
					$positive_sj{$chr}->{"$pos[$i]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[$i]:$flank_end"});
					$model_sj_reads .= ",".($positive_sj{$chr}->{"$pos[$i]:$flank_end"});
				}else{
					$positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"} = 0 if not defined($positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$model_sj_reads .= ",".($positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$i++;
				}
			}
			$line[7]=$model_sj_reads;
			$line[5]=$alt_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$flank_start:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$flank_end"});
			my $alt_sj_reads = $negative_sj{$chr}->{"$flank_start:$flank_end"};
			my $model_sj_reads = "";
			for(my $i=0;$i<=$#pos;$i++){
				if($i==0){
					$negative_sj{$chr}->{"$flank_start:$pos[$i]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[$i]"});
					$model_sj_reads = $negative_sj{$chr}->{"$flank_start:$pos[$i]"};
				}elsif($i==$#pos){
					$negative_sj{$chr}->{"$pos[$i]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[$i]:$flank_end"});
					$model_sj_reads .= ",".($negative_sj{$chr}->{"$pos[$i]:$flank_end"});
				}else{
					$negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"} = 0 if not defined($negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$model_sj_reads .= ",".($negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$i++;
				}
			}
			$line[7]=$model_sj_reads;
			$line[5]=$alt_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_CE_gff($gff_file)
sub load_CE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$flank_start:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$flank_end"});
			my $model_sj_reads = $positive_sj{$chr}->{"$flank_start:$flank_end"};
			my $alt_sj_reads = "";
			for(my $i=0;$i<=$#pos;$i++){
				if($i==0){
					$positive_sj{$chr}->{"$flank_start:$pos[$i]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[$i]"});
					$alt_sj_reads = $positive_sj{$chr}->{"$flank_start:$pos[$i]"};
				}elsif($i==$#pos){
					$positive_sj{$chr}->{"$pos[$i]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[$i]:$flank_end"});
					$alt_sj_reads .= ",".($positive_sj{$chr}->{"$pos[$i]:$flank_end"});
				}else{
					$positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"} = 0 if not defined($positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$alt_sj_reads .= ",".($positive_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$i++;
				}
			}
			$line[5]=$alt_sj_reads;
			$line[7]=$model_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$flank_start:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$flank_end"});
			my $model_sj_reads = $negative_sj{$chr}->{"$flank_start:$flank_end"};
			my $alt_sj_reads = "";
			for(my $i=0;$i<=$#pos;$i++){
				if($i==0){
					$negative_sj{$chr}->{"$flank_start:$pos[$i]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[$i]"});
					$alt_sj_reads = $negative_sj{$chr}->{"$flank_start:$pos[$i]"};
				}elsif($i==$#pos){
					$negative_sj{$chr}->{"$pos[$i]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[$i]:$flank_end"});
					$alt_sj_reads .= ",".($negative_sj{$chr}->{"$pos[$i]:$flank_end"});
				}else{
					$negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"} = 0 if not defined($negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$alt_sj_reads .= ",".($negative_sj{$chr}->{"$pos[$i]:$pos[$i+1]"});
					$i++;
				}
			}
			$line[5]=$alt_sj_reads;
			$line[7]=$model_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_MXE_gff($gff_file)
sub load_MXE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[0]"});
			$positive_sj{$chr}->{"$pos[1]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[1]:$flank_end"});
			$positive_sj{$chr}->{"$flank_start:$pos[2]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[2]"});
			$positive_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type1_sj_reads = ($positive_sj{$chr}->{"$flank_start:$pos[0]"}).",".($positive_sj{$chr}->{"$pos[1]:$flank_end"});
			my $type2_sj_reads = ($positive_sj{$chr}->{"$flank_start:$pos[2]"}).",".($positive_sj{$chr}->{"$pos[3]:$flank_end"});
			if($line[8]=~/3\-4\^\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}elsif($line[8]=~/1\-2\^\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
			# print OUT join("\t",@line[0..8])," ID \"MXE_$count\"\;\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[0]"});
			$negative_sj{$chr}->{"$pos[1]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[1]:$flank_end"});
			$negative_sj{$chr}->{"$flank_start:$pos[2]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[2]"});
			$negative_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type1_sj_reads = ($negative_sj{$chr}->{"$flank_start:$pos[0]"}).",".($negative_sj{$chr}->{"$pos[1]:$flank_end"});
			my $type2_sj_reads = ($negative_sj{$chr}->{"$flank_start:$pos[2]"}).",".($negative_sj{$chr}->{"$pos[3]:$flank_end"});
			if($line[8]=~/3\-4\^\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}elsif($line[8]=~/1\-2\^\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_AFE_gff($gff_file)
sub load_AFE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"null\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_end = $1;
			my $splice_position = $2;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[3]:$flank_end"});
			$positive_sj{$chr}->{"$pos[1]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[1]:$flank_end"});
			my $type1_sj_reads = $positive_sj{$chr}->{"$pos[1]:$flank_end"};
			my $type2_sj_reads = $positive_sj{$chr}->{"$pos[3]:$flank_end"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"null\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_end = $1;
			my $splice_position = $2;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[3]:$flank_end"});
			$negative_sj{$chr}->{"$pos[1]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[1]:$flank_end"});
			my $type1_sj_reads = $negative_sj{$chr}->{"$pos[1]:$flank_end"};
			my $type2_sj_reads = $negative_sj{$chr}->{"$pos[3]:$flank_end"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_ALE_gff($gff_file)
sub load_ALE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,null\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $splice_position = $2;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[0]"});
			$positive_sj{$chr}->{"$flank_start:$pos[2]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[2]"});
			my $type1_sj_reads = $positive_sj{$chr}->{"$flank_start:$pos[2]"};
			my $type2_sj_reads = $positive_sj{$chr}->{"$flank_start:$pos[0]"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,null\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $splice_position = $2;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[0]"});
			$negative_sj{$chr}->{"$flank_start:$pos[2]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[2]"});
			my $type1_sj_reads = $negative_sj{$chr}->{"$flank_start:$pos[2]"};
			my $type2_sj_reads = $negative_sj{$chr}->{"$flank_start:$pos[0]"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_A3SS_gff($gff_file)
sub load_A3SS_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,null\-\"\;.*splice_chain\s\"(\d+)\-\,(\d+)\-\S*\"\;/;
			my $flank_start = $1;
			my $alt_position = $2;
			my $model_position = $3;
			$positive_sj{$chr}->{"$flank_start:$model_position"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$model_position"});
			$positive_sj{$chr}->{"$flank_start:$alt_position"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$alt_position"});
			my $type2_sj_reads = $positive_sj{$chr}->{"$flank_start:$model_position"};
			my $type1_sj_reads = $positive_sj{$chr}->{"$flank_start:$alt_position"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,null\-\"\;.*splice_chain\s\"(\d+)\-\,(\d+)\-\S*\"\;/;
			my $flank_start = $1;
			my $alt_position = $2;
			my $model_position = $3;
			$negative_sj{$chr}->{"$flank_start:$model_position"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$model_position"});
			$negative_sj{$chr}->{"$flank_start:$alt_position"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$alt_position"});
			my $type2_sj_reads = $negative_sj{$chr}->{"$flank_start:$model_position"};
			my $type1_sj_reads = $negative_sj{$chr}->{"$flank_start:$alt_position"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_A5SS_gff($gff_file)
sub load_A5SS_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"null\^\,(\d+)\-\"\;.*splice_chain\s\"(\d+)\^\,(\d+)\^\"\;/;
			my $flank_end = $1;
			my $alt_position = $2;
			my $model_position = $3;
			$positive_sj{$chr}->{"$model_position:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$model_position:$flank_end"});
			$positive_sj{$chr}->{"$alt_position:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$alt_position:$flank_end"});
			my $type2_sj_reads = $positive_sj{$chr}->{"$model_position:$flank_end"};
			my $type1_sj_reads = $positive_sj{$chr}->{"$alt_position:$flank_end"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"null\^\,(\d+)\-\"\;.*splice_chain\s\"(\d+)\^,(\d+)\^\"\;/;
			my $flank_end = $1;
			my $alt_position = $2;
			my $model_position = $3;
			$negative_sj{$chr}->{"$model_position:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$model_position:$flank_end"});
			$negative_sj{$chr}->{"$alt_position:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$alt_position:$flank_end"});
			my $type2_sj_reads = $negative_sj{$chr}->{"$model_position:$flank_end"};
			my $type1_sj_reads = $negative_sj{$chr}->{"$alt_position:$flank_end"};
			$line[5]=$type1_sj_reads;
			$line[7]=$type2_sj_reads;
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_II_gff($gff_file)
sub load_II_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/splice_position\s\"(\S+)\"\;/;
			my $splice_position = $1;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$pos[0]:$pos[3]"} = 0 if not defined($positive_sj{$chr}->{"$pos[0]:$pos[3]"});
			$positive_sj{$chr}->{"$pos[1]:$pos[2]"} = 0 if not defined($positive_sj{$chr}->{"$pos[1]:$pos[2]"});
			my $type1_sj_reads = $positive_sj{$chr}->{"$pos[0]:$pos[3]"};
			my $type2_sj_reads = $positive_sj{$chr}->{"$pos[1]:$pos[2]"};
			if($line[8]=~/1\^4\-\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}elsif($line[8]=~/2\^3\-\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/splice_position\s\"(\S+)\"\;/;
			my $splice_position = $1;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$pos[0]:$pos[3]"} = 0 if not defined($negative_sj{$chr}->{"$pos[0]:$pos[3]"});
			$negative_sj{$chr}->{"$pos[1]:$pos[2]"} = 0 if not defined($negative_sj{$chr}->{"$pos[1]:$pos[2]"});
			my $type1_sj_reads = $negative_sj{$chr}->{"$pos[0]:$pos[3]"};
			my $type2_sj_reads = $negative_sj{$chr}->{"$pos[1]:$pos[2]"};
			if($line[8]=~/1\^4\-\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}elsif($line[8]=~/2\^3\-\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_EE_gff($gff_file)
sub load_EE_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[0]"});
			$positive_sj{$chr}->{"$pos[2]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[2]:$flank_end"});
			$positive_sj{$chr}->{"$flank_start:$pos[1]"} = 0 if not defined($positive_sj{$chr}->{"$flank_start:$pos[1]"});
			$positive_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($positive_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type1_sj_reads = ($positive_sj{$chr}->{"$flank_start:$pos[0]"}).",".($positive_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type2_sj_reads = ($positive_sj{$chr}->{"$flank_start:$pos[1]"}).",".($positive_sj{$chr}->{"$pos[2]:$flank_end"});
			if($line[8]=~/1\-4\^\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}elsif($line[8]=~/2\-3\^\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/flanks\s\"(\d+)\^\,(\d+)\-\"\;.*splice_position\s\"(\S+)\"\;/;
			my $flank_start = $1;
			my $flank_end = $2;
			my $splice_position = $3;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$flank_start:$pos[0]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[0]"});
			$negative_sj{$chr}->{"$pos[2]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[2]:$flank_end"});
			$negative_sj{$chr}->{"$flank_start:$pos[1]"} = 0 if not defined($negative_sj{$chr}->{"$flank_start:$pos[1]"});
			$negative_sj{$chr}->{"$pos[3]:$flank_end"} = 0 if not defined($negative_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type1_sj_reads = ($negative_sj{$chr}->{"$flank_start:$pos[0]"}).",".($negative_sj{$chr}->{"$pos[3]:$flank_end"});
			my $type2_sj_reads = ($negative_sj{$chr}->{"$flank_start:$pos[1]"}).",".($negative_sj{$chr}->{"$pos[2]:$flank_end"});
			if($line[8]=~/1\-4\^\,/){
				$line[5]=$type1_sj_reads;
				$line[7]=$type2_sj_reads;
			}elsif($line[8]=~/2\-3\^\,/){
				$line[5]=$type2_sj_reads;
				$line[7]=$type1_sj_reads;
			}
			my $alt_sj_avg=&avgsj($line[5]);
			my $model_sj_avg=&avgsj($line[7]);
			next if ($alt_sj_avg+$model_sj_avg==0);
			my $alt_ratio = sprintf("%.2f",$alt_sj_avg/($alt_sj_avg+$model_sj_avg));
			# print OUT join("\t",@line[0..8])," ID \"SE_$count\"\;\n";
			print OUT join("\t",@line[0..8]),"\t",$alt_ratio,"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

#子过程-----load_EE_gff($gff_file)
sub load_IR_gff {
	my ($gff_file) = @_;
	my $count = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		my $chr = $line[0];
		$count++;
		if($line[6] eq "+"){
			$line[8]=~/splice_position\s\"(\S+)\"\;/;
			my $splice_position = $1;
			my @pos = split(",",$splice_position);
			$positive_sj{$chr}->{"$pos[0]:$pos[1]"} = 0 if not defined($positive_sj{$chr}->{"$pos[0]:$pos[1]"});
			my $sj_reads = $positive_sj{$chr}->{"$pos[0]:$pos[1]"};
			$line[7]=$sj_reads;
			print OUT join("\t",@line[0..8]),"\n";
		}
		if($line[6] eq "-"){
			$line[8]=~/splice_position\s\"(\S+)\"\;/;
			my $splice_position = $1;
			my @pos = split(",",$splice_position);
			$negative_sj{$chr}->{"$pos[0]:$pos[1]"} = 0 if not defined($negative_sj{$chr}->{"$pos[0]:$pos[1]"});
			my $sj_reads = $negative_sj{$chr}->{"$pos[0]:$pos[1]"};
			$line[7]=$sj_reads;
			print OUT join("\t",@line[0..8]),"\n";
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
}

###子过程-----load_all_sj($sj_file)
sub load_all_sj {
	my ($sj_file) = @_;
	open SJ, $sj_file || die;
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
			if(defined($positive_sj{$chr}->{"$start:$end"})){
				$positive_sj{$chr}->{"$start:$end"} += $readsNum;
			}else{
				$positive_sj{$chr}->{"$start:$end"} = $readsNum;
			}
		}
		else {
			$start = $right - $rightsize + 1;
			$end   = $left + $leftsize;
			if(defined($negative_sj{$chr}->{"$start:$end"})){
				$negative_sj{$chr}->{"$start:$end"} += $readsNum;
			}else{
				$negative_sj{$chr}->{"$start:$end"} = $readsNum;
			}
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

sub sumsj {
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	return sum(@num);

}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	my $avg = (sum(@num)) / (scalar(@num)) ;
	return $avg;
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