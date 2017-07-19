#!/usr/bin/perl -w
# 
# Copyright (c)   AB_Life 2011
# Writer:         xuxiong <xuxiong19880610@163.com>
# Program Date:   2012.03.16
# Modifier:       xuxiong <xuxiong19880610@163.com>
# Last Modified:  2014.08.27
# 
# 2014.08.27  修复如果某个样品没有某种AS类型的事件，结果统计会错误的bug。
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my %opts;
GetOptions(\%opts,"i=s","p=s","f=s","o=s");

if(!defined($opts{i})  || !defined($opts{p}) || !defined($opts{o}) )
{
	print <<"	Usage End.";

	Description:This programme is used for summarizing the statistics result of RNA-seq or ChIP-seq pipeline

	Example:perl $0 -i /data2/human/ChIP-seq_whu/A3_D2_H3K36me3_0319 -p _0319 -o stat_out_0319

		Version: $ver

	Usage:perl $0

    -i        indir(project directory)     must be given

    -p        postfix                     must be given

    -f        as result file                default is AS_events_result

    -o        outfile                     must be given

	Usage End.

	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################
my $current_dir = `pwd`;chomp($current_dir);
#my $outdir = defined $opts{outdir} ? $opts{outdir} : "./";
#`mkdir -p $outdir` if (!-d $outdir) ;
#$outdir = "$current_dir/$outdir" if ($outdir!~/^\/|\~/) ;

my $indir=&AbsolutePath("dir",$opts{i});
my $postfix=$opts{p};
my $result_file="AS_events_result";
$result_file=$opts{f} if defined $opts{f};
my $outfile = $opts{o};

# my $clean_dir=(glob("$indir/clean.sh.*.qsub"))[0];

chdir($indir);
my @sample=();
&load_matched_file($indir,\@sample,"^$postfix\\S*");
# print "^\\S+$postfix\$";
# map {$_=~s/$postfix$//g;} @sample;

print join("\t",@sample),"\n";
my %map_stat=();
# my @qsub_stdout=();
# &load_matched_file($clean_dir,\@qsub_stdout,'^clean_\d+.sh.o\d+');
# &read_qsub_stdout(\@qsub_stdout,\%map_stat);

# my @qsub_stderr=();
# &load_matched_file($clean_dir,\@qsub_stderr,'^clean_\d+.sh.e\d+');
# &read_qsub_stderr(\@qsub_stderr,\%map_stat);

my @map_result=();
my @map_region=();
print "\n$sample[0]\n";
if (-d "$sample[0]") {
	@map_result=sort glob("$indir/$postfix*/AS/$result_file");
	# @map_region=sort glob("$indir/*$postfix/*.region");
}
else{
	&load_matched_file($indir,\@map_result,"^$postfix\\S*$result_file\$");
}

# print join("\t",@map_result),"\n";
&read_map_result(\@map_result,\%map_stat);

my %map_region_stat=();
# &read_mapped_region(\@map_region,\%map_region_stat) if (@map_region);

chdir($current_dir);
open (OUT,">$outfile") or die $!;
# print OUT join("\t",qw(Sample raw_data hight_quality adapter_dimer clean_reads after_clip_adapter too_short)),"\n";
# for (my $i=0;$i<@{$map_stat{"raw_data"}};$i++) {
# 	$map_stat{"raw_data"}->[$i]=0 if (not defined $map_stat{"raw_data"}->[$i]) ;
# 	$map_stat{"hight_quality"}->[$i]=0 if (not defined $map_stat{"hight_quality"}->[$i]) ;
# 	$map_stat{"adapter_dimer"}->[$i]=0 if (not defined $map_stat{"adapter_dimer"}->[$i]) ;
# 	$map_stat{"clean_reads"}->[$i]=0 if (not defined $map_stat{"clean_reads"}->[$i]) ;
# 	$map_stat{"after_clip_adapter"}->[$i]=0 if (not defined $map_stat{"after_clip_adapter"}->[$i]) ;
# 	$map_stat{"too_short"}->[$i]=0 if (not defined $map_stat{"too_short"}->[$i]) ;
# 	print OUT $sample[$i],"\t",$map_stat{"raw_data"}->[$i],"\t",$map_stat{"hight_quality"}->[$i],"\t",
# 		$map_stat{"adapter_dimer"}->[$i],"\t",$map_stat{"clean_reads"}->[$i],"\t",
# 		$map_stat{"after_clip_adapter"}->[$i],"\t",$map_stat{"too_short"}->[$i],"\n",
# }
print OUT join("\t",qw(Sample 3pMXE 5pMXE A3SS A3SS&ES A5SS A5SS&ES cassetteExon ES IntronR MXE)),"\n";
# print "haha",$map_stat{"clean_reads"}[0];
for (my $i=0;$i<@sample;$i++) {
	$map_stat{"3pMXE"}->[$i]=0 if (not defined $map_stat{"3pMXE"}->[$i]) ;
	$map_stat{"5pMXE"}->[$i]=0 if (not defined $map_stat{"5pMXE"}->[$i]) ;
	$map_stat{"A3SS"}->[$i]=0 if (not defined $map_stat{"A3SS"}->[$i]) ;
	$map_stat{"A3SS_ES"}->[$i]=0 if (not defined $map_stat{"A3SS_ES"}->[$i]) ;
	$map_stat{"A5SS"}->[$i]=0 if (not defined $map_stat{"A5SS"}->[$i]) ;
	$map_stat{"A5SS_ES"}->[$i]=0 if (not defined $map_stat{"A5SS_ES"}->[$i]) ;
	$map_stat{"cassetteExon"}->[$i]=0 if (not defined $map_stat{"cassetteExon"}->[$i]) ;
	$map_stat{"ES"}->[$i]=0 if (not defined $map_stat{"ES"}->[$i]) ;
	$map_stat{"IntronR"}->[$i]=0 if (not defined $map_stat{"IntronR"}->[$i]) ;
	$map_stat{"MXE"}->[$i]=0 if (not defined $map_stat{"MXE"}->[$i]) ;
	print OUT $sample[$i],"\t",$map_stat{"3pMXE"}->[$i],"\t",$map_stat{"5pMXE"}->[$i],"\t",
		$map_stat{"A3SS"}->[$i],"\t",$map_stat{"A3SS_ES"}->[$i],"\t",
		$map_stat{"A5SS"}->[$i],"\t",$map_stat{"A5SS_ES"}->[$i],"\t",
		$map_stat{"cassetteExon"}->[$i],"\t",$map_stat{"ES"}->[$i],"\t",
		$map_stat{"IntronR"}->[$i],"\t",$map_stat{"MXE"}->[$i],"\n";
}

if (%map_region_stat) {
	print OUT "\n",join("\t",qw(Sample 5'UTR 3'UTR CDS Nc_exon Introns Intergenic)),"\n";
	for (my $i=0;$i<@{$map_region_stat{"intergenic"}};$i++) {
		$map_region_stat{"five_prime_UTR"}->[$i]=0 if (not defined $map_region_stat{"five_prime_UTR"}->[$i]) ;
		$map_region_stat{"three_prime_UTR"}->[$i]=0 if (not defined $map_region_stat{"three_prime_UTR"}->[$i]) ;
		$map_region_stat{"CDS"}->[$i]=0 if (not defined $map_region_stat{"CDS"}->[$i]) ;
		$map_region_stat{"noncoding_exon"}->[$i]=0 if (not defined $map_region_stat{"noncoding_exon"}->[$i]) ;
		$map_region_stat{"intron"}->[$i]=0 if (not defined $map_region_stat{"intron"}->[$i]) ;
		$map_region_stat{"intergenic"}->[$i]=0 if (not defined $map_region_stat{"intergenic"}->[$i]) ;
		print OUT $sample[$i],"\t",$map_region_stat{"five_prime_UTR"}->[$i],"\t",
			$map_region_stat{"three_prime_UTR"}->[$i],"\t",$map_region_stat{"CDS"}->[$i],"\t",
			$map_region_stat{"noncoding_exon"}->[$i],"\t",$map_region_stat{"intron"}->[$i],"\t",
			$map_region_stat{"intergenic"}->[$i],"\n";
	}
}
close OUT;
sub load_matched_file{
	my ($INDIR,$filenames_ref,$match)=@_;
	opendir(DIR,$INDIR) or die $!;
	my $tmp;
	while ($tmp=readdir(DIR)) {
		chomp $tmp;
		# print "$INDIR/$tmp/AS/$result_file\n";
		push @{$filenames_ref},$tmp if ($tmp=~/$match/) && ((-f "$INDIR/$tmp/AS/$result_file") || (-f "$INDIR/$postfix\S+$result_file"));
	}
	@{$filenames_ref} = sort @{$filenames_ref};
#	print join ("\t",@{$filenames_ref}),"\n";
	close(DIR);
}

# sub read_qsub_stdout{
# 	my ($file_list,$hash)=@_;
# 	chdir($clean_dir);
# 	foreach my $filename (@{$file_list}) {
# 		open(IN,$filename) or die $!;
# 		while (<IN>) {
# 			chomp;
# 			if (/^raw:\s*(\d+)$/) {
# 				push @{$hash->{"raw_data"}},$1;
# 			}
# 			elsif(/^no_NN\s*(\d+)$/){
# 				push @{$hash->{"no_NN"}},$1;
# 			}
# 		}
# 		close IN;
# 	}
# 	chdir($indir);
# }

# sub read_qsub_stderr{
# 	my ($file_list,$hash)=@_;
# 	chdir($clean_dir);
# 	foreach my $filename (@{$file_list}) {
# 		open(IN,$filename) or die $!;
# 		while (<IN>) {
# 			chomp;
# 			if (/^Input:\s*(\d+)\s*reads\.$/) {
# 				push @{$hash->{"hight_quality"}},$1;
# 			}
# 			elsif(/^Output:\s*(\d+)\s*reads\.$/){
# 				push @{$hash->{"after_clip_adapter"}},$1;
# 			}
# 			elsif(/^discarded\s*(\d+)\s*too-short\s*reads\.$/){
# 				push @{$hash->{"too_short"}},$1;
# 			}
# 			elsif(/^discarded\s*(\d+)\s*adapter-only\s*reads\.$/){
# 				push @{$hash->{"adapter_dimer"}},$1;
# 			}
# 		}
# 		close IN;
# 	}
# 	chdir($indir);
# }


sub read_map_result{
	my ($file_list,$hash) = @_;
	foreach my $filename (@{$file_list}) {
		print $filename,"\n";
		open(IN,$filename) or die $!;
		my %tmphash = ();
		while (<IN>) {
			chomp;
			print $_,"\n";
			if (/(\d+)\s3pMXE$/) {
				push @{$hash->{"3pMXE"}},$1;
				$tmphash{"3pMXE"}=1;
			}
			elsif(/(\d+)\s5pMXE$/){
				push @{$hash->{"5pMXE"}},$1;
				$tmphash{"5pMXE"}=1;
			}
			elsif(/(\d+)\sA3SS$/){
				push @{$hash->{"A3SS"}},$1;
				$tmphash{"A3SS"}=1;
			}
			elsif(/(\d+)\sA3SS\&ES$/){
				push @{$hash->{"A3SS_ES"}},$1;
				$tmphash{"A3SS_ES"}=1;
			}
			elsif(/(\d+)\sA5SS$/){
				push @{$hash->{"A5SS"}},$1;
				$tmphash{"A5SS"}=1;
			}
			elsif(/(\d+)\sA5SS\&ES$/){
				push @{$hash->{"A5SS_ES"}},$1;
				$tmphash{"A5SS_ES"}=1;
			}
			elsif(/(\d+)\scassetteExon$/){
				push @{$hash->{"cassetteExon"}},$1;
				$tmphash{"cassetteExon"}=1;
			}
			elsif(/(\d+)\sES$/){
				push @{$hash->{"ES"}},$1;
				$tmphash{"ES"}=1;
			}
			elsif(/(\d+)\sIntronR$/){
				push @{$hash->{"IntronR"}},$1;
				$tmphash{"IntronR"}=1;
			}
			elsif(/(\d+)\sMXE$/){
				push @{$hash->{"MXE"}},$1;
				$tmphash{"MXE"}=1;
			}
		}
		if(not defined($tmphash{"MXE"})){
			push @{$hash->{"MXE"}},0;
		}
		if(not defined($tmphash{"IntronR"})){
			push @{$hash->{"IntronR"}},0;
		}
		if(not defined($tmphash{"ES"})){
			push @{$hash->{"ES"}},0;
		}
		if(not defined($tmphash{"cassetteExon"})){
			push @{$hash->{"cassetteExon"}},0;
		}
		if(not defined($tmphash{"A5SS_ES"})){
			push @{$hash->{"A5SS_ES"}},0;
		}
		if(not defined($tmphash{"A5SS"})){
			push @{$hash->{"A5SS"}},0;
		}
		if(not defined($tmphash{"A3SS_ES"})){
			push @{$hash->{"A3SS_ES"}},0;
		}
		if(not defined($tmphash{"A3SS"})){
			push @{$hash->{"A3SS"}},0;
		}
		if(not defined($tmphash{"5pMXE"})){
			push @{$hash->{"5pMXE"}},0;
		}
		if(not defined($tmphash{"3pMXE"})){
			push @{$hash->{"3pMXE"}},0;
		}
		close IN;
	}
}

sub read_mapped_region{
	my ($file_list,$hash)=@_;
	foreach my $filename (@{$file_list}) {
		my $temp_sum=0;
		open(IN,$filename) or die $!;
		while (<IN>) {
			chomp;
			$temp_sum+=(split(/\s+/,$_))[1];
		}
		close IN;

		open(IN,$filename) or die $!;
		while (<IN>) {
			chomp;
			my @line=split(/\s+/,$_);
			push @{$hash->{$line[0]}},$line[1]."(".sprintf("%3.2f%%",$line[1]/$temp_sum*100).")";
		}
		close IN;
	}
#	print Dumper %{$hash};
}

sub AbsolutePath{		#获取指定目录或文件的决定路径
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
###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
