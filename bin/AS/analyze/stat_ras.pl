#!/usr/bin/perl -w
# 
# Copyright (c)   AB_Life 2011
# Writer:         xuxiong <xuxiong19880610@163.com>
# Program Date:   2012.03.16
# Modifier:       xuxiong <xuxiong19880610@163.com>
# Last Modified:  2011.03.16
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my %opts;
GetOptions(\%opts,"i=s","p=s","m=s","o=s");

if(!defined($opts{i}) || !defined($opts{p}) || !defined($opts{o}) )
{
	print <<"	Usage End.";

	Description:This programme is used for summarizing the statistics result of RNA-seq or ChIP-seq pipeline

	Example:perl $0 -i /data2/human/ChIP-seq_whu/A3_D2_H3K36me3_0319 -p _0319 -o stat_out_0319

		Version: $ver

	Usage:perl $0

    -i        indir(project directory)     must be given

    -p        postfix                     must be given

    -m        match                       default is "_map_result"

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
my $match= "_result";
$match=$opts{m} if defined($opts{m});
my $outfile = $opts{o};

chdir($indir);
my @sample=();
&load_matched_file($indir,\@sample,"^$postfix\\S+$match\$");
# print "^\\S+$postfix\$";
# map {$_=~s/$postfix$//g;} @sample;

# print join("\t",@sample),"\n";
my %map_stat=();

my @map_result=();
# print "\n$sample[0]$postfix\n";
if (-d "$sample[0]$postfix") {
	@map_result=sort glob("$indir/*$postfix/*$match");
	# print "yes","\n";
}
else{
	&load_matched_file($indir,\@map_result,"^$postfix\\S+$match\$");
}
# print join("\t",@map_result),"\n";
&read_map_result(\@map_result,\%map_stat);

chdir($current_dir);

open (OUT,">$outfile") or die $!;


for (my $i=0;$i<@sample;$i++) {
	if($sample[$i]=~/_IR_fisher_test/){
		my $tmp = $map_stat{"IntronR"}{$sample[$i]};
		my $tmpsample=$sample[$i];
		$tmpsample=~s/_IR_fisher_test/_fisher_test/;
		$map_stat{"IntronR"}{$tmpsample}=$tmp;
	}
}

print OUT "Sample";
foreach my $key (sort keys %map_stat){
	next if $key=~/AS_TYPE/;
	print OUT "\t$key";
}
print OUT "\n";

# print "haha",$map_stat{"clean_reads"}[0];
my $thissample = "$postfix\_fisher_test_result";
print OUT "Total AS Number";
foreach my $key (sort keys %map_stat){
	next if $key=~/AS_TYPE/;
	$map_stat{$key}->{$thissample}=0 if not defined($map_stat{$key}->{$thissample});
	print OUT "\t",$map_stat{$key}->{$thissample};
}
print OUT "\n";

$thissample = "$postfix\_fisher_test_diff_up_result";
print OUT "Up";
foreach my $key (sort keys %map_stat){
	next if $key=~/AS_TYPE/;
	$map_stat{$key}->{$thissample}=0 if not defined($map_stat{$key}->{$thissample});
	print OUT "\t",$map_stat{$key}->{$thissample};
}
print OUT "\n";

$thissample = "$postfix\_fisher_test_diff_down_result";
print OUT "Down";
foreach my $key (sort keys %map_stat){
	next if $key=~/AS_TYPE/;
	$map_stat{$key}->{$thissample}=0 if not defined($map_stat{$key}->{$thissample});
	print OUT "\t",$map_stat{$key}->{$thissample};
}
print OUT "\n";

$thissample = "$postfix\_fisher_test_diff_result";
print OUT "Total changed AS Number";
foreach my $key (sort keys %map_stat){
	next if $key=~/AS_TYPE/;
	$map_stat{$key}->{$thissample}=0 if not defined($map_stat{$key}->{$thissample});
	print OUT "\t",$map_stat{$key}->{$thissample};
}
print OUT "\n";

close OUT;

&tras("$indir/$outfile");

sub load_matched_file{
	my ($INDIR,$filenames_ref,$match)=@_;
	# print $INDIR;
	opendir(DIR,$INDIR) or die $!;
	my $tmp;
	while ($tmp=readdir(DIR)) {
		chomp $tmp;
		push @{$filenames_ref},$tmp if ($tmp=~/$match/) ;
	}
	@{$filenames_ref} = sort @{$filenames_ref};
#	print join ("\t",@{$filenames_ref}),"\n";
	close(DIR);
}


sub read_map_result{
	my ($file_list,$hash) = @_;
	my $i=0;
	foreach my $filename (@{$file_list}) {
		open(IN,$filename) or die $!;
		while (<IN>) {
			chomp;
			next if(/^Start\s+Time/);
			next if(/^End\s+Time/);
			$_=~s/^\s+//;
			if(/\s+/){
				my @line=split(/\s+/,$_,2);
				$hash->{$line[1]}->{$sample[$i]}=$line[0];
				# print $j,"\n";
			}
		}
		close IN;
		$i++;
	}
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

sub tras{
	my ($file)=@_;
	print $file,"\n";
	my @matrix;
	my $max_len = 0;
	open(IN ,"$file") or die $!;
	while(<IN>){
	  chomp;
	  print $_,"\n";
	  s/\t+$//g;
	  my @fields = split /\t/, $_;
	  my $len = @fields;
	  $max_len = $max_len > $len ? $max_len : $len;
	  push @matrix, [@fields];
	}
	close IN;
	open IN ,">$file"  or die $!;;
	for my $col (0..$max_len - 1){
	  for my $line (@matrix){
	    print IN $line->[$col] || ' ', "\t";
	  }
	  print IN "\n";
	}
	close IN;
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
