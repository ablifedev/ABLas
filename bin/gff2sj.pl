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
GetOptions( \%opts, "gff=s", "o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff        annotation file        must be given;

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
my $out_file = $opts{o};
open OUT, ">$out_file" || die;


###读取gff的内容，记录每个转录本的splice site。
my %gff = ();
&load_gff(\%gff,$gff_file);


my %sj=();
###以每个基因作为model找可变剪接事件。
foreach my $chr (sort keys %gff){
	foreach my $strand (sort keys %{$gff{$chr}}){
		foreach my $gene (sort keys %{$gff{$chr}->{$strand}}){
			my @thisgene = @{$gff{$chr}->{$strand}->{$gene}};
			# print $gene,"\n";
			&gene2sj($chr,$strand,\@thisgene);   ##find Mutually exclusive exons events
		}
	}
}

close OUT;



###子过程-----gene2sj($chr,$strand,$gene,\@thisgene)
sub gene2sj {
	my ($chr,$strand,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			
			my $sj_start = $genesite[2*$f+1];
			my $sj_end   = $genesite[2*$f+2];
			# print "$strand\t$sj_start\t$sj_end\n";
			my $junc = "$sj_start\t$sj_end";
			if(defined($sj{$chr}->{$strand}->{$junc})){
				next;
			}else{
				$sj{$chr}->{$strand}->{$junc} = 1;
				print OUT "$chr\t$junc\t$strand\n";
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			
			my $sj_start = $genesite[2*$f+1];
			my $sj_end   = $genesite[2*$f+2];
			# print "$strand\t$sj_start\t$sj_end\n";
			my $junc = "$sj_start\t$sj_end";
			if(defined($sj{$chr}->{$strand}->{$junc})){
				next;
			}else{
				$sj{$chr}->{$strand}->{$junc} = 1;
				print OUT "$chr\t$junc\t$strand\n";
			}
		}
	}
}

#子过程-----load_gff(\%gff,$gff_file)
sub load_gff {
	my ($gff, $gff_file) = @_;
	my $gene = "";
	my $this_chr = "";
	my $this_strand = "";
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my (
		$chr, $source, $feature, $start, $end,
		$other1, $strand, $other2, $info
		)= split /\t/;
		next if ( $feature =~ /chromosome/  );
		if($feature =~ /^gene$/){
			$gene = "";
			next;
		}
		if ( $feature =~ /^(mRNA|miRNA|mRNA_TE_gene|ncRNA|rRNA|snoRNA|snRNA|tRNA|transcript)$/ ) {
			#将上一个gene的各位点排序，正向的按从小到大，负向的按从大到小。
			if($gene ne "" ){
				if($this_strand eq "+"){
					if(defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
						@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$a<=>$b} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
						}
				}else{
					if(defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
						@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$b<=>$a} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
						}
				}
			}
			#开始记录新的基因，以ID为转录本（gene）名称
			$this_chr = $chr;
			$this_strand = $strand;
			$info =~ m/^ID=([\w|\.|\%|\:|\_|\-|\s]+);/;
			$gene = $1;
		}
		if ( $feature =~ /^exon$/ ){
			#将exon的起始和终止位置放进gene数组。
			next if $gene eq "";
			push @{$gff->{$this_chr}->{$this_strand}->{$gene}},$start;
			push @{$gff->{$this_chr}->{$this_strand}->{$gene}},$end;
		}
	}
	if($gene ne ""){
		if($this_strand eq "+"){
			if(defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
				@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$a<=>$b} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
				}
		}else{
			if(defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
				@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$b<=>$a} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
				}
		}
	}
	close GFF;
	#调试信息
	print "done reading GFF ...........", "\n\n";
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