#!/usr/bin/perl -w
#
# Copyright (c)   A_B_Life 2013
# Writer:         ChengChao <chaocheng@ablife.cc>
# Program Date:   2013.8.21
# Modifier:       
# Last Modified:  
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your program you must write the detailed time discriptions parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions( \%opts, "config=s", "h" );
if ( !defined( $opts{config} ) ) {
	print <<"	Usage End.";
	Description:This program is used for do NAS and RAS for each sample

		Version: $ver

	Usage:perl $0

		-config               config file                               must be given

	Usage End.
	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $config = $opts{config};
my %ConfigOption = ();
&Get_config_option($config,\%ConfigOption);

my $indir = $ConfigOption{"INDIR"};
my $species = $ConfigOption{"SPECIES"};
my $gff = $ConfigOption{"GFF"};
my $fasta = $ConfigOption{"FASTA"};
my $RAS = $ConfigOption{"RAS"};
my $outdir = $ConfigOption{"OUTDIR"};

my @samplename = ();
my @samplebed  = ();
my @samplesam  = ();
my @samplebam  = ();

my $current_dir = `pwd`;
chomp($current_dir);
`mkdir -p $outdir` if ( !-d $outdir );
$outdir = "$current_dir/$outdir" if ( $outdir !~ /^\/|\~/ );
my $tmpdir = "";
my $SH = "";


##对每个sample do NAS

chdir($outdir);
$SH="$outdir/NAS.sh";
# print $SH;
open(SH ,">$SH") || die $!;
foreach my $sample (sort keys %{$ConfigOption{"SAMPLE"}}) {
	my $bedfile = $ConfigOption{"SAMPLE"}->{$sample}->[0];
	push @samplebed,"$indir/$bedfile";
	my $samfile = $ConfigOption{"SAMPLE"}->{$sample}->[1];
	push @samplesam,"$indir/$samfile";
	my $bamfile = $ConfigOption{"SAMPLE"}->{$sample}->[2];
	push @samplebam,"$indir/$bamfile";
	my $sample_name = $ConfigOption{"SAMPLE"}->{$sample}->[3];
	push @samplename,$sample_name;

	##判断文件是否存在
	open(IN,"<$indir/$bedfile") or die "can't open $bedfile!\n";
	close IN;
	open(IN,"<$indir/$samfile") or die "can't open $samfile!\n";
	close IN;

	print SH "cd $outdir && python2.7 $Bin/pipeline.py -gff $gff -allsj $indir/$bedfile -sam $indir/$samfile -bam $indir/$bamfile -fa $fasta -outDir $outdir/$sample_name && \n" ;
}
close SH;

###detect genome known AS events
chdir($outdir);
$tmpdir = $outdir."/genome_known_as";
`mkdir -p $tmpdir` if ( !-d $tmpdir );
$SH="$outdir/detect_known_as.sh";
open(SH ,">$SH") || die $!;
my $bed_str = join(" ",@samplebed);
print SH "cd $outdir/genome_known_as && perl $Bin/bin/gff2sj.pl -gff $gff -o gffsj && perl $Bin/bin/ABlas.pl -gff $gff -sj gffsj -allsj gffsj -g 1 -o genome_known_as &&\n";
print SH "cat $bed_str > all_sj.bed && sh $Bin/bin/known_as_uniq.sh && sh $Bin/bin/known_as_sj.sh $Bin/bin/ABlas_quantify.pl all_sj.bed&& sh $Bin/bin/known_as_filter.sh";
close SH;
# `sh $SH`;

`sh detect_known_as.sh`;
`perl /public/bin/qsub-sge.pl --queue all.q --resource vf=15.0G --maxproc 30 NAS.sh ` ;

##END 对每个sample do NAS



# ####统计
# chdir($outdir);
# $tmpdir = $outdir."/Stat";
# `mkdir -p $tmpdir` if ( !-d $tmpdir );
# $tmpdir = $outdir."/Stat/EXP";
# `mkdir -p $tmpdir` if ( !-d $tmpdir );
# $tmpdir = $outdir."/Stat/SJ";
# `mkdir -p $tmpdir` if ( !-d $tmpdir );
# $tmpdir = $outdir."/Stat/AS";
# `mkdir -p $tmpdir` if ( !-d $tmpdir );
# $SH="$outdir/stat.sh";
# open(SH ,">$SH") || die $!;


# ##EXP
# print SH "cd $outdir/Stat/EXP && ";
# foreach my $name (@samplename){
# 	print SH "more $outdir/$name/Expression/EXP_exon | awk '\$7>0' | wc -l > $name"."_EXP_exon_num && more $outdir/$name/Expression/EXP_exon | awk '{print \$7}' > $name"."_EXP_exon_base &&" ;
# }

# print SH "paste ".join("\ ",(map {$_."_EXP_exon_base"} @samplename))." > ".join("_",@samplename)."_EXP_exon && \n";

# ##SJ
# print SH "cd $outdir/Stat/SJ && ";
# print SH "cat ".join("\ ",@samplebed)." > raw_all_sj.bed && perl $Bin/bin/SJ/get_totaluniqsj.pl -i raw_all_sj.bed -o all_sj.bed && perl $Bin/bin/SJ/get_known_novel_sj.pl -asj $outdir/$samplename[0]/temp/gffSj.bed -allsj all_sj.bed -o all_sj_novel.bed -o2 all_sj_known.bed && more all_sj.bed | wc -l > all_sj_num && more all_sj_novel.bed | wc -l > all_sj_novel_num && more all_sj_known.bed | wc -l > all_sj_known_num && \n";

# ##AS
# print SH "cd $outdir/Stat/AS && ";
# print SH "cat ".join("\ ",(map {"$outdir/".$_."/AS/AS_events"} @samplename))." > AS_events_cat && perl $Bin/bin/AS/analyze/get_totalas.pl -i AS_events_cat -o total_AS_events && ";
# print SH "cat ".join("\ ",(map {"$outdir/".$_."/AS/AS_events_known"} @samplename))." > AS_events_known_cat && perl $Bin/bin/AS/analyze/get_totalas.pl -i AS_events_known_cat -o total_AS_events_known && ";
# print SH "cat ".join("\ ",(map {"$outdir/".$_."/AS/AS_events_novel"} @samplename))." > AS_events_novel_cat && perl $Bin/bin/AS/analyze/get_totalas.pl -i AS_events_novel_cat -o total_AS_events_novel && \n";

# close SH;
# `perl /public/bin/qsub-sge.pl --queue all.q --resource vf=10.0G --maxproc 30 $SH ` ;

# ####END 统计


# #### DO RAS


if ($RAS == 1) {
	chdir($outdir);
	$tmpdir = $outdir."/RAS";
	`mkdir -p $tmpdir` if ( !-d $tmpdir );

	$SH="$outdir/ras.sh";
	open(SH ,">$SH") || die $!;

	foreach my $group (sort keys %{$ConfigOption{"GROUP"}}) {
		my $sample1 = $ConfigOption{"GROUP"}->{$group}->[0];
		my $sample2 = $ConfigOption{"GROUP"}->{$group}->[1];
		my $sample1_bed = "$indir/".$ConfigOption{'SAMPLE'}->{$sample1}->[0];
		my $sample2_bed = "$indir/".$ConfigOption{'SAMPLE'}->{$sample2}->[0];
		my $sample1_bam = "$indir/".$ConfigOption{'SAMPLE'}->{$sample1}->[2];
		my $sample2_bam = "$indir/".$ConfigOption{'SAMPLE'}->{$sample2}->[2];
		# print $sample1_bed,"\n";
		# print $sample2_bed,"\n";
		# print SH "cd $outdir/RAS && ";
		# print SH "perl $Bin/bin/AS/ras/ras_v3.pl -i1 $outdir/$sample2/AS/availabe_as_events -i2 $outdir/$sample1/AS/availabe_as_events -j1 $indir/$sample2_bed -j2 $indir/$sample1_bed -o $sample1"."_vs_"."$sample2 && R --slave </users/chengc/work/AS/CASC/bin/AS/ras/Fisher_exact_test.r --args $sample1"."_vs_"."$sample2"." $sample1"."_vs_"."$sample2"."_fisher_test ./ && less -S $sample1"."_vs_"."$sample2"."_fisher_test |perl -ne 'chomp;\@line=split;if(\$line[-3]>=0.2 || \$line[-3]<=-0.2){if(\$line[-1]<=0.05){print \$_,\"\\n\";}}' > $sample1"."_vs_"."$sample2"."_fisher_test_diff && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff | awk '{print \$13}'|perl -ne 'chomp;\$_=~s/\\.\\S+//;print \$_,\"\\n\";' |sort|uniq> $sample1"."_vs_"."$sample2"."_fisher_test_diff_gene && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff | awk '{print \$5}' | awk -F '|' '{print \$1}'|sort|uniq -c > $sample1"."_vs_"."$sample2"."_fisher_test_diff_result && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff |awk '\$12>0'|more > $sample1"."_vs_"."$sample2"."_fisher_test_diff_up && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff |awk '\$12<0'|more > $sample1"."_vs_"."$sample2"."_fisher_test_diff_down && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff_up | awk '{print \$13}'|perl -ne 'chomp;\$_=~s/\\.\\S+//;print \$_,\"\\n\";' |sort|uniq> $sample1"."_vs_"."$sample2"."_fisher_test_diff_up_gene && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff_up | awk '{print \$5}' | awk -F '|' '{print \$1}'|sort|uniq -c > $sample1"."_vs_"."$sample2"."_fisher_test_diff_up_result && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff_down | awk '{print \$13}'|perl -ne 'chomp;\$_=~s/\\.\\S+//;print \$_,\"\\n\";' |sort|uniq> $sample1"."_vs_"."$sample2"."_fisher_test_diff_down_gene && less -S $sample1"."_vs_"."$sample2"."_fisher_test_diff_down | awk '{print \$5}' | awk -F '|' '{print \$1}'|sort|uniq -c > $sample1"."_vs_"."$sample2"."_fisher_test_diff_down_result &&\n";
		$tmpdir = $outdir."/RAS/$sample2"."_vs_$sample1";
		`mkdir -p $tmpdir` if ( !-d $tmpdir );
		print SH "cd $tmpdir && ";
		print SH "sh $Bin/bin/ras.sh $outdir $sample1 $sample2 $Bin/bin/ABlas_quantify.pl $sample1_bed $sample2_bed $Bin/bin/Fisher_exact_test.r $Bin/bin/intronR_quantify.pl $sample1_bam $sample2_bam $fasta &&\n";
	}
	close SH;
	# `bash $SH`;
	`perl /public/bin/qsub-sge.pl --queue all.q --resource vf=10.0G --maxproc 30 $SH ` ;
}


sub Get_config_option {
	my ($file,$hash_ref) = @_;
	#print $file,"\n";
	open (CONFIG,"<$file") or die "Can't open $file!\n";
	while (<CONFIG>) {
		chomp;
		next if(/^\#/);
		next if(/^\s*$/);
		s/\r//g;		#chomp the \r character
		my @line = split(/\t/);
		if ($line[0] =~ /sample/i) {
			my $s_name = (split(":",$line[1]))[3];
			@{$hash_ref->{"SAMPLE"}->{$s_name}} = split(":",$line[1]);

			# @{$hash_ref->{"SAMPLE"}->{$line[0]}} = split(":",$line[1]);
		} elsif ($line[0] =~ /group/i) {
			@{$hash_ref->{"GROUP"}->{$line[0]}} = split(":",$line[1]);
		} elsif ($line[0] =~ /coras/i) {
			@{$hash_ref->{"CORAS"}->{$line[0]}} = split(":",$line[1]);
		} else {
			$hash_ref->{$line[0]} = $line[1];
		}
	}
	close(CONFIG);
}

# sub DEG_calculate {
# 	my ($option_ref, $reads_file, $RPKM_file) = @_;
# 	foreach my $group (sort keys %{$option_ref->{"GROUP"}}) {
# 		my $sample1 = $option_ref->{"GROUP"}->{$group}->[0];
# 		my $sample2 = $option_ref->{"GROUP"}->{$group}->[1];
# 		my $index1 = &select_index(\@sample_name,$sample1)+3;
# 		my $index2 = &select_index(\@sample_name,$sample2)+3;
# 		print $sample1,"\t",$sample2,"\n";
# 		`cat /users/chend/work/R/RNA-Seq/linear_models_single.r | R --slave --args $RPKM_file $log $index1 $index2 $outdir`;
# 		my $dir = $sample1."_".$sample2."_".$option_ref->{"METHOD"}."_".$option_ref->{"PVALUE"}."_".$option_ref->{"FC"};
# 		`mkdir -p $outdir/$dir`;	#make the out direcotry
# 		$dir = &AbsolutePath("dir",$dir);		#the absolute path of out directory
# 		if ($option_ref->{"METHOD"} eq "edgeR") {
# 			`cat /users/chend/work/R/RNA-Seq/edgeR_single.r | R --slave --args $reads_file $RPKM_file $option_ref->{"PVALUE"} $option_ref->{"FC"} $index1 $index2 $dir`;
# 		} elsif ($option_ref->{"METHOD"} eq "DESeq") {
			
# 		} else {
# 			die"You choose the wrong method to calculate DEG!\n";
# 		}
# 		#chdir($dir);
# 		if (exists $ConfigOption{"SYMBOL"} ) {
# 			&AddSymbolInf($dir,"\.txt",\%gene_symbol);
# 		}
# 	}
# }

# sub Co_DEG_analysis {    #&Select_co_DEG($method,$outdir,$file1,$file2,$i,$j,\%DEG_method)
# 	my ($option_ref, $outdir) = @_;
# 	`mkdir -p co_DEG` if ( !-d "co_DEG"); 
# 	foreach my $coDEG (sort keys %{$option_ref->{"CODEG"}}) {
# 		my @Sig_file = ();
# 		chdir($outdir);
# 		#print join(":",@{$option_ref->{"CODEG"}->{$coDEG}}),"\n";
# 		#print `pwd`;
# 		$outdir = &AbsolutePath("dir",$outdir);
# 		opendir(DH,$outdir) or die "$!\n";
# 		foreach my $element (readdir DH) {
# 			if (-d $element) {
# 				foreach my $group (@{$option_ref->{"CODEG"}->{$coDEG}}) {
# 					if ($element =~ /$group/i ) {
# 						opendir(DIR,$element) or die "$!\n";
# 						print $outdir."/".$element,"\n";
# 						foreach my $file (readdir DIR) {
# 							if ($file =~ /Sig_DEG.txt_symbol/i) {
# 								push @Sig_file, $outdir."/".$element."/".$file;
# 								#print $outdir."/".$element."/".$file,"\n";
# 							} else {
# 								next;
# 							}
# 						}
# 					}
# 				}
# 			} else {
# 				next;
# 			}
# 		}
# 		close(DH);
# 		print join("\n",@Sig_file),"\n";
# 		chdir("co_DEG");
# 		print `pwd`;
# 		my $co_dir = join("_and_",@{$option_ref->{"CODEG"}->{$coDEG}})."_"."co_DEG";
# 		`mkdir -p $co_dir` if ( !-d $co_dir);
# 		$co_dir = &AbsolutePath("dir",$co_dir); chdir($co_dir);
# 		print `pwd`;
# 		`perl /users/chend/work/Perl/RNA-Seq/edgeR_Select_co_DEG.pl -i_1 $Sig_file[0] -i_2 $Sig_file[1] >co_DEG_log.txt`;
# 	}
# }

#sub select_index{
#	my ($array_ref,$sample) = @_;
#	my $index = 0;
#	#print scalar(@{$array_ref}),"\t",$sample,"\n";
#	for (my $i=0; $i<scalar(@{$array_ref}); $i++) {
#		#print $array_ref->[$i],"\t",$sample,"\n";
#		if ($array_ref->[$i] eq $sample ) {		#$array_ref->[$i] =~ /$sample/i
#			$index = $i;
#		}
#	}
#	print $index,"\n";
#	return $index;
#}

sub AbsolutePath{		# Get the absolute path of the target directory or file
	my ($type,$input) = @_;
	my $return;
	if ($type eq "dir"){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	} elsif($type eq 'file'){
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
