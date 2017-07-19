#!/usr/bin/perl -w

my $ver = "1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your program you must write the detailed time discriptions parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions( \%opts, "config=s", "conir=s","coir=s", "reads=s", "o=s","h" );
if ( !defined( $opts{config} ) ) {
	print <<"	Usage End.";
	Description:This program is used for do ras between repetitions.

		Version: $ver

	Usage:perl $0

		-config               config file                               must be given
		-conir                co nir总表                                  must be given
		-coir                 co ir总表                                  must be given
		-o                    输出文件夹                                 default is ./

	Usage End.
	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################

my %nir_alt_num = (); #number of sample name
my %nir_model_num = (); #number of sample name
my %nir_ratio_num = (); #number of sample name

my %ir_alt_num = (); #number of sample name
my %ir_model_num = (); #number of sample name
my %ir_ratio_num = (); #number of sample name

my $config = $opts{config};
my $conir   = $opts{conir};
$conir = &AbsolutePath("file",$conir);
my $coir   = $opts{coir};
$coir = &AbsolutePath("file",$coir);

my %ConfigOption = ();
&Get_config_option($config,\%ConfigOption);

$ConfigOption{"pvalue"} ||= 0.05;
$ConfigOption{"diff"} ||= 0.2;

my $outdir = $opts{o} || "./";;

my $current_dir = `pwd`;
chomp($current_dir);
$outdir = "$current_dir/$outdir" if ( $outdir !~ /^\/|\~/ );
print $outdir,"\n";
`mkdir -p $outdir`;
# chdir($outdir)

my @sample_name = ();
my $sample_num = get_Sample_number( \%ConfigOption);

&RAS_calculate(\%ConfigOption,$conir,$coir);


if (exists($ConfigOption{"CORAS"})) {
	&Co_RAS_analysis(\%ConfigOption, $outdir);
}

if (exists($ConfigOption{"DIRAS"})) {
	&Di_RAS_analysis(\%ConfigOption, $outdir);
}

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

##################################

sub Get_config_option {
	my ($file,$hash_ref) = @_;
	my $start_flag = 0;
	#print $file,"\n";
	open (CONFIG,"<$file") or die "Can't open $file!\n";
	while (<CONFIG>) {
		chomp;
		next if(/^\#/);
		next if(/^\/\/\//);        # 注释行
		next if(/^\s*$/);
		s/\/\/\/.*//g;        # 行内注释
		s/\r//g;		#chomp the \r character
		if(/^\[\[ablife:(\w+)\]\]/){
			$start_flag=0;
			my $tag = $1;
			if($tag eq "config" || $tag eq "sample" || $tag eq "as_ras"){
				print $_,"\n";
				$start_flag=1;
			}
			next;
		}
		if($start_flag==0){
			next;
		}
		print $_,"\n";
		my @line = split(/\t/);
		if ($line[0] =~ /^sample/i) {
			$line[0]=~s/SAMPLE//i;
			@{$hash_ref->{"SAMPLE"}->{$line[0]}} = split(":",$line[1]);
			my $sample_name = $hash_ref->{"SAMPLE"}->{$line[0]}->[1];
			$nir_alt_num{$sample_name}=&select_index($conir,"$sample_name\:Alt_average_sj") + 1;
			$ir_alt_num{$sample_name}=&select_index($coir,"$sample_name\:intron_ave_base") + 1;
			$nir_model_num{$sample_name}=&select_index($conir,"$sample_name\:Model_average_sj") + 1;
			$ir_model_num{$sample_name}=&select_index($coir,"$sample_name\:ends_exon_ave_base") + 1;
			$nir_ratio_num{$sample_name}=&select_index($conir,"$sample_name\:Ratio") + 1;
			$ir_ratio_num{$sample_name}=&select_index($coir,"$sample_name\:Ratio") + 1;
		} elsif ($line[0] =~ /^group/i) {
			@{$hash_ref->{"GROUP"}->{$line[0]}} = split(":",$line[1]);
		} elsif ($line[0] =~ /^coRAS/i) {
			@{$hash_ref->{"CORAS"}->{$line[0]}} = split(":",$line[1]);
		} elsif ($line[0] =~ /^diRAS/i) {
			@{$hash_ref->{"DIRAS"}->{$line[0]}} = split(":",$line[1]);
		} else {
			$hash_ref->{$line[0]} = $line[1];
		}
	}
	close(CONFIG);
}

sub get_Sample_number {    #&Load_RNA_seq(\@RNA_arrReads,$outdir,$EG_files)
	my ($option_ref) = @_;

	my $sample_num = 0;

	foreach my $sample (sort keys %{$option_ref->{"SAMPLE"}}) {    #load the RNA seq of every sample
		$sample_num++;
	}

	return $sample_num;
}


sub RAS_calculate {
	my ($option_ref, $conir_file, $coir_file) = @_;

	#`R --slave </users/chend/work/R/RNA-Seq/Sample_correlation.r --args $RPKM_file $index1 $index2 $outdir`;
	foreach my $group (sort keys %{$option_ref->{"GROUP"}}) {
		chdir($outdir);
		my $sample1 = $option_ref->{"GROUP"}->{$group}->[0];
		my $sample2 = $option_ref->{"GROUP"}->{$group}->[1];
		my $nir_fishertest_list = $nir_alt_num{$sample1}.",".$nir_model_num{$sample1}.",".$nir_alt_num{$sample2}.",".$nir_model_num{$sample2};
		my $ir_fishertest_list = $ir_alt_num{$sample1}.",".$ir_model_num{$sample1}.",".$ir_alt_num{$sample2}.",".$ir_model_num{$sample2};

		my $nir_sample1_ratio = $nir_ratio_num{$sample1};
		my $nir_sample2_ratio = $nir_ratio_num{$sample2};
		my $ir_sample1_ratio = $ir_ratio_num{$sample1};
		my $ir_sample2_ratio = $ir_ratio_num{$sample2};

		print $sample1,"_vs_",$sample2,"\n",$nir_fishertest_list,"\n";
		print $sample1,"_vs_",$sample2,"\n",$ir_fishertest_list,"\n";
		print $sample1,"_vs_",$sample2,"\n",$nir_sample1_ratio,"\n";
		print $sample1,"_vs_",$sample2,"\n",$nir_sample2_ratio,"\n";

		# `cat /users/chend/work/R/RNA-Seq/Cor_single.r | R --slave --args $RPKM_file $log $index1 $index2 $outdir`;
		my $vs = $sample1."_vs_".$sample2;

		`mkdir -p $outdir/$vs`;	#make the out direcotry
		my $dir = &AbsolutePath("dir","$outdir/$vs");		#the absolute path of out directory
		# nir ras
		&make_ras_test_script($nir_fishertest_list,$nir_sample1_ratio,$nir_sample2_ratio,"_conir.r",$dir);
		print "nir ras cmd: cat $dir/_conir.r | R --slave --args $conir_file _$vs\_conir_test_result.txt $dir \n";
		`cat $dir/_conir.r | R --slave --args $conir_file _$vs\_conir_test_result.txt $dir `;
		# ir ras
		&make_ras_test_script($ir_fishertest_list,$ir_sample1_ratio,$ir_sample2_ratio,"_coir.r",$dir);
		print "ir ras cmd: cat $dir/_coir.r | R --slave --args $coir_file _$vs\_coir_test_result.txt $dir \n";
		`cat $dir/_coir.r | R --slave --args $coir_file _$vs\_coir_test_result.txt $dir `;
		chdir($dir);
		my $pf = $ConfigOption{"pvalue"};
		my $df = $ConfigOption{"diff"};
		`perl /users/ablife/ablife-perl/public/get_sub_file_filter_by_titleitem.pl _$vs\_conir_test_result.txt 'ID|Chr|Ass|Strand|Mss|Type|Isoform|KorN|$sample1|$sample2|Diff|Pvalue' > _$vs\_nir_test_result.txt && cat _$vs\_nir_test_result.txt | perl -ne 'chomp;if(/^X\\.ID/){print \$_,"\\tDiffType\\n";next;}\@line=split(/\\t/);if(\$line[-1]<=$pf && \$line[-2]<=-$df){print \$_,"\\tdown\\n";}if(\$line[-1]<=$pf && \$line[-2]>=$df){print \$_,"\\tup\\n";}' > $vs\_NIR_RAS_p$pf\_d$df.txt `;
		`perl /users/ablife/ablife-perl/public/get_sub_file_filter_by_titleitem.pl _$vs\_coir_test_result.txt 'ID|Chr|Ass|Strand|Mss|Type|Isoform|KorN|$sample1|$sample2|Diff|Pvalue' > _$vs\_ir_test_result.txt && cat _$vs\_ir_test_result.txt | perl -ne 'chomp;if(/^X\\.ID/){print \$_,"\\tDiffType\\n";next;}\@line=split(/\\t/);if(\$line[-1]<=$pf && \$line[-2]<=-$df){print \$_,"\\tdown\\n";}if(\$line[-1]<=$pf && \$line[-2]>=$df){print \$_,"\\tup\\n";}' > $vs\_IR_RAS_p$pf\_d$df.txt `;
	}
	chdir($outdir);
}

sub make_ras_test_script {
	my ($list,$ratio1,$ratio2,$rfilename,$dir)=@_;
	my $r_file = $dir."/$rfilename";
	open ROUT ,">$r_file" or die "can't open $r_file!\n";
	print ROUT <<EOF;
rm(list=ls())
library('edgeR')
args <- commandArgs(trailingOnly = TRUE)	#get the arguments from the command line
if (length(args) < 3) {
	stop('Your input arguments is wrong!\\n
		args1:\\t co_as file\\n
		args2:\\t out  file\\n
		args3:\\t out directory\\n'
	)
} else {
	infile <- args[1]
	outfile <- args[2]
	outdir <- args[3]
}

setwd(outdir)
Data <- read.delim(file=infile,header=T,stringsAsFactors=FALSE)
Pvalue <-c()
Diff <- c()
for (i in 1:nrow(Data)) {
	diff <- 0
	x <-c()
	xnum <- 0
	for (j in c($list)) {
	    if(!is.na(as.numeric(Data[i,j]))){
	      x <- c(x,as.numeric(Data[i,j]))
	      xnum <- xnum+1
	    }
	}
	if(!is.null(x) && xnum==4){
		compare <- matrix(x,nr=2)
		compare <- round(compare)
		result <- fisher.test(compare)		#,alternative = "greater"
		Pvalue <- c(Pvalue,result\$p.value)
		diff <- as.numeric(Data[i,$ratio1])-as.numeric(Data[i,$ratio2])
		Diff <- c(Diff,diff)
	}else{
	    Pvalue <- c(Pvalue,1)
	    Diff <- c(Diff,diff)
	}
}

Data <- cbind(Data,Diff,Pvalue)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')
EOF
}

sub Co_RAS_analysis {    #&Select_co_RAS($method,$outdir,$file1,$file2,$i,$j,\%RAS_method)
	my ($option_ref, $outdir) = @_;
	chdir($outdir);
	`mkdir -p co_RAS` if ( !-d "co_RAS"); 
	my $coRAS_dir = $outdir."/co_RAS";
	foreach my $coRAS (sort keys %{$option_ref->{"CORAS"}}) {
		my @nir_file = ();
		my @ir_file = ();
		my @groups = ();
		chdir($outdir);
		#print join(":",@{$option_ref->{"CORAS"}->{$coRAS}}),"\n";
		#print `pwd`;
		$outdir = &AbsolutePath("dir",$outdir);
		foreach my $group (@{$option_ref->{"CORAS"}->{$coRAS}}) {
			push @groups,$group;
			opendir(DH,$outdir) or die "$!\n";
			foreach my $element (readdir DH) {
				if (-d $element) {
					if ($element =~ /$group/i ) {
						opendir(DIR,$element) or die "$!\n";
						# print $outdir."/".$element,"\n";
						foreach my $file (readdir DIR) {
							next if $file =~ /^_/;
							print $file,"\n";
							if ($file =~ /NIR_RAS/) {
								push @nir_file, $outdir."/".$element."/".$file;
								#print $outdir."/".$element."/".$file,"\n";
							}elsif ($file =~ /IR_RAS/i) {
								push @ir_file, $outdir."/".$element."/".$file;
								#print $outdir."/".$element."/".$file,"\n";
							} else {
								next;
							}
						}
					}
				} else {
					next;
				}
			}
			close(DH);
		}
		print join("\n",@nir_file),"\n";
		print join("\n",@ir_file),"\n";
		chdir("co_RAS");
		my $co_dir = join("_and_",@{$option_ref->{"CORAS"}->{$coRAS}})."_"."co_RAS";
		`mkdir -p $co_dir` if ( !-d $co_dir);
		$co_dir = &AbsolutePath("dir",$co_dir); chdir($co_dir);
		print `pwd`;
		
		my $coRAS_name = join("_and_",@{$option_ref->{"CORAS"}->{$coRAS}})."_"."co_RAS";
		
		&Do_coras_nir(\@nir_file,\@groups,$conir,$coRAS_name);
		&Do_coras_ir(\@ir_file,\@groups,$coir,$coRAS_name);
		
		# my $venn_f = join(",",@ir_file).",".join(",",@Down_file);
		# `Rscript /users/ablife/ablife-R/Venn/latest/Venn.r -f $venn_f -n $coRAS_name\_co_RAS_venn`;
	}
	`cd $coRAS_dir && perl /users/ablife/ablife-perl/public/pipeline_stat.pl -i ./ -p _co_RAS -m coras.log -o coRAS_statics.xls`;
	chdir($outdir);
}

sub Do_coras_nir {
	my ($nirs,$vsgroups,$conir_file,$vs) = @_;
	my $nir_number = scalar(@{$nirs});
	my %co_up_number = ();
	my %co_down_number = ();
	my %diffinfo = ();
	my $i =0;
	#记录每个nir ras中的up和down的asid出现的个数，并记录pvalue，tvalue，difftype信息。
	foreach my $nir_ras_file (@{$nirs}){
		open IN,$nir_ras_file;
		my $upidfile = $vsgroups->[$i]."_nir_up";
		my $downidfile = $vsgroups->[$i]."_nir_down";
		while(<IN>){
			chomp;
			next if(/^X\.ID/);
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				$co_up_number{$line[0]} = 0 if not defined($co_up_number{$line[0]});
				$co_up_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}="$line[-3]\t$line[-2]\t$line[-1]";
				`echo $line[0] >> $upidfile`;
			}
			if ($line[-1] eq "down"){
				$co_down_number{$line[0]} = 0 if not defined($co_down_number{$line[0]});
				$co_down_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}="$line[-3]\t$line[-2]\t$line[-1]";
				`echo $line[0] >> $downidfile`;
			}
		}
		close IN;
		$i++;
	}
	my $venn_f = `ls *_nir_up | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn_notitle.r -f $venn_f -n coRAS_NIR_up_venn`;
	$venn_f = `ls *_nir_down | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn_notitle.r -f $venn_f -n coRAS_NIR_down_venn`;

	#遍历%co_up_number和%co_down_number，如果asid出现的次数等于组的个数，则输出。
	open UP,">coup_nir.txt";
	open DOWN,">codown_nir.txt";
	open LOG,">coras.log";
	foreach my $id (sort keys %co_up_number){
		if ($co_up_number{$id} < $nir_number){
			delete $co_up_number{$id};
		}
	}
	foreach my $id (sort keys %co_down_number){
		if ($co_down_number{$id} < $nir_number){
			delete $co_down_number{$id};
		}
	}
	#读取conir_file，分别找出up和down的结果
	open NIR,$conir_file;
	while(<NIR>){
		chomp;
		if(/^#ID/){
			print UP $_;
			print DOWN $_;
			foreach my $group (@{$vsgroups}){
				print UP "\t$group:Pvalue\t$group:Tvalue\t$group:DiffType";
				print DOWN "\t$group:Pvalue\t$group:Tvalue\t$group:DiffType";
			}
			print UP "\n";
			print DOWN "\n";
			next;
		}
		my @line=split(/\t/);
		if (defined($co_up_number{$line[0]})){
			print UP $_;
			foreach my $group (@{$vsgroups}){
				print UP "\t",$diffinfo{$line[0]}{$group};
			}
			print UP "\n";
		}
		if (defined($co_down_number{$line[0]})){
			print DOWN $_;
			foreach my $group (@{$vsgroups}){
				print DOWN "\t",$diffinfo{$line[0]}{$group};
			}
			print DOWN "\n";
		}
	}
	close NIR;
	close UP;
	close DOWN;
	#输出统计结果
	my $coupnum = scalar(keys %co_up_number);
	my $codownnum = scalar(keys %co_down_number);
	print LOG "NIR coup as:$coupnum\n";
	print LOG "NIR codown as:$codownnum\n";
	close LOG;

}


sub Do_coras_ir {
	my ($irs,$vsgroups,$coir_file,$vs) = @_;
	my $ir_number = scalar(@{$irs});
	my %co_up_number = ();
	my %co_down_number = ();
	my %diffinfo = ();
	my $i =0;
	#记录每个ir ras中的up和down的asid出现的个数，并记录pvalue，tvalue，difftype信息。
	foreach my $ir_ras_file (@{$irs}){
		open IN,$ir_ras_file;
		my $upidfile = $vsgroups->[$i]."_ir_up";
		my $downidfile = $vsgroups->[$i]."_ir_down";
		while(<IN>){
			chomp;
			next if(/^X\.ID/);
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				$co_up_number{$line[0]} = 0 if not defined($co_up_number{$line[0]});
				$co_up_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}="$line[-3]\t$line[-2]\t$line[-1]";
				`echo $line[0] >> $upidfile`;
			}
			if ($line[-1] eq "down"){
				$co_down_number{$line[0]} = 0 if not defined($co_down_number{$line[0]});
				$co_down_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}="$line[-3]\t$line[-2]\t$line[-1]";
				`echo $line[0] >> $downidfile`;
			}
		}
		close IN;
		$i++;
	}
	my $venn_f = `ls *_ir_up | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn_notitle.r -f $venn_f -n coRAS_IR_up_venn`;
	$venn_f = `ls *_ir_down | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn_notitle.r -f $venn_f -n coRAS_IR_down_venn`;

	#遍历%co_up_number和%co_down_number，如果asid出现的次数等于组的个数，则输出。
	open UP,">coup_ir.txt";
	open DOWN,">codown_ir.txt";
	open LOG,">>coras.log";
	foreach my $id (sort keys %co_up_number){
		if ($co_up_number{$id} < $ir_number){
			delete $co_up_number{$id};
		}
	}
	foreach my $id (sort keys %co_down_number){
		if ($co_down_number{$id} < $ir_number){
			delete $co_down_number{$id};
		}
	}
	#读取coir_file，分别找出up和down的结果
	open IR,$coir_file;
	while(<IR>){
		chomp;
		if(/^#ID/){
			print UP $_;
			print DOWN $_;
			foreach my $group (@{$vsgroups}){
				print UP "\t$group:Pvalue\t$group:Tvalue\t$group:DiffType";
				print DOWN "\t$group:Pvalue\t$group:Tvalue\t$group:DiffType";
			}
			print UP "\n";
			print DOWN "\n";
			next;
		}
		my @line=split(/\t/);
		if (defined($co_up_number{$line[0]})){
			print UP $_;
			foreach my $group (@{$vsgroups}){
				print UP "\t",$diffinfo{$line[0]}{$group};
			}
			print UP "\n";
		}
		if (defined($co_down_number{$line[0]})){
			print DOWN $_;
			foreach my $group (@{$vsgroups}){
				print DOWN "\t",$diffinfo{$line[0]}{$group};
			}
			print DOWN "\n";
		}
	}
	close IR;
	close UP;
	close DOWN;
	#输出统计结果
	my $coupnum = scalar(keys %co_up_number);
	my $codownnum = scalar(keys %co_down_number);
	print LOG "IR coup as:$coupnum\n";
	print LOG "IR codown as:$codownnum\n";
	close LOG;

}

sub Di_RAS_analysis {    #&Select_co_RAS($method,$outdir,$file1,$file2,$i,$j,\%RAS_method)
	my ($option_ref, $outdir) = @_;
	chdir($outdir);
	`mkdir -p di_RAS` if ( !-d "di_RAS"); 
	my $diRAS_dir = $outdir."/di_RAS";
	foreach my $diRAS (sort keys %{$option_ref->{"DIRAS"}}) {
		my @nir_file = ();
		my @ir_file = ();
		my @groups = ();
		chdir($outdir);
		#print join(":",@{$option_ref->{"DIRAS"}->{$DIRAS}}),"\n";
		#print `pwd`;
		$outdir = &AbsolutePath("dir",$outdir);
		foreach my $group (@{$option_ref->{"DIRAS"}->{$diRAS}}) {
			push @groups,$group;
			opendir(DH,$outdir) or die "$!\n";
			foreach my $element (readdir DH) {
				if (-d $element) {
					if ($element =~ /$group/i ) {
						opendir(DIR,$element) or die "$!\n";
						# print $outdir."/".$element,"\n";
						foreach my $file (readdir DIR) {
							next if $file =~ /^_/;
							print $file,"\n";
							if ($file =~ /NIR_RAS/) {
								push @nir_file, $outdir."/".$element."/".$file;
								#print $outdir."/".$element."/".$file,"\n";
							}elsif ($file =~ /IR_RAS/i) {
								push @ir_file, $outdir."/".$element."/".$file;
								#print $outdir."/".$element."/".$file,"\n";
							} else {
								next;
							}
						}
					}
				} else {
					next;
				}
			}
			close(DH);
		}
		print join("\n",@nir_file),"\n";
		print join("\n",@ir_file),"\n";
		chdir("di_RAS");
		my $di_dir = join("_and_",@{$option_ref->{"DIRAS"}->{$diRAS}})."_"."di_RAS";
		`mkdir -p $di_dir` if ( !-d $di_dir);
		$di_dir = &AbsolutePath("dir",$di_dir); chdir($di_dir);
		print `pwd`;
		
		my $diRAS_name = join("_and_",@{$option_ref->{"DIRAS"}->{$diRAS}})."_"."di_RAS";
		
		&Do_diras_nir(\@nir_file,\@groups,$conir,$diRAS_name);
		&Do_diras_ir(\@ir_file,\@groups,$coir,$diRAS_name);

	}
	`cd $diRAS_dir && cat *_di_RAS/diras.log|sed -e 's/:\\s*/\\t/' > diRAS_statics.xls`;
	chdir($outdir);
}

sub Do_diras_nir {
	my ($nirs,$vsgroups,$conir_file,$vs) = @_;
	my $nir_number = scalar(@{$nirs});
	my %co_up_number = ();
	my %co_down_number = ();
	my %diffinfo = ();
	my $i =0;
	#记录每个nir ras中的up和down的asid出现的个数，并记录pvalue，tvalue，difftype信息。
	foreach my $nir_ras_file (@{$nirs}){
		open IN,$nir_ras_file;
		my $upidfile = $vsgroups->[$i]."_nir_up";
		my $downidfile = $vsgroups->[$i]."_nir_down";
		while(<IN>){
			chomp;
			next if(/^X\.ID/);
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				$co_up_number{$line[0]} = 0 if not defined($co_up_number{$line[0]});
				$co_up_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}=$_;
				`echo $line[0] >> $upidfile`;
			}
			if ($line[-1] eq "down"){
				$co_down_number{$line[0]} = 0 if not defined($co_down_number{$line[0]});
				$co_down_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}=$_;
				`echo $line[0] >> $downidfile`;
			}
		}
		close IN;
		$i++;
	}
	my $venn_f = `ls *_nir_up | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn.r -f $venn_f -n diRAS_NIR_up_venn`;
	$venn_f = `ls *_nir_down | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn.r -f $venn_f -n diRAS_NIR_down_venn`;

	#再读一遍，如果在一个样品中%co_up_number和%co_down_number为1，则为该样品的di
	open LOG,">diras.log";
	print LOG "#$vs\n";
	$i=0;
	foreach my $nir_ras_file (@{$nirs}){
		open DIUP,">".$vsgroups->[$i]."_diRas_up_nir.txt";
		open DIDOWN,">".$vsgroups->[$i]."_diRas_down_nir.txt";
		open IN,$nir_ras_file;
		my $up_n = 0;
		my $down_n = 0;
		while(<IN>){
			chomp;
			if(/^X\.ID/){
				print DIUP $_,"\n";
				print DIDOWN $_,"\n";
				next;
			}
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				if ($co_up_number{$line[0]}==1){
					print DIUP $_,"\n";
					$up_n++;
				}
			}
			if ($line[-1] eq "down"){
				if ($co_down_number{$line[0]}==1){
					print DIDOWN $_,"\n";
					$down_n++;
				}
			}
		}
		close IN;
		close DIUP;
		close DIDOWN;
		print LOG $vsgroups->[$i]," NIR diup as:$up_n\n";
		print LOG $vsgroups->[$i]," NIR didown as:$down_n\n";
		$i++;
	}
	close LOG;
}


sub Do_diras_ir {
	my ($irs,$vsgroups,$coir_file,$vs) = @_;
	my $ir_number = scalar(@{$irs});
	my %co_up_number = ();
	my %co_down_number = ();
	my %diffinfo = ();
	my $i =0;
	#记录每个ir ras中的up和down的asid出现的个数，并记录pvalue，tvalue，difftype信息。
	foreach my $ir_ras_file (@{$irs}){
		open IN,$ir_ras_file;
		my $upidfile = $vsgroups->[$i]."_ir_up";
		my $downidfile = $vsgroups->[$i]."_ir_down";
		while(<IN>){
			chomp;
			next if(/^X\.ID/);
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				$co_up_number{$line[0]} = 0 if not defined($co_up_number{$line[0]});
				$co_up_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}=$_;
				`echo $line[0] >> $upidfile`;
			}
			if ($line[-1] eq "down"){
				$co_down_number{$line[0]} = 0 if not defined($co_down_number{$line[0]});
				$co_down_number{$line[0]}++;
				$diffinfo{$line[0]}{$vsgroups->[$i]}=$_;
				`echo $line[0] >> $downidfile`;
			}
		}
		close IN;
		$i++;
	}
	my $venn_f = `ls *_ir_up | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn.r -f $venn_f -n diRAS_IR_up_venn`;
	$venn_f = `ls *_ir_down | cat | perl -ne 'chomp;print \$_,",";'|sed 's/,\$//'`;
	`Rscript /users/ablife/ablife-R/Venn/latest/Venn.r -f $venn_f -n diRAS_IR_down_venn`;

	#再读一遍，如果在一个样品中%co_up_number和%co_down_number为1，则为该样品的di
	open LOG,">>diras.log";
	$i=0;
	foreach my $ir_ras_file (@{$irs}){
		open DIUP,">".$vsgroups->[$i]."_diRas_up_ir.txt";
		open DIDOWN,">".$vsgroups->[$i]."_diRas_down_ir.txt";
		open IN,$ir_ras_file;
		my $up_n = 0;
		my $down_n = 0;
		while(<IN>){
			chomp;
			if(/^X\.ID/){
				print DIUP $_,"\n";
				print DIDOWN $_,"\n";
				next;
			}
			my @line =split (/\t/);
			if ($line[-1] eq "up"){
				if ($co_up_number{$line[0]}==1){
					print DIUP $_,"\n";
					$up_n++;
				}
			}
			if ($line[-1] eq "down"){
				if ($co_down_number{$line[0]}==1){
					print DIDOWN $_,"\n";
					$down_n++;
				}
			}
		}
		close IN;
		close DIUP;
		close DIDOWN;
		print LOG $vsgroups->[$i]," IR diup as:$up_n\n";
		print LOG $vsgroups->[$i]," IR didown as:$down_n\n";
		$i++;
	}
	close LOG;
}


sub select_index{
	my ($file,$column) = @_;
	my $index = 0;
	open IN,$file;
	my $title = <IN>;
	close IN;
	chomp($title);
	my @line=split(/\t/,$title);
	for(my $i=0;$i<scalar(@line);$i++){
		if($line[$i] eq $column){
			$index = $i;
			return $index;
		}
	}
	return -1;
}


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



###############Sub_format_datetime
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf( "%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec );
}
