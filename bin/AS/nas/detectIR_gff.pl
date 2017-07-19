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
GetOptions( \%opts, "gff=s", "sj=s", "allsj=s", "g=s", "o=s", "h" );

if (   !defined( $opts{gff} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{allsj} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-gff        annotation file        must be given;

		-sj         splice junction used to detect as events        must be given;

		-allsj         all splice junction         must be given;

		-g          use gff sj? 1 or 0     default is 0;

		-o          out file prefix               must be given;

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
# my $mod_file  = $opts{mod};
my $allsj_file  = $opts{allsj};
my $gffsj = 0;
$gffsj = $opts{g} if defined($opts{g});
my $out_file = $opts{o};
# open OUT_SE, ">$out_file\_SE" || die;
# open OUT_CE, ">$out_file\_CE" || die;
# open OUT_MXE, ">$out_file\_MXE" || die;
# open OUT_ALE, ">$out_file\_ALE" || die;
# open OUT_AFE, ">$out_file\_AFE" || die;
# open OUT_A3SS, ">$out_file\_A3SS" || die;
# open OUT_A5SS, ">$out_file\_A5SS" || die;
# open OUT_II, ">$out_file\_II" || die;
open OUT_IR, ">$out_file" || die;
# open OUT_EE, ">$out_file\_EE" || die;

###读取junctions.bed中sj的信息
my %negative_start = ();    #负向junction起始位点信息记录hash
my %positive_start = ();    #正向同上
my %negative_end   = ();    #负向junction结束位点信息记录hash
my %positive_end   = ();    #正向同上

my %positive_sj    = ();    #正向sj的reads数
my %negative_sj    = ();    #负向同上
if($gffsj==1){
	&load_gffsj($sj_file);
	&load_all_gffsj($allsj_file);
}elsif($gffsj==0){
	&load_sj($sj_file);
	&load_all_sj($allsj_file);
}else{
	die "-g must set to 0 or 1\n";
}


###读取gff的内容，记录每个转录本的splice site。
my %gff = ();
my %exon_hash = ();
my %exon_negative_start = ();    #负向exon起始位点信息记录hash
my %exon_positive_start = ();    #正向同上
my %exon_negative_end   = ();    #负向exon结束位点信息记录hash
my %exon_positive_end   = ();    #正向同上
&load_gff(\%gff,$gff_file);

my %usedjunc = ();


# my %gene_model = ();

# #将gene model保存进哈希以便查询
# open MOD, $mod_file || die;
# while (<MOD>) {
# 	chomp;
# 	next if ( $_ =~ /chrom/ );
# 	my @line = split;
# 	$gene_model{$line[-1]} = 1;
# }
# close MOD;


###以每个基因作为model找可变剪接事件。
foreach my $chr (sort keys %gff){
	foreach my $strand (sort keys %{$gff{$chr}}){
		foreach my $gene (sort keys %{$gff{$chr}->{$strand}}){
			my @thisgene = @{$gff{$chr}->{$strand}->{$gene}};
			# &find_se($chr,$strand,$gene,\@thisgene);   ##find skipped exon events
			# # print "done find SE in $gene\n";
			# &find_ce($chr,$strand,$gene,\@thisgene) if $gffsj==0;   ##find cassette exon events
			# # print "done find CE in $gene\n";
			# &find_mxe($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# # print "done find MXE in $gene\n";
			# &find_afe_ale($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# # print "done find AFE and ALE in $gene\n";
			# &find_ii($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find II in $gene\n";
			&find_ir($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find IR in $gene\n";
			# &find_ee($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find EE in $gene\n";
			# print "done $gene\n";
		}
	}
}

# ###以每个基因作为model找可变剪接事件(A3SS,A5SS)。
# foreach my $chr (sort keys %gff){
# 	foreach my $strand (sort keys %{$gff{$chr}}){
# 		foreach my $gene (sort keys %{$gff{$chr}->{$strand}}){
# 			my @thisgene = @{$gff{$chr}->{$strand}->{$gene}};
# 			&find_a3ss($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
# 			# print "done find A3SS in $gene\n";
# 			&find_a5ss($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
# 			# print "done find A5SS in $gene\n";
# 			# print "done $gene\n";
# 		}
# 	}
# }

# close OUT_SE;
# close OUT_CE;
# close OUT_A3SS;
# close OUT_A5SS;
# close OUT_MXE;
# close OUT_II;
close OUT_IR;
# close OUT_EE;
# close OUT_AFE;
# close OUT_ALE;


###子过程-----find_ir($chr,$strand,$gene,\@thisgene)
sub find_ir {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f<$#genesite;$f++){
			for(my $t = $genesite[2*$f];$t<=$genesite[2*$f+1];$t++){
				if ( defined( $positive_start{$chr}->{$t} ) ){
					my @start_junction = split /\s+/, $positive_start{$chr}->{$t};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position, $genename ) = (split /\:/, $i,2);
						if($position>=$genesite[2*$f] && $position<=$genesite[2*$f+1]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+1];
							# my $structure = "1\^2\-\,0";
							# my $splice_chain = "$t\^$position\-\,";
							# my $splice_position = "$t,$position";
							# my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							# print OUT_IR "$chr\t$gene\tIR\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							print OUT_IR "$chr\t$t\t$position\t$strand\tIntronR\t$genename\n";
						}
					}
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f<$#genesite;$f++){
			for(my $t = $genesite[2*$f];$t>=$genesite[2*$f+1];$t--){
				if ( defined( $negative_start{$chr}->{$t} ) ){
					my @start_junction = split /\s+/, $negative_start{$chr}->{$t};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position, $genename ) = (split /\:/, $i,2);
						if($position<=$genesite[2*$f] && $position>=$genesite[2*$f+1]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+1];
							# my $structure = "1\^2\-\,0";
							# my $splice_chain = "$t\^$position\-\,";
							# my $splice_position = "$t,$position";
							# my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							# print OUT_IR "$chr\t$gene\tIR\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							print OUT_IR "$chr\t$position\t$t\t$strand\tIntronR\t$genename\n";
						}
					}
				}
			}
		}
	}
}


###子过程-----load_sj($sj_file)
sub load_sj {
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
			$positive_sj{$chr}->{"$start:$end"} = $readsNum;
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
			$negative_sj{$chr}->{"$start:$end"} = $readsNum;
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
	print "done reading SJ ...........", "\n\n";
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
			$positive_sj{$chr}->{"$start:$end"} = $readsNum;
		}
		else {
			$start = $right - $rightsize + 1;
			$end   = $left + $leftsize;
			$negative_sj{$chr}->{"$start:$end"} = $readsNum;
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

###子过程-----load_gffsj($sj_file)
sub load_gffsj {
	my ($sj_file) = @_;
	open SJ, $sj_file || die;
	while (<SJ>) {              #把所有的junction添加进hash
		#gff的sj
		#chr20   199843  204678  -
		chomp;
		next if ( $_ =~ /^track/i );
		my ( $chr, $start, $end, $strand,$genename) = split /\t/ ;
		if ( $strand eq "+" ) {
			$positive_sj{$chr}->{"$start:$end"} = $genename;
			if ( !defined( $positive_start{$chr}->{$start} ) ) {
				$positive_start{$chr}->{$start} = "$end\:$genename";
			}
			else {
				$positive_start{$chr}->{$start} .= "\t$end\:$genename";

			}
			if ( !defined( $positive_end{$chr}->{$end} ) ) {
				$positive_end{$chr}->{$end} = "$start\:$genename";
			}
			else {
				$positive_end{$chr}->{$end} .= "\t$start\:$genename";

				#调试信息
				#print $positive_end{$chr}->{$end},"\n";
			}
		}
		else {
			$negative_sj{$chr}->{"$start:$end"} = $genename;
			if ( !defined( $negative_start{$chr}->{$start} ) ) {
				$negative_start{$chr}->{$start} = "$end\:$genename";
			}
			else {
				$negative_start{$chr}->{$start} .= "\t$end\:$genename";
			}
			if ( !defined( $negative_end{$chr}->{$end} ) ) {
				$negative_end{$chr}->{$end} = "$start\:$genename";
			}
			else {
				$negative_end{$chr}->{$end} .= "\t$start\:$genename";
			}
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

###子过程-----load_all_gffsj($sj_file)
sub load_all_gffsj {
	my ($sj_file) = @_;
	open SJ, $sj_file || die;
	while (<SJ>) {              #把所有的junction添加进hash
		#gff的sj
		#chr20   199843  204678  -
		chomp;
		next if ( $_ =~ /^track/i );
		my ( $chr, $start, $end, $strand) = split /\t/ ;
		if ( $strand eq "+" ) {
			$positive_sj{$chr}->{"$start:$end"} = 3;
		}
		else {
			$negative_sj{$chr}->{"$start:$end"} = 3;
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

#子过程-----load_gff(\%gff,$gff_file)
sub load_gff {
	my ($gff, $gff_file) = @_;
	my $gene = "";
	my $this_chr = "";
	my $this_strand = "";
	my $flag = 0;
	open GFF, $gff_file    || die;
	while (<GFF>) {
		chomp;
		next if ( $_ =~ /^\#/ );
		my (
		$chr, $source, $feature, $start, $end,
		$other1, $strand, $other2, $info
		)= split /\t/;
		next if ( $feature =~ /chromosome/ || $_ =~ /match/ || $_ =~ /region/);
		if($feature =~ /^gene$/){
			$flag = 0;
			next;
		}
		if ( $feature =~ /^(mRNA|miRNA|mRNA_TE_gene|ncRNA|rRNA|snoRNA|snRNA|tRNA|transcript)$/ ) {
			#将上一个gene的各位点排序，正向的按从小到大，负向的按从大到小。
			if(defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
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
			$flag = 1;
		}
		if ( $feature =~ /^exon$/ ){
			next if $flag == 0;
			#将exon的起始和终止位置放进gene数组。
			push @{$gff->{$this_chr}->{$this_strand}->{$gene}},$start;
			push @{$gff->{$this_chr}->{$this_strand}->{$gene}},$end;
			if($this_strand eq "+"){
				$exon_hash{$this_chr}->{$this_strand}->{"$start:$end"} = 1;
				if(!defined($exon_positive_start{$this_chr}->{$this_strand}->{$start})){
					$exon_positive_start{$this_chr}->{$this_strand}->{$start} = $end;
				}else{
					$exon_positive_start{$this_chr}->{$this_strand}->{$start} .= "\t$end";
				}
				if(!defined($exon_positive_end{$this_chr}->{$this_strand}->{$end})){
					$exon_positive_end{$this_chr}->{$this_strand}->{$end} = $start;
				}else{
					$exon_positive_end{$this_chr}->{$this_strand}->{$end} .= "\t$start";
				}

			}elsif($this_strand eq "-"){
				$exon_hash{$this_chr}->{$this_strand}->{"$end:$start"} = 1;
				if(!defined($exon_negative_start{$this_chr}->{$this_strand}->{$end})){
					$exon_negative_start{$this_chr}->{$this_strand}->{$end} = $start;
				}else{
					$exon_negative_start{$this_chr}->{$this_strand}->{$end} .= "\t$start";
				}
				if(!defined($exon_negative_end{$this_chr}->{$this_strand}->{$start})){
					$exon_negative_end{$this_chr}->{$this_strand}->{$start} = $end;
				}else{
					$exon_negative_end{$this_chr}->{$this_strand}->{$start} .= "\t$end";
				}
			}
		}
	}
	if($gene ne ""){
		if($this_strand eq "+"){
			@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$a<=>$b} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
		}else{
			@{$gff->{$this_chr}->{$this_strand}->{$gene}} = sort {$b<=>$a} @{$gff->{$this_chr}->{$this_strand}->{$gene}};
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