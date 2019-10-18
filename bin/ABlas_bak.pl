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
my $allsj_file  = $opts{allsj};
my $gffsj = 0;
$gffsj = $opts{g} if defined($opts{g});
my $out_file = $opts{o};
open OUT_SE, ">$out_file\_SE" || die;
open OUT_CE, ">$out_file\_CE" || die;
open OUT_MXE, ">$out_file\_MXE" || die;
open OUT_ALE, ">$out_file\_ALE" || die;
open OUT_AFE, ">$out_file\_AFE" || die;
open OUT_A3SS, ">$out_file\_A3SS" || die;
open OUT_A5SS, ">$out_file\_A5SS" || die;
open OUT_II, ">$out_file\_II" || die;
open OUT_IR, ">$out_file\_IR" || die;
open OUT_EE, ">$out_file\_EE" || die;

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


###以每个基因作为model找可变剪接事件。
foreach my $chr (sort keys %gff){
	foreach my $strand (sort keys %{$gff{$chr}}){
		foreach my $gene (sort keys %{$gff{$chr}->{$strand}}){
			my @thisgene = @{$gff{$chr}->{$strand}->{$gene}};
			&find_se($chr,$strand,$gene,\@thisgene);   ##find skipped exon events
			# print "done find SE in $gene\n";
			&find_ce($chr,$strand,$gene,\@thisgene) if $gffsj==0;   ##find cassette exon events
			# print "done find CE in $gene\n";
			&find_mxe($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find MXE in $gene\n";
			&find_afe_ale($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find AFE and ALE in $gene\n";
			&find_ii($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find II in $gene\n";
			&find_ir($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find IR in $gene\n";
			&find_ee($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find EE in $gene\n";
			# print "done $gene\n";
		}
	}
}

###以每个基因作为model找可变剪接事件(A3SS,A5SS)。
foreach my $chr (sort keys %gff){
	foreach my $strand (sort keys %{$gff{$chr}}){
		foreach my $gene (sort keys %{$gff{$chr}->{$strand}}){
			my @thisgene = @{$gff{$chr}->{$strand}->{$gene}};
			&find_a3ss($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find A3SS in $gene\n";
			&find_a5ss($chr,$strand,$gene,\@thisgene);   ##find Mutually exclusive exons events
			# print "done find A5SS in $gene\n";
			# print "done $gene\n";
		}
	}
}

close OUT_SE;
close OUT_CE;
close OUT_A3SS;
close OUT_A5SS;
close OUT_MXE;
close OUT_II;
close OUT_IR;
close OUT_EE;
close OUT_AFE;
close OUT_ALE;

###子过程-----find_ee($chr,$strand,$gene,\@thisgene)
sub find_ee {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	return 0 if @genesite < 6;
	if($strand eq "+"){
		for(my $f=0;2*$f<$#genesite-3;$f++){
			##exon truncation (both)
			if ( defined( $positive_start{$chr}->{$genesite[2*$f+1]} ) ){
				my @start_junction = split /\s+/, $positive_start{$chr}->{$genesite[2*$f+1]};
				#查看以该位点为起点的每个junction的end
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					if($position_left>$genesite[2*$f+2] && $position_left<$genesite[2*$f+3]){
						if ( defined( $positive_end{$chr}->{$genesite[2*$f+4]} ) ){
							my @end_junction = split /\s+/, $positive_end{$chr}->{$genesite[2*$f+4]};
							#查看以该位点为起点的每个junction的end
							foreach my $j (@end_junction) {
								my ( $position_right, $reads_right ) = split /\:/, $j;
								if($position_right>$genesite[2*$f+2] && $position_right<$genesite[2*$f+3]){
									next if $position_right<=$position_left;  ##有交叉
									if($gffsj==1){
										next if not defined($exon_hash{$chr}->{$strand}->{"$position_left:$position_right"});
									}
									my $flanks_left = $genesite[2*$f+1];
									my $flanks_right= $genesite[2*$f+4];
									my $flanks_exon_no = "$f,".($f+2);
									my $structure = "2\-3\^\,1\-4\^";
									my $splice_chain = "$position_left\-$position_right\^\,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
									my $splice_position = "$genesite[2*$f+2],$position_left,$position_right,$genesite[2*$f+3]";
									my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; splice_position \"$splice_position\"\;";
									print OUT_EE "$chr\t$gene\tEE\t$flanks_left\t$flanks_right\t$reads_left\t$strand\t$reads_right\t$info\n";
									my $junc = "$genesite[2*$f+1]:$position_left";
									$usedjunc{$chr}->{"+"}->{$junc} = 1;
									$junc = "$position_right:$genesite[2*$f+4]";
									$usedjunc{$chr}->{"+"}->{$junc} = 1;
								}
							}
						}
					}
				}
			}
			##exon extension (both)
			if($gffsj==0){
				if ( defined( $positive_start{$chr}->{$genesite[2*$f+1]} ) ){
					my @start_junction = split /\s+/, $positive_start{$chr}->{$genesite[2*$f+1]};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position_left, $reads_left ) = split /\:/, $i;
						if($position_left<$genesite[2*$f+2]){
							if ( defined( $positive_end{$chr}->{$genesite[2*$f+4]} ) ){
								my @end_junction = split /\s+/, $positive_end{$chr}->{$genesite[2*$f+4]};
								#查看以该位点为起点的每个junction的end
								foreach my $j (@end_junction) {
									my ( $position_right, $reads_right ) = split /\:/, $j;
									if($position_right>$genesite[2*$f+3]){
										next if $position_right-$position_left>400;  ##扩展后的exon大于400过滤掉
										my $flanks_left = $genesite[2*$f+1];
										my $flanks_right= $genesite[2*$f+4];
										my $flanks_exon_no = "$f,".($f+2);
										my $structure = "1\-4\^\,2\-3\^";
										my $splice_chain = "$position_left\-$position_right\^\,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
										my $splice_position = "$position_left,$genesite[2*$f+2],$genesite[2*$f+3],$position_right";
										my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; splice_position \"$splice_position\"\;";
										print OUT_EE "$chr\t$gene\tEE\t$flanks_left\t$flanks_right\t$reads_left\t$strand\t$reads_right\t$info\n";
										my $junc = "$genesite[2*$f+1]:$position_left";
										$usedjunc{$chr}->{"+"}->{$junc} = 1;
										$junc = "$position_right:$genesite[2*$f+4]";
										$usedjunc{$chr}->{"+"}->{$junc} = 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f<$#genesite-3;$f++){
			##exon truncation (both)
			if ( defined( $negative_start{$chr}->{$genesite[2*$f+1]} ) ){
				my @start_junction = split /\s+/, $negative_start{$chr}->{$genesite[2*$f+1]};
				#查看以该位点为起点的每个junction的end
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					if($position_left<$genesite[2*$f+2] && $position_left>$genesite[2*$f+3]){
						if ( defined( $negative_end{$chr}->{$genesite[2*$f+4]} ) ){
							my @end_junction = split /\s+/, $negative_end{$chr}->{$genesite[2*$f+4]};
							#查看以该位点为起点的每个junction的end
							foreach my $j (@end_junction) {
								my ( $position_right, $reads_right ) = split /\:/, $j;
								if($position_right<$genesite[2*$f+2] && $position_right>$genesite[2*$f+3]){
									next if $position_right>=$position_left;  ##有交叉
									if($gffsj==1){
										next if not defined($exon_hash{$chr}->{$strand}->{"$position_left:$position_right"});
									}
									my $flanks_left = $genesite[2*$f+1];
									my $flanks_right= $genesite[2*$f+4];
									my $flanks_exon_no = "$f,".($f+2);
									my $structure = "2\-3\^\,1\-4\^";
									my $splice_chain = "$position_left\-$position_right\^\,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
									my $splice_position = "$genesite[2*$f+2],$position_left,$position_right,$genesite[2*$f+3]";
									my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; splice_position \"$splice_position\"\;";
									print OUT_EE "$chr\t$gene\tEE\t$flanks_left\t$flanks_right\t$reads_left\t$strand\t$reads_right\t$info\n";
									my $junc = "$genesite[2*$f+1]:$position_left";
									$usedjunc{$chr}->{"-"}->{$junc} = 1;
									$junc = "$position_right:$genesite[2*$f+4]";
									$usedjunc{$chr}->{"-"}->{$junc} = 1;
								}
							}
						}
					}
				}
			}
			##exon extension (both)
			if($gffsj==0){
				if ( defined( $negative_start{$chr}->{$genesite[2*$f+1]} ) ){
					my @start_junction = split /\s+/, $negative_start{$chr}->{$genesite[2*$f+1]};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position_left, $reads_left ) = split /\:/, $i;
						if($position_left>$genesite[2*$f+2]){
							if ( defined( $negative_end{$chr}->{$genesite[2*$f+4]} ) ){
								my @end_junction = split /\s+/, $negative_end{$chr}->{$genesite[2*$f+4]};
								#查看以该位点为起点的每个junction的end
								foreach my $j (@end_junction) {
									my ( $position_right, $reads_right ) = split /\:/, $j;
									if($position_right<$genesite[2*$f+3]){
										next if $position_left-$position_right>400;  ##扩展后的exon大于400过滤掉
										my $flanks_left = $genesite[2*$f+1];
										my $flanks_right= $genesite[2*$f+4];
										my $flanks_exon_no = "$f,".($f+2);
										my $structure = "1\-4\^\,2\-3\^";
										my $splice_chain = "$position_left\-$position_right\^\,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
										my $splice_position = "$position_left,$genesite[2*$f+2],$genesite[2*$f+3],$position_right";
										my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; splice_position \"$splice_position\"\;";
										print OUT_EE "$chr\t$gene\tEE\t$flanks_left\t$flanks_right\t$reads_left\t$strand\t$reads_right\t$info\n";
										my $junc = "$genesite[2*$f+1]:$position_left";
										$usedjunc{$chr}->{"-"}->{$junc} = 1;
										$junc = "$position_right:$genesite[2*$f+4]";
										$usedjunc{$chr}->{"-"}->{$junc} = 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

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
						my ( $position, $reads ) = split /\:/, $i;
						if($position>=$genesite[2*$f] && $position<=$genesite[2*$f+1]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+1];
							my $structure = "1\^2\-\,0";
							my $splice_chain = "$t\^$position\-\,";
							my $splice_position = "$t,$position";
							my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							print OUT_IR "$chr\t$gene\tIR\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
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
						my ( $position, $reads ) = split /\:/, $i;
						if($position<=$genesite[2*$f] && $position>=$genesite[2*$f+1]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+1];
							my $structure = "1\^2\-\,0";
							my $splice_chain = "$t\^$position\-\,";
							my $splice_position = "$t,$position";
							my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							print OUT_IR "$chr\t$gene\tIR\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
						}
					}
				}
			}
		}
	}
}

###子过程-----find_ii($chr,$strand,$gene,\@thisgene)
sub find_ii {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			##intron extension (both)
			for(my $t = $genesite[2*$f]+1;$t<$genesite[2*$f+1];$t++){
				if ( defined( $positive_start{$chr}->{$t} ) ){
					my @start_junction = split /\s+/, $positive_start{$chr}->{$t};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position, $reads ) = split /\:/, $i;
						if($position>$genesite[2*$f+2] && $position<$genesite[2*$f+3]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+3];
							my $structure = "1\^4\-\,2\^3\-";
							my $splice_chain = "$t\^$position\-\,$genesite[2*$f+1]\^$genesite[2*$f+2]\-";
							my $splice_position = "$t,$genesite[2*$f+1],$genesite[2*$f+2],$position";
							my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							print OUT_II "$chr\t$gene\tII\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
						}
					}
				}
			}
			##intron truncation (both)
			if($gffsj==0){
				for(my $t = $genesite[2*$f+1]+1;$t<$genesite[2*$f+2];$t++){
					last if($t-$genesite[2*$f+1]>400);  #两侧exon延长后也应该不能过长，这里限制为400
					if ( defined( $positive_start{$chr}->{$t} ) ){
						my @start_junction = split /\s+/, $positive_start{$chr}->{$t};
						#查看以该位点为起点的每个junction的end
						foreach my $i (@start_junction) {
							my ( $position, $reads ) = split /\:/, $i;
							if($position>$genesite[2*$f+1] && $position<$genesite[2*$f+2] && $genesite[2*$f+2]-$position<400){
								my $flanks_left = $genesite[2*$f];
								my $flanks_right= $genesite[2*$f+3];
								my $structure = "2\^3\-\,1\^4\-";
								my $splice_chain = "$t\^$position\-\,$genesite[2*$f+1]\^$genesite[2*$f+2]\-";
								my $splice_position = "$genesite[2*$f+1],$t,$position,$genesite[2*$f+2]";
								my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
								print OUT_II "$chr\t$gene\tII\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							}
						}
					}
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			##intron extension (both)
			for(my $t = $genesite[2*$f]-1;$t>$genesite[2*$f+1];$t--){
				if ( defined( $negative_start{$chr}->{$t} ) ){
					my @start_junction = split /\s+/, $negative_start{$chr}->{$t};
					#查看以该位点为起点的每个junction的end
					foreach my $i (@start_junction) {
						my ( $position, $reads ) = split /\:/, $i;
						if($position<$genesite[2*$f+2] && $position>$genesite[2*$f+3]){
							my $flanks_left = $genesite[2*$f];
							my $flanks_right= $genesite[2*$f+3];
							my $structure = "1\^4\-\,2\^3\-";
							my $splice_chain = "$t\^$position\-\,$genesite[2*$f+1]\^$genesite[2*$f+2]\-";
							my $splice_position = "$t,$genesite[2*$f+1],$genesite[2*$f+2],$position";
							my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
							print OUT_II "$chr\t$gene\tII\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
						}
					}
				}
			}
			##intron truncation (both)
			if($gffsj==0){
				for(my $t = $genesite[2*$f+1]-1;$t>$genesite[2*$f+2];$t--){
					last if($genesite[2*$f+1]-$t>400);  #两侧exon延长后也应该不能过长，这里限制为400
					if ( defined( $negative_start{$chr}->{$t} ) ){
						my @start_junction = split /\s+/, $negative_start{$chr}->{$t};
						#查看以该位点为起点的每个junction的end
						foreach my $i (@start_junction) {
							my ( $position, $reads ) = split /\:/, $i;
							if($position<$genesite[2*$f+1] && $position>$genesite[2*$f+2] && $position-$genesite[2*$f+2]<400 ){
								my $flanks_left = $genesite[2*$f];
								my $flanks_right= $genesite[2*$f+3];
								my $structure = "2\^3\-\,1\^4\-";
								my $splice_chain = "$t\^$position\-\,$genesite[2*$f+1]\^$genesite[2*$f+2]\-";
								my $splice_position = "$genesite[2*$f+1],$t,$position,$genesite[2*$f+2]";
								my $info = "flanks \"$flanks_left\-\,$flanks_right\^\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; splice_position \"$splice_position\"\;";
								print OUT_II "$chr\t$gene\tII\t$flanks_left\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							}
						}
					}
				}
			}
		}
	}
}


sub reverse_structure {
	my ($str) = @_;
	my $r_str = "";
	if($str=~m/^(\d+)\S+(\d+)\-$/){
		my $left = $1;
		my $right = $2;
		my $max = $left>$right?$left:$right;
		if($left>$right){
			$r_str="1\^,2\^";
			for(my $i=3;$i<=$max;$i++){
				if($i%2==1){
					$r_str.="$i\-";
				}else{
					$r_str.="$i\^";
				}
			}
		}else{
			$r_str="2\^,1\^";
			for(my $i=3;$i<=$max;$i++){
				if($i%2==1){
					$r_str.="$i\-";
				}else{
					$r_str.="$i\^";
				}
			}
		}
	}
	if($str=~m/^(\d+)\S+(\d+)\[$/){
		my $left = $1;
		my $right = $2;
		my $max = $left>$right?$left:$right;
		if($left>$right){
			$r_str="1\^,2\[";
			for(my $i=3;$i<=$max;$i++){
				if($i%2==1){
					$r_str.="$i\^";
				}else{
					$r_str.="$i\-";
				}
			}
		}
	}
	return $r_str;
}

###子过程-----find_a5ss($chr,$strand,$gene,\@thisgene)
sub find_a5ss {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=1;2*$f<$#genesite;$f++){
			my $junc_end = $genesite[2*$f];
			if ( defined( $positive_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $positive_end{$chr}->{$junc_end};
				#查看以该位点为起点的每个junction的end
				LABEL1:foreach my $i (@start_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					my $junc = "$position:$junc_end";
					next if defined($usedjunc{$chr}->{"+"}->{$junc});
					next if $position==$genesite[2*$f-1];
					my $flanks_left = "null";
					my $flanks_right= $junc_end;
					my $flank_exon_start="$genesite[2*$f]";
					my $flank_exon_end="$genesite[2*$f+1]";
					my $splice_site = 1;
					my $structure = "";
					my $splice_chain = "";
					my $splice_position = "";
					my $alt_exon_start = "";
					for(my $m=$f-1;2*$m>=0;$m--){
						if($genesite[2*$m+1]<=$position){
							$structure = "$splice_site\-\,".$structure.($splice_site+1)."\-";
							$structure = &reverse_structure($structure);
							$splice_chain = "$position\^\,".$genesite[2*$m+1]."\^".$splice_chain;
							$splice_position = "$genesite[2*$m+1],$position,".$splice_position;
							$alt_exon_start = $genesite[2*$m];
							$splice_position=~s/,$//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
							print OUT_A5SS "$chr\t$gene\tA5SS\t$genesite[2*$m+1]\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;

						}elsif($genesite[2*$m]<=$position){
							$structure .= "$splice_site\-";
							$structure = ($splice_site+1)."\-\,".$structure;
							$structure = &reverse_structure($structure);
							$splice_chain = "$genesite[2*$m+1]\^".$splice_chain;
							$splice_chain = $position."\^\,".$splice_chain;
							$splice_position = "$position,$genesite[2*$m+1],".$splice_position;
							$alt_exon_start = $genesite[2*$m];
							$splice_position=~s/,$//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
							print OUT_A5SS "$chr\t$gene\tA5SS\t$position\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;
						}
						if(2*$m==0){
							$structure .= "$splice_site\-".($splice_site+1)."\[";
							$splice_chain = "$genesite[2*$m]\[$genesite[2*$m+1]\^".$splice_chain;
							$splice_position = "$genesite[2*$m],$genesite[2*$m+1],".$splice_position;
						}else{
							$structure .= "$splice_site\-".($splice_site+1)."\^";
							$splice_chain = "$genesite[2*$m]\-$genesite[2*$m+1]\^".$splice_chain;
							$splice_position = "$genesite[2*$m],$genesite[2*$m+1],".$splice_position;
						}
						
						$splice_site = $splice_site+2;
					}
					$structure = "$splice_site\-\,".$structure;
					$structure = &reverse_structure($structure);
					$splice_chain = "$position\^\,".$splice_chain;
					$splice_position = "$position,".$splice_position;
					$alt_exon_start = $genesite[0];
					$splice_position=~s/,$//;
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
					print OUT_A5SS "$chr\t$gene\tA5SS\t$position\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
					next LABEL1;
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=1;2*$f<$#genesite;$f++){
			my $junc_end = $genesite[2*$f];
			if ( defined( $negative_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $negative_end{$chr}->{$junc_end};
				#查看以该位点为起点的每个junction的end
				LABEL1:foreach my $i (@start_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					my $junc = "$position:$junc_end";
					next if defined($usedjunc{$chr}->{"-"}->{$junc});
					next if $position==$genesite[2*$f-1];
					my $flanks_left = "null";
					my $flanks_right= $junc_end;
					my $flank_exon_start="$genesite[2*$f]";
					my $flank_exon_end="$genesite[2*$f+1]";
					my $splice_site = 1;
					my $structure = "";
					my $splice_chain = "";
					my $splice_position = "";
					my $alt_exon_start = "";
					for(my $m=$f-1;2*$m>=0;$m--){
						if($genesite[2*$m+1]>=$position){
							$structure = "$splice_site\-\,".$structure.($splice_site+1)."\-";
							$structure = &reverse_structure($structure);
							$splice_chain = "$position\^\,".$genesite[2*$m+1]."\^".$splice_chain;
							$splice_position = "$genesite[2*$m+1],$position,".$splice_position;
							$alt_exon_start = $genesite[2*$m];
							$splice_position=~s/,$//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
							print OUT_A5SS "$chr\t$gene\tA5SS\t$genesite[2*$m+1]\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;

						}elsif($genesite[2*$m]>=$position){
							$structure .= "$splice_site\-";
							$structure = ($splice_site+1)."\-\,".$structure;
							$structure = &reverse_structure($structure);
							$splice_chain = "$genesite[2*$m+1]\^".$splice_chain;
							$splice_chain = $position."\^\,".$splice_chain;
							$splice_position = "$position,$genesite[2*$m+1],".$splice_position;
							$alt_exon_start = $genesite[2*$m];
							$splice_position=~s/,$//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
							print OUT_A5SS "$chr\t$gene\tA5SS\t$position\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;
						}
						if(2*$m==0){
							$structure .= "$splice_site\-".($splice_site+1)."\[";
							$splice_chain = "$genesite[2*$m]\[$genesite[2*$m+1]\^".$splice_chain;
							$splice_position = "$genesite[2*$m],$genesite[2*$m+1],".$splice_position;
						}else{
							$structure .= "$splice_site\-".($splice_site+1)."\^";
							$splice_chain = "$genesite[2*$m]\-$genesite[2*$m+1]\^".$splice_chain;
							$splice_position = "$genesite[2*$m],$genesite[2*$m+1],".$splice_position;
						}
						
						$splice_site = $splice_site+2;
					}
					$structure = "$splice_site\-\,".$structure;
					$structure = &reverse_structure($structure);
					$splice_chain = "$position\^\,".$splice_chain;
					$splice_position = "$position,".$splice_position;
					$alt_exon_start = $genesite[0];
					$splice_position=~s/,$//;
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; splice_position \"$splice_position\"\;";
					print OUT_A5SS "$chr\t$gene\tA5SS\t$position\t$flanks_right\t\.\t$strand\t$reads\t$info\n";
					next LABEL1;
				}
			}
		}
	}
}


###子过程-----find_a3ss($chr,$strand,$gene,\@thisgene)
sub find_a3ss {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			my $junc_start = $genesite[2*$f+1];
			if ( defined( $positive_start{$chr}->{$junc_start} ) ){
				my @end_junction = split /\s+/, $positive_start{$chr}->{$junc_start};
				#查看以该位点为起点的每个junction的end
				LABEL1:foreach my $i (@end_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					my $junc = "$junc_start:$position";
					next if defined($usedjunc{$chr}->{"+"}->{$junc});
					next if $position==$genesite[2*$f+2];
					my $flanks_left = $junc_start;
					my $flanks_right= "null";
					my $flank_exon_start="$genesite[2*$f]";
					my $flank_exon_end="$genesite[2*$f+1]";
					my $splice_site = 1;
					my $structure = "";
					my $splice_chain = "";
					my $splice_position = "";
					my $alt_exon_end = "";
					for(my $m=$f+1;2*$m<=$#genesite;$m++){
						if($genesite[2*$m]>=$position){
							$structure = "$splice_site\-\,".$structure.($splice_site+1)."\-";
							$splice_chain = "$position\-\,".$splice_chain.$genesite[2*$m]."\-";
							$splice_position = $splice_position.",$position,$genesite[2*$m]";
							$alt_exon_end = $genesite[2*$m+1];
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$genesite[2*$m]\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;

						}elsif($genesite[2*$m+1]>=$position){
							$structure .= "$splice_site\-";
							$structure = ($splice_site+1)."\-\,".$structure;
							$splice_chain .= "$genesite[2*$m]\-";
							$splice_chain = $position."\-\,".$splice_chain;
							$splice_position = $splice_position.",$genesite[2*$m],$position";
							$alt_exon_end = $genesite[2*$m+1];
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$position\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;
						}
						if(2*$m+1==$#genesite){
							$structure .= "$splice_site\-".($splice_site+1)."\]";
							$splice_chain .= "$genesite[2*$m]\-$genesite[2*$m+1]\]";
							$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
						}else{
							$structure .= "$splice_site\-".($splice_site+1)."\^";
							$splice_chain .= "$genesite[2*$m]\-$genesite[2*$m+1]\^";
							$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
						}
						
						$splice_site = $splice_site+2;
					}
					$structure = "$splice_site\-\,".$structure;
					$splice_chain = "$position\-\,".$splice_chain;
					$splice_position = $splice_position.",$position";
					$alt_exon_end = $genesite[-1];
					$splice_position=~s/^,//;
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$position\t\.\t$strand\t$reads\t$info\n";
					next LABEL1;
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			my $junc_start = $genesite[2*$f+1];
			if ( defined( $negative_start{$chr}->{$junc_start} ) ){
				my @end_junction = split /\s+/, $negative_start{$chr}->{$junc_start};
				#查看以该位点为起点的每个junction的end
				LABEL1:foreach my $i (@end_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					my $junc = "$junc_start:$position";
					next if defined($usedjunc{$chr}->{"-"}->{$junc});
					next if $position==$genesite[2*$f+2];
					my $flanks_left = $junc_start;
					my $flanks_right= "null";
					my $flank_exon_start="$genesite[2*$f]";
					my $flank_exon_end="$genesite[2*$f+1]";
					my $splice_site = 1;
					my $structure = "";
					my $splice_chain = "";
					my $splice_position = "";
					my $alt_exon_end = "";
					for(my $m=$f+1;2*$m<=$#genesite;$m++){
						if($genesite[2*$m]<=$position){
							$structure = "$splice_site\-\,".$structure.($splice_site+1)."\-";
							$splice_chain = "$position\-\,".$splice_chain.$genesite[2*$m]."\-";
							$splice_position = $splice_position.",$position,$genesite[2*$m]";
							$alt_exon_end = $genesite[2*$m+1];
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$genesite[2*$m]\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;

						}elsif($genesite[2*$m+1]<=$position){
							$structure .= "$splice_site\-";
							$structure = ($splice_site+1)."\-\,".$structure;
							$splice_chain .= "$genesite[2*$m]\-";
							$splice_chain = $position."\-\,".$splice_chain;
							$splice_position = $splice_position.",$genesite[2*$m],$position";
							$alt_exon_end = $genesite[2*$m+1];
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$position\t\.\t$strand\t$reads\t$info\n";
							next LABEL1;
						}
						if(2*$m+1==$#genesite){
							$structure .= "$splice_site\-".($splice_site+1)."\]";
							$splice_chain .= "$genesite[2*$m]\-$genesite[2*$m+1]\]";
							$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
						}else{
							$structure .= "$splice_site\-".($splice_site+1)."\^";
							$splice_chain .= "$genesite[2*$m]\-$genesite[2*$m+1]\^";
							$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
						}
						
						$splice_site = $splice_site+2;
					}
					$structure = "$splice_site\-\,".$structure;
					$splice_chain = "$position\-\,".$splice_chain;
					$splice_position = $splice_position.",$position";
					$alt_exon_end = $genesite[-1];
					$splice_position=~s/^,//;
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_A3SS "$chr\t$gene\tA3SS\t$flanks_left\t$position\t\.\t$strand\t$reads\t$info\n";
					next LABEL1;
				}
			}
		}
		
	}
}


###子过程-----find_mxe($chr,$strand,$gene,\@thisgene)
sub find_afe_ale {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	return 0 if @genesite < 4;
	if($strand eq "+"){
		if ( defined( $positive_end{$chr}->{$genesite[2]} ) ){
			my @end_junction = split /\s+/, $positive_end{$chr}->{$genesite[2]};
			foreach my $j (@end_junction) {
				my ( $position, $reads ) = split /\:/, $j;
				if($position<$genesite[0]){
					#flanks "null,831045-"; structure "1[2^,3[4^"; splice_chain "830120[830204^,830375[830922^";AFE
					my $flanks_left = "null";
					my $flanks_right= $genesite[2];
					my $flank_exon_start="$genesite[2]";
					my $flank_exon_end="$genesite[3]";
					my $alt_exon_start="unknown,$genesite[0]";
					my $alt_exon_end="$position,$genesite[1]";
					my $structure = "1\[2\^\,3\[4\^";
					$positive_sj{$chr}->{"$genesite[1]:$genesite[2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[1]:$genesite[2]"});
					my $model_sj_reads = $positive_sj{$chr}->{"$genesite[1]:$genesite[2]"};
					my $splice_chain = "unknown\[$position\^,$genesite[0]\[$genesite[1]\^";
					my $splice_position = "unknown,$position,$genesite[0],$genesite[1]";
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_AFE "$chr\t$gene\tAFE\t$position\t$flanks_right\t$model_sj_reads\t$strand\t$reads\t$info\n";
					my $junc = "$position:$genesite[2]";
					$usedjunc{$chr}->{"+"}->{$junc} = 1;
				}
			}
		}
		if ( defined( $positive_start{$chr}->{$genesite[-3]} ) ){
			my @start_junction = split /\s+/, $positive_start{$chr}->{$genesite[-3]};
			foreach my $j (@start_junction) {
				my ( $position, $reads ) = split /\:/, $j;
				if($position>$genesite[-1]){
					#flanks "831045^,null"; structure "3-4],1-2]"; splice_chain "830120-830204],830375-830922]";ALE
					my $flanks_left = $genesite[-3];
					my $flanks_right= "null";
					my $flank_exon_start="$genesite[-4]";
					my $flank_exon_end="$genesite[-3]";
					my $alt_exon_start="$genesite[-2],$position";
					my $alt_exon_end="$genesite[-1],unknown";
					my $structure = "3\-4\]\,1\-2\]";
					$positive_sj{$chr}->{"$genesite[-3]:$genesite[-2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[-3]:$genesite[-2]"});
					my $model_sj_reads = $positive_sj{$chr}->{"$genesite[-3]:$genesite[-2]"};
					my $splice_chain = "$position\-unknown\],$genesite[-2]\-$genesite[-1]\]";
					my $splice_position = "$genesite[-2],$genesite[-1],$position,unknown";
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_ALE "$chr\t$gene\tALE\t$flanks_left\t$position\t$model_sj_reads\t$strand\t$reads\t$info\n";
					my $junc = "$genesite[-3]:$position";
					$usedjunc{$chr}->{"+"}->{$junc} = 1;
				}
			}
		}
	}
	if($strand eq "-"){
		if ( defined( $negative_end{$chr}->{$genesite[2]} ) ){
			my @end_junction = split /\s+/, $negative_end{$chr}->{$genesite[2]};
			foreach my $j (@end_junction) {
				my ( $position, $reads ) = split /\:/, $j;
				if($position>$genesite[0]){
					#flanks "null,831045-"; structure "1[2^,3[4^"; splice_chain "830120[830204^,830375[830922^";AFE
					my $flanks_left = "null";
					my $flanks_right= $genesite[2];
					my $flank_exon_start="$genesite[2]";
					my $flank_exon_end="$genesite[3]";
					my $alt_exon_start="unknown,$genesite[0]";
					my $alt_exon_end="$position,$genesite[1]";
					my $structure = "1\[2\^\,3\[4\^";
					$negative_sj{$chr}->{"$genesite[1]:$genesite[2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[1]:$genesite[2]"});
					my $model_sj_reads = $negative_sj{$chr}->{"$genesite[1]:$genesite[2]"};
					my $splice_chain = "unknown\[$position\^,$genesite[0]\[$genesite[1]\^";
					my $splice_position = "unknown,$position,$genesite[0],$genesite[1]";
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_AFE "$chr\t$gene\tAFE\t$position\t$flanks_right\t$model_sj_reads\t$strand\t$reads\t$info\n";
					my $junc = "$position:$genesite[2]";
					$usedjunc{$chr}->{"-"}->{$junc} = 1;
				}
			}
		}
		if ( defined( $negative_start{$chr}->{$genesite[-3]} ) ){
			my @start_junction = split /\s+/, $negative_start{$chr}->{$genesite[-3]};
			foreach my $j (@start_junction) {
				my ( $position, $reads ) = split /\:/, $j;
				if($position<$genesite[-1]){
					#flanks "831045^,null"; structure "3-4],1-2]"; splice_chain "830120-830204],830375-830922]";ALE
					my $flanks_left = $genesite[-3];
					my $flanks_right= "null";
					my $flank_exon_start="$genesite[-4]";
					my $flank_exon_end="$genesite[-3]";
					my $alt_exon_start="$genesite[-2],$position";
					my $alt_exon_end="$genesite[-1],unknown";
					my $structure = "3\-4\]\,1\-2\]";
					$negative_sj{$chr}->{"$genesite[-3]:$genesite[-2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[-3]:$genesite[-2]"});
					my $model_sj_reads = $negative_sj{$chr}->{"$genesite[-3]:$genesite[-2]"};
					my $splice_chain = "$position\-unknown\],$genesite[-2]\-$genesite[-1]\]";
					my $splice_position = "$genesite[-2],$genesite[-1],$position,unknown";
					my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; alt_exon_start \"$alt_exon_start\"\; alt_exon_end \"$alt_exon_end\"\; splice_position \"$splice_position\"\;";
					print OUT_ALE "$chr\t$gene\tALE\t$flanks_left\t$position\t$model_sj_reads\t$strand\t$reads\t$info\n";
					my $junc = "$genesite[-3]:$position";
					$usedjunc{$chr}->{"-"}->{$junc} = 1;
				}
			}
		}
	}
}

###子过程-----find_mxe($chr,$strand,$gene,\@thisgene)
sub find_mxe {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f+1<$#genesite-2;$f++){
			my $junc_start = $genesite[2*$f+1];
			my $junc_end   = $genesite[2*$f+4];
			if ( defined( $positive_start{$chr}->{$junc_start} ) && defined( $positive_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $positive_start{$chr}->{$junc_start};
				my @end_junction = split /\s+/, $positive_end{$chr}->{$junc_end};
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					next if $position_left>= $junc_end; #junction的终点越过了下一个exon的起点
					foreach my $j (@end_junction) {
						my ( $position_right, $reads_right ) = split /\:/, $j;
						next if $position_right<=$junc_start; #junction的起点越过了上一个exon的终点
						next if $position_right<=$position_left; #两个junction有交叉
						if($gffsj==1){
							next if not defined($exon_hash{$chr}->{$strand}->{"$position_left:$position_right"});
						}else{
							next if $position_right-$position_left>400; #cassette exon的长度不应太长，这里取400（人类的exon平均长度为200bp左右）；
						}
						if($position_left>$junc_start && $position_left<$genesite[2*$f+2] && $position_right>$junc_start && $position_right<$genesite[2*$f+2]){
							my $flanks_left = $junc_start;
							my $flanks_right= $junc_end;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+4]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+5]";
							my $mxe_start="$position_left,$genesite[2*$f+2]";
							my $mxe_end="$position_right,$genesite[2*$f+3]";
							my $structure = "1\-2\^,3\-4\^";
							my $splice_chain = "$position_left\-$position_right\^,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
							my $splice_position = "$position_left,$position_right,$genesite[2*$f+2],$genesite[2*$f+3]";
							my $flanks_exon_no = "$f,".($f+2);
							my $type1_sj_reads = $reads_left.",".$reads_right;
							$positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							$positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $type2_sj_reads = ($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"}).",".($positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; mxe_start \"$mxe_start\"\; mxe_end \"$mxe_end\"\; splice_position \"$splice_position\"\;";
							print OUT_MXE "$chr\t$gene\tMXE\t$flanks_left\t$flanks_right\t$type1_sj_reads\t$strand\t$type2_sj_reads\t$info\n";
							my $junc = "$junc_start:$position_left";
							$usedjunc{$chr}->{"+"}->{$junc} = 1;
							$junc = "$position_right:$junc_end";
							$usedjunc{$chr}->{"+"}->{$junc} = 1;
						}elsif($position_left>$genesite[2*$f+3] && $position_left<$junc_end && $position_right>$genesite[2*$f+3] && $position_right<$junc_end){
							my $flanks_left = $junc_start;
							my $flanks_right= $junc_end;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+4]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+5]";
							my $mxe_start="$genesite[2*$f+2],$position_left";
							my $mxe_end="$genesite[2*$f+3],$position_right";
							my $structure = "3\-4\^,1\-2\^";
							my $splice_chain = "$position_left\-$position_right\^,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
							my $splice_position = "$genesite[2*$f+2],$genesite[2*$f+3],$position_left,$position_right";
							my $flanks_exon_no = "$f,".($f+2);
							my $type2_sj_reads = $reads_left.",".$reads_right;
							$positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							$positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $type1_sj_reads = ($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"}).",".($positive_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; mxe_start \"$mxe_start\"\; mxe_end \"$mxe_end\"\; splice_position \"$splice_position\"\;";
							print OUT_MXE "$chr\t$gene\tMXE\t$flanks_left\t$flanks_right\t$type1_sj_reads\t$strand\t$type2_sj_reads\t$info\n";
							my $junc = "$junc_start:$position_left";
							$usedjunc{$chr}->{"+"}->{$junc} = 1;
							$junc = "$position_right:$junc_end";
							$usedjunc{$chr}->{"+"}->{$junc} = 1;
						}
					}
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f+1<$#genesite-2;$f++){
			my $junc_start = $genesite[2*$f+1];
			my $junc_end   = $genesite[2*$f+4];
			if ( defined( $negative_start{$chr}->{$junc_start} ) && defined( $negative_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $negative_start{$chr}->{$junc_start};
				my @end_junction = split /\s+/, $negative_end{$chr}->{$junc_end};
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					next if $position_left<= $junc_end; #junction的终点越过了下一个exon的起点
					foreach my $j (@end_junction) {
						my ( $position_right, $reads_right ) = split /\:/, $j;
						next if $position_right>=$junc_start; #junction的起点越过了上一个exon的终点
						next if $position_right>=$position_left; #两个junction有交叉
						if($gffsj==1){
							next if not defined($exon_hash{$chr}->{$strand}->{"$position_left:$position_right"});
						}else{
							next if $position_left-$position_right>400; #cassette exon的长度不应太长，这里取400（人类的exon平均长度为200bp左右）；
						}
						if($position_left<$junc_start && $position_left>$genesite[2*$f+2] && $position_right<$junc_start && $position_right>$genesite[2*$f+2]){
							my $flanks_left = $junc_start;
							my $flanks_right= $junc_end;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+4]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+5]";
							my $mxe_start="$position_left,$genesite[2*$f+2]";
							my $mxe_end="$position_right,$genesite[2*$f+3]";
							my $structure = "1\-2\^,3\-4\^";
							my $splice_chain = "$position_left\-$position_right\^,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
							my $splice_position = "$position_left,$position_right,$genesite[2*$f+2],$genesite[2*$f+3]";
							my $flanks_exon_no = "$f,".($f+2);
							my $type1_sj_reads = $reads_left.",".$reads_right;
							$negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							$negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $type2_sj_reads = ($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"}).",".($negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; mxe_start \"$mxe_start\"\; mxe_end \"$mxe_end\"\; splice_position \"$splice_position\"\;";
							print OUT_MXE "$chr\t$gene\tMXE\t$flanks_left\t$flanks_right\t$type1_sj_reads\t$strand\t$type2_sj_reads\t$info\n";
							my $junc = "$junc_start:$position_left";
							$usedjunc{$chr}->{"-"}->{$junc} = 1;
							$junc = "$position_right:$junc_end";
							$usedjunc{$chr}->{"-"}->{$junc} = 1;
						}elsif($position_left<$genesite[2*$f+3] && $position_left>$junc_end && $position_right<$genesite[2*$f+3] && $position_right>$junc_end){
							my $flanks_left = $junc_start;
							my $flanks_right= $junc_end;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+4]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+5]";
							my $mxe_start="$genesite[2*$f+2],$position_left";
							my $mxe_end="$genesite[2*$f+3],$position_right";
							my $structure = "3\-4\^,1\-2\^";
							my $splice_chain = "$position_left\-$position_right\^,$genesite[2*$f+2]\-$genesite[2*$f+3]\^";
							my $splice_position = "$genesite[2*$f+2],$genesite[2*$f+3],$position_left,$position_right";
							my $flanks_exon_no = "$f,".($f+2);
							my $type2_sj_reads = $reads_left.",".$reads_right;
							$negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							$negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $type1_sj_reads = ($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"}).",".($negative_sj{$chr}->{"$genesite[2*$f+3]:$genesite[2*$f+4]"});
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; mxe_start \"$mxe_start\"\; mxe_end \"$mxe_end\"\; splice_position \"$splice_position\"\;";
							print OUT_MXE "$chr\t$gene\tMXE\t$flanks_left\t$flanks_right\t$type1_sj_reads\t$strand\t$type2_sj_reads\t$info\n";
							my $junc = "$junc_start:$position_left";
							$usedjunc{$chr}->{"-"}->{$junc} = 1;
							$junc = "$position_right:$junc_end";
							$usedjunc{$chr}->{"-"}->{$junc} = 1;
						}
					}
				}
			}
		}

	}
}


###子过程-----find_ce($chr,$strand,$gene,\@thisgene)
sub find_ce {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	if($strand eq "+"){
		for(my $f=0;2*$f+1<$#genesite;$f++){   #从第一个exon到倒数第二个exon
			my $junc_start = $genesite[2*$f+1];
			my $junc_end   = $genesite[2*$f+2];
			if ( defined( $positive_start{$chr}->{$junc_start} ) && defined( $positive_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $positive_start{$chr}->{$junc_start};
				my @end_junction = split /\s+/, $positive_end{$chr}->{$junc_end};
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					next if $position_left>= $junc_end; #junction的终点越过了下一个exon的起点
					foreach my $j (@end_junction) {
						my ( $position_right, $reads_right ) = split /\:/, $j;
						next if $position_right<=$junc_start; #junction的起点越过了上一个exon的终点
						next if $position_right<=$position_left; #两个junction有交叉
						next if $position_right-$position_left>400; #cassette exon的长度不应太长，这里取400（人类的exon平均长度为200bp左右）；
						# my $skipped_exon_len = $position_right-$position_left;
						my $flanks_left = $junc_start;
						my $flanks_right= $junc_end;
						my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+2]";
						my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+3]";
						my $skipped_exon_start="$position_left";
						my $skipped_exon_end="$position_right";
						my $structure = "1\-2\^,0";
						my $splice_chain = "$position_left\-$position_right\^,";
						my $splice_position = "$position_left,$position_right";
						my $flanks_exon_no = "$f,".($f+1);
						$positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
						my $exclusive_sj_reads = $positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"};
						my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; skipped_exon_start \"$skipped_exon_start\"\; skipped_exon_end \"$skipped_exon_end\"\; splice_position \"$splice_position\"\;";
						print OUT_CE "$chr\t$gene\tCE\t$flanks_left\t$flanks_right\t$reads_left,$reads_right\t$strand\t$exclusive_sj_reads\t$info\n";
						my $junc = "$junc_start:$position_left";
						$usedjunc{$chr}->{"+"}->{$junc} = 1;
						$junc = "$position_right:$junc_end";
						$usedjunc{$chr}->{"+"}->{$junc} = 1;
					}
				}
			}
		}
	}
	if($strand eq "-"){
		for(my $f=0;2*$f+1<$#genesite;$f++){
			my $junc_start = $genesite[2*$f+1];
			my $junc_end   = $genesite[2*$f+2];
			if ( defined( $negative_start{$chr}->{$junc_start} ) && defined( $negative_end{$chr}->{$junc_end} ) ){
				my @start_junction = split /\s+/, $negative_start{$chr}->{$junc_start};
				my @end_junction = split /\s+/, $negative_end{$chr}->{$junc_end};
				foreach my $i (@start_junction) {
					my ( $position_left, $reads_left ) = split /\:/, $i;
					next if $position_left<= $junc_end; #junction的终点越过了下一个exon的起点
					foreach my $j (@end_junction) {
						my ( $position_right, $reads_right ) = split /\:/, $j;
						next if $position_right>=$junc_start; #junction的起点越过了上一个exon的终点
						next if $position_right>=$position_left; #两个junction有交叉
						next if $position_left-$position_right>400; #cassette exon的长度不应太长，这里取400（人类的exon平均长度为200bp左右）；
						# my $skipped_exon_len = $position_left-$position_right;
						my $flanks_left = $junc_start;
						my $flanks_right= $junc_end;
						my $flank_exon_start="$genesite[2*$f],$genesite[2*$f+2]";
						my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$f+3]";
						my $skipped_exon_start="$position_left";
						my $skipped_exon_end="$position_right";
						my $structure = "1\-2\^,0";
						my $splice_chain = "$position_left\-$position_right\^,";
						my $splice_position = "$position_left,$position_right";
						my $flanks_exon_no = "$f,".($f+1);
						$negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
						my $exclusive_sj_reads = $negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"};
						my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; flanks_exon_no \"$flanks_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; skipped_exon_start \"$skipped_exon_start\"\; skipped_exon_end \"$skipped_exon_end\"\; splice_position \"$splice_position\"\;";
						print OUT_CE "$chr\t$gene\tCE\t$flanks_left\t$flanks_right\t$reads_left,$reads_right\t$strand\t$exclusive_sj_reads\t$info\n";
						my $junc = "$junc_start:$position_left";
						$usedjunc{$chr}->{"-"}->{$junc} = 1;
						$junc = "$position_right:$junc_end";
						$usedjunc{$chr}->{"-"}->{$junc} = 1;
					}
				}
			}
		}
	}
}


###子过程-----find_se($chr,$strand,$gene,\@thisgene)
sub find_se {
	my ($chr,$strand,$gene,$thisgene) = @_;
	my @genesite = @{$thisgene};
	#正向基因通过正向的junction判断，从左到右看每个以“^”site为起点的junction的终点是不是属于“-”site。
	if($strand eq "+"){
		for(my $f=0;2*$f+3<$#genesite;$f++){ #从第一个exon到倒数第三个exon
			my $junc_start = $genesite[2*$f+1];
			if ( defined( $positive_start{$chr}->{$junc_start} ) ){
				my @start_junction = split /\s+/, $positive_start{$chr}->{$junc_start};
				#查看以该位点为起点的每个junction的end
				foreach my $i (@start_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					#查看该junction的end是否与后面的某个"-"site相同
					for(my $r=$f+2;2*$r<$#genesite;$r++){
						if($position==$genesite[2*$r]){
							my $flanks_left = $junc_start;
							my $flanks_right= $position;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$r]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$r+1]";
							my $skipped_exon_start="";
							my $skipped_exon_end="";
							my $splice_site = 1;
							my $structure = "0,";
							my $splice_chain = ",";
							my $skipped_exon_no = "";
							my $splice_position = "";
							$positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							my $inclusive_sj_reads = $positive_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"};
							for(my $m=$f+1;$m<$r;$m++){
								$positive_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"} = 0 if not defined($positive_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"});
								$inclusive_sj_reads = $inclusive_sj_reads.",".($positive_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"});
								if($m==$f+1){
									$skipped_exon_no = $skipped_exon_no.$m;
									$skipped_exon_start="$genesite[2*$m]";
									$skipped_exon_end="$genesite[2*$m+1]";
								}else{
									$skipped_exon_no = $skipped_exon_no.",".$m;
									$skipped_exon_start=$skipped_exon_start.","."$genesite[2*$m]";
									$skipped_exon_end=$skipped_exon_end.","."$genesite[2*$m+1]";
								}
								$structure=$structure.$splice_site."-".($splice_site+1)."\^";
								$splice_chain = $splice_chain.$genesite[2*$m]."-".$genesite[2*$m+1]."\^";
								$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
								my $junc = "$genesite[2*$m-1]:$genesite[2*$m]";
								$usedjunc{$chr}->{"+"}->{$junc} = 1;
								$splice_site+=2;
							}
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; skipped_exon_no \"$skipped_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; skipped_exon_start \"$skipped_exon_start\"\; skipped_exon_end \"$skipped_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_SE "$chr\t$gene\tSE\t$flanks_left\t$flanks_right\t$inclusive_sj_reads\t$strand\t$reads\t$info\n";
							my $junc = "$junc_start:$position";
							$usedjunc{$chr}->{"+"}->{$junc} = 1;
						}
					}
				}
			}
		}
	}
	#负向基因通过负向的junction判断，从右到左看每个以“^”site为起点的junction的终点是不是属于“-”site。
	if($strand eq "-"){
		for(my $f=0;2*$f+3<$#genesite;$f++){
			my $junc_start = $genesite[2*$f+1];
			if ( defined( $negative_start{$chr}->{$junc_start} ) ){
				my @start_junction = split /\s+/, $negative_start{$chr}->{$junc_start};
				#查看以该位点为起点的每个junction的end
				foreach my $i (@start_junction) {
					my ( $position, $reads ) = split /\:/, $i;
					#查看该junction的end是否与后面的某个"-"site相同
					for(my $r=$f+2;2*$r<$#genesite;$r++){
						if($position==$genesite[2*$r]){
							my $flanks_left = $junc_start;
							my $flanks_right= $position;
							my $flank_exon_start="$genesite[2*$f],$genesite[2*$r]";
							my $flank_exon_end="$genesite[2*$f+1],$genesite[2*$r+1]";
							my $skipped_exon_start="";
							my $skipped_exon_end="";
							my $splice_site = 1;
							my $structure = "0,";
							my $splice_chain = ",";
							my $skipped_exon_no = "";
							my $splice_position = "";
							$negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"});
							my $inclusive_sj_reads = $negative_sj{$chr}->{"$genesite[2*$f+1]:$genesite[2*$f+2]"};
							for(my $m=$f+1;$m<$r;$m++){
								$negative_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"} = 0 if not defined($negative_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"});
								$inclusive_sj_reads = $inclusive_sj_reads.",".($negative_sj{$chr}->{"$genesite[2*$m+1]:$genesite[2*$m+2]"});
								if($m==$f+1){
									$skipped_exon_no = $skipped_exon_no.$m;
									$skipped_exon_start="$genesite[2*$m]";
									$skipped_exon_end="$genesite[2*$m+1]";
								}else{
									$skipped_exon_no = $skipped_exon_no.",".$m;
									$skipped_exon_start=$skipped_exon_start.","."$genesite[2*$m]";
									$skipped_exon_end=$skipped_exon_end.","."$genesite[2*$m+1]";
								}
								$structure=$structure.$splice_site."-".($splice_site+1)."\^";
								$splice_chain = $splice_chain.$genesite[2*$m]."-".$genesite[2*$m+1]."\^";
								$splice_position = $splice_position.",$genesite[2*$m],$genesite[2*$m+1]";
								my $junc = "$genesite[2*$m-1]:$genesite[2*$m]";
								$usedjunc{$chr}->{"-"}->{$junc} = 1;
								$splice_site+=2;
							}
							$splice_position=~s/^,//;
							my $info = "flanks \"$flanks_left\^\,$flanks_right\-\"\; structure \"$structure\"\; splice_chain \"$splice_chain\"\; skipped_exon_no \"$skipped_exon_no\"\; flank_exon_start \"$flank_exon_start\"\; flank_exon_end \"$flank_exon_end\"\; skipped_exon_start \"$skipped_exon_start\"\; skipped_exon_end \"$skipped_exon_end\"\; splice_position \"$splice_position\"\;";
							print OUT_SE "$chr\t$gene\tSE\t$flanks_left\t$flanks_right\t$inclusive_sj_reads\t$strand\t$reads\t$info\n";
							my $junc = "$junc_start:$position";
							$usedjunc{$chr}->{"-"}->{$junc} = 1;
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
		my ( $chr, $start, $end, $strand) = split /\t/ ;
		if ( $strand eq "+" ) {
			$positive_sj{$chr}->{"$start:$end"} = 3;
			if ( !defined( $positive_start{$chr}->{$start} ) ) {
				$positive_start{$chr}->{$start} = "$end\:3";
			}
			else {
				$positive_start{$chr}->{$start} .= "\t$end\:3";

			}
			if ( !defined( $positive_end{$chr}->{$end} ) ) {
				$positive_end{$chr}->{$end} = "$start\:3";
			}
			else {
				$positive_end{$chr}->{$end} .= "\t$start\:3";

				#调试信息
				#print $positive_end{$chr}->{$end},"\n";
			}
		}
		else {
			$negative_sj{$chr}->{"$start:$end"} = 3;
			if ( !defined( $negative_start{$chr}->{$start} ) ) {
				$negative_start{$chr}->{$start} = "$end\:3";
			}
			else {
				$negative_start{$chr}->{$start} .= "\t$end\:3";
			}
			if ( !defined( $negative_end{$chr}->{$end} ) ) {
				$negative_end{$chr}->{$end} = "$start\:3";
			}
			else {
				$negative_end{$chr}->{$end} .= "\t$start\:3";
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
			$gene = "";
			next;
		}
		if ( $feature =~ /^(mRNA|miRNA|mRNA_TE_gene|ncRNA|rRNA|snoRNA|snRNA|tRNA|transcript)$/ ) {
			#将上一个gene的各位点排序，正向的按从小到大，负向的按从大到小。
			if($gene ne "" and defined(@{$gff->{$this_chr}->{$this_strand}->{$gene}})){
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