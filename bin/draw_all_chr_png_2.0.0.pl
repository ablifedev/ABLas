#!/usr/bin/perl

# Copyright (c)   AB_Life 2012
# Writer:         XUXIONG
# Program Date:   2011.
# Modifier:       CHENGCHAO
# Last Modified:  2012.

use strict;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use GD::Simple;
use warnings;

#Before writing your programmeyou must write the detailed time /discriptions /parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

if (scalar(@ARGV)<4) {
	print "  When you run the program ,the sam|bam file needs to be sorted,the bam file can use the samtools(sort) to sort,the commadline:\"samtools sort infile.bam infile_sort &\",if the file format is sam,you can use the commmandline:\"sort -k 3,3 -k 4,4n infile.sam > infile_sort.sam &\" \n";
	print "Example1:\n  perl $0 /data2/verticillium_dahliae_vdls/gene/verticillium_dahliae_vdls.17_1_transcripts.gff /data2/human/seq_20120222_filter_reads/verticillium_dahliae_vdls/RNA_0224/sequence_3_clean.fq_0224/accepted_hits.sam 0228 5000 0 1\n";
	print "Example2:\n  perl $0 /data2/human/hg18.gff /data1/test/xux/test_software/human_1026/SpliceMap/A3_L1_1.fq/good_hits.sam 0227 5000 5000 1 /data2/human/kgXref.txt\n";
	print "Example3:\n  perl $0 /data2/arabidopsis/TAIR10/GFF/TAIR10_GFF3_genes.gff /data1/test/xux/arabidopsis_seq_20120614_seq_20120618/arabidopsis_0711/Ninanjie_mRNA_W_small_clean.fq_0710/accepted_hits.sam.uniq  Ninanjie_mRNA_W_small 10 10 0 /data2/arabidopsis/TAIR10/function_annotation/gene_aliases.20101027\n";
	print "\nUsage:\n  perl $0 gff sam|bam outfile promoter_length downstream_length merge_dentical_reads[0|1] kgXref ......\n";
	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################
my $current_dir = `pwd`;chomp($current_dir);

my $gff = $ARGV[0];
my $infile = $ARGV[1];
my $outfile = $ARGV[2];
my $promoter=$ARGV[3] ;
my $downstream=$ARGV[4] ;
my $merge_dentical_reads= $ARGV[5];
my $kgXrefinfile = $ARGV[6] if (defined $ARGV[6] ) ;
my $asfile = $ARGV[7] if (defined $ARGV[7] ) ;

my %kgXref=();
&load_kgXref($kgXrefinfile,\%kgXref) if (defined $ARGV[6] ) ;
my %gene=();
&load_gff($gff,\%gene,\%kgXref);
#print Dumper %gene;
my %cluster_gene=();
&cluster_gene(\%gene,\%cluster_gene,$promoter,$downstream);
&output_clustered_gene(\%cluster_gene,$outfile);
##############test#####
#undef(%gene);
#my @array_test;
#for (my $i=1;$i<=50;$i++) {
#&draw(\@array_test,$cluster_gene{"Chr1"}->[$i],$cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[6] eq "+" ? $cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[3]-$promoter : $cluster_gene{"Chr1"}->[$i]->[0]->[0]->[0]->[3],max(map {$_->[0]->[0]->[4]} @{$cluster_gene{"Chr1"}->[$i]}),"Chr1_cluster$i.png")
#}
#exit;
#######################
my %as = ();
my $direction_bool=0;
&load_sorted_as($asfile,\%cluster_gene,\%as,$promoter,$downstream) if (defined $ARGV[7] );


my %distribution=();

chdir($current_dir);
if ($infile=~/\.sam/) {
	&load_sorted_sam($infile,$outfile,\%cluster_gene,\%distribution,\%as,$direction_bool,$promoter,$downstream);
}
elsif($infile=~/\.bam$/){
	$infile="samtools view $infile|";
	&load_sorted_sam($infile,$outfile,\%cluster_gene,\%distribution,\%as,$direction_bool,$promoter,$downstream);
}
#elsif($infile=~/\.bed$/){
#	
#}
#elsif(){
#	
#}

sub load_kgXref{
	my ($infile,$hash_ref)=@_;
	open (IN ,"<$infile") || die "Can't open $!";
	while (<IN>) {
		chomp;
		my @line=split(/\t+/,$_);
		my $geneID=$line[0];
		if ($geneID=~/([\w\.]+)/) {
			$hash_ref->{$1}=join("\t",@line);
		}
		else{
			$hash_ref->{$geneID}=join("\t",@line);
		}
	}
	close IN ;
	print "Done load kgXref\n";
}

sub load_gff{
	my ($GFF_INFILE,$hash,$kgXref)=@_;
	my @gene_content=();
	my $RNA_index=0;
	my $gene_start=0;
	my $CHROMOSOME="";
	open (IN,"$GFF_INFILE") || die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/) ;
		my @line=split(/\s+/,$_);
		if ($line[2]=~/^gene|pseudogene|transposable_element_gene$/ ) {
			if (@gene_content){
				foreach my $mRNA (@gene_content[1..$#gene_content]) {
					my $mRNAID="";
					if ($mRNA->[0]->[-1]=~/^ID=([\w\.]+)\;/) {
						$mRNAID=$1 ;
						push @{$mRNA->[0]},$kgXref->{$mRNAID} if (exists($kgXref->{$mRNAID})) ;
					}
					@{$mRNA}[1..$#{$mRNA}]=sort {$a->[3]<=>$b->[3]} sort {$a->[2] cmp $b->[2]}  @{$mRNA}[1..$#{$mRNA}];
				}
#				&build_index(\@gene_content,$promoter,10000,$hash,$CHROMOSOME);
				push @{$hash->{$CHROMOSOME}},[@gene_content];
			}
			$RNA_index=0;
			($CHROMOSOME,$gene_start)=@line[0,3];
			@gene_content=();
			push @{$gene_content[$RNA_index]},[@line];
		}
		elsif($line[2]=~/^mRNA|miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/){
			$RNA_index++;
			push @{$gene_content[$RNA_index]},[@line];
		}
		elsif($line[2]=~/^chromosome|region$/i){
			next;
		}
		else{
			push @{$gene_content[$RNA_index]},[@line];
		}
		if (eof(IN)) {
			if (@gene_content){
				foreach my $mRNA (@gene_content[1..$#gene_content]) {
					@{$mRNA}[1..$#{$mRNA}]=sort {$a->[3]<=>$b->[3]} sort {$a->[2] cmp $b->[2]}  @{$mRNA}[1..$#{$mRNA}];
				}
#				&build_index(\@gene_content,$promoter,10000,$hash,$CHROMOSOME);
				push @{$hash->{$CHROMOSOME}},[@gene_content];
			}
		}
	}
	close(IN);
	&sort_gene($hash);
	print "Done load gff\n";
}

sub sort_gene{
	my $hash_ref=shift;
	foreach my $TEMP_CHRO (keys %{$hash_ref}) {
		@{$hash_ref->{$TEMP_CHRO}}= sort {$a->[0]->[0]->[3]<=>$b->[0]->[0]->[3] or $a->[0]->[0]->[4]<=>$b->[0]->[0]->[4] } @{$hash_ref->{$TEMP_CHRO}};
	}
}

sub build_index{
	my ($GENE_CONTENT,$PROMOTER_LENGTH,$SIZE,$HASH,$CHROM)=@_;
	my ($blockNum_start,$blockNum_end);
	if($GENE_CONTENT->[0]->[0]->[6]eq"+"){
		$blockNum_start=int(($GENE_CONTENT->[0]->[0]->[3]-$PROMOTER_LENGTH)/$SIZE);
		$blockNum_end=int(($GENE_CONTENT->[0]->[0]->[4])/$SIZE);
	}
	else{
		$blockNum_start=int($GENE_CONTENT->[0]->[0]->[3]/$SIZE);
		$blockNum_end=int(($GENE_CONTENT->[0]->[0]->[4]+$PROMOTER_LENGTH)/$SIZE);
	}
	map {push(@{$HASH->{$CHROM}{$_}},[@{$GENE_CONTENT}])}($blockNum_start..$blockNum_end);
}

sub cluster_gene{
	my ($gene_ref,$gene_cluster_ref,$PROMOTER_LENGTH,$DOWNSTREAM_LENGTH)=@_;
	foreach my $TEMP_CHRO (sort keys %{$gene_ref}) {
		my $cluster_index=0;
		my ($Start,$End,$current_Start,$current_End)=(0,0,0,0);
		foreach my $temp_RNA (@{$gene_ref->{$TEMP_CHRO}}) {
			if ($temp_RNA->[0]->[0]->[6] eq "+") {
				($Start,$End)=($temp_RNA->[0]->[0]->[3]-$PROMOTER_LENGTH,$temp_RNA->[0]->[0]->[4]+$DOWNSTREAM_LENGTH);
				$Start=1 if ($Start<0) ;
			}else{
				($Start,$End)=($temp_RNA->[0]->[0]->[3]-$DOWNSTREAM_LENGTH,$temp_RNA->[0]->[0]->[4]+$PROMOTER_LENGTH);
				$Start=1 if ($Start<0) ;
			}
			my @current_region=&region_combination($Start,$End,$current_Start,$current_End);
			if (@current_region) {
				($current_Start,$current_End)=@current_region;
			}else{
				($current_Start,$current_End)=($Start,$End);
				$cluster_index++;
			}
			push @{$gene_cluster_ref->{$TEMP_CHRO}->[$cluster_index-1]},[@{$temp_RNA}];
		}
	}
	print "Done cluster gene\n";
}

sub output_clustered_gene{
	my ($cluster_ref,$OUTFILE)=@_;
	open(OUT,">".$OUTFILE."cluster.gff") || die $!;
	foreach my $TEMP_CHRO (sort keys %{$cluster_ref}) {
		for(my $cluster_index=0;$cluster_index<@{$cluster_ref->{$TEMP_CHRO}};$cluster_index++) {
			print OUT ">$TEMP_CHRO\_cluster$cluster_index\n";
			foreach my $gene (@{$cluster_ref->{$TEMP_CHRO}->[$cluster_index]}) {
				foreach my $RNA (@{$gene}) {
					foreach my $temp_line (@{$RNA}) {
						print OUT join("\t",@{$temp_line}),"\n";
					}
				}
			}
		}
	}
	close OUT;
	print "Done output cluster\n";
}

sub region_combination{
	my ($start1,$end1,$start2,$end2)=@_;
	if ($end2<$start1) {
		return;
	}else{
		return(min($start1,$start2),max($end1,$end2));
	}
}

sub load_sorted_as {
	my ($AS_INFILE,$CLUSTER_GENE,$ashash,$promoter_length,$downstream_length)=@_;
	my @as_content=();
	my $as_cluster_index=0;
	my $gene_start=0;
	my $CHROMOSOME="";
	my $current_chr=(sort keys(%{$CLUSTER_GENE}))[0];
	my $cluster_index=0;
	my $cluster_start=min(map {$_->[0]->[0]->[3]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
	$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
	$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
	my $tag_chromosome_end=0;
	open (IN,"$AS_INFILE") || die $!;
	while (<IN>) {
		chomp;
		next if (/^\@|^\#/) ;
		my @line=split(/\t+/,$_);
		if ($line[0] ne $current_chr) {
			if ($tag_chromosome_end==0) {
				push @{$ashash->{$current_chr}->[$cluster_index]},[@as_content] if (scalar(@as_content)>0) ;
				# foreach my $lineds (@as_content){
				# 	print join("\t",@{$lineds}),"\n";	
				# }
			}

			$current_chr=$line[0];
			@as_content=();
			$cluster_index=0;
#			$cluster_end=(sort {$a->[0]->[0]->[4]<=>$b->[0]->[0]->[4]} map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]})[-1]->[0]->[0]->[4];
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[0]}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$line[0]}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
		}
		if ($line[1]>$cluster_end ) {
			if ($cluster_index+1>$#{$CLUSTER_GENE->{$current_chr}}) {
				if ($tag_chromosome_end==0) {
					push @{$ashash->{$current_chr}->[$cluster_index]},[@as_content] if (scalar(@as_content)>0) ;
				# 	foreach my $lineds (@as_content){
				# 	print join("\t",@{$lineds}),"\n";	
				# }	
				}
				$tag_chromosome_end=1;
				next;
			}else{
				$tag_chromosome_end=0;
			}
			push @{$ashash->{$current_chr}->[$cluster_index]},[@as_content] if (scalar(@as_content)>0) ;
			# foreach my $lineds (@as_content){
			# 		print join("\t",@{$lineds}),"\n";	
			# 	}
			@as_content=();
			$cluster_index++;
			
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			while ($line[1]>$cluster_end) {
				$cluster_index++;
				if ($cluster_index>$#{$CLUSTER_GENE->{$current_chr}}){
					$cluster_index--;
					$tag_chromosome_end=1;
					last;
				}
				$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
				$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
				$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
				$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			}
		}
		$current_chr=$line[0];
		push @as_content,[@line];
		
		if (eof(IN)) {
			if ($tag_chromosome_end==0) {
				push @{$ashash->{$current_chr}->[$cluster_index]},[@as_content] if (scalar(@as_content)>0) ;
				# foreach my $lineds (@as_content){
				# 	print join("\t",@{$lineds}),"\n";	
				# }
			}
			print "start as cluster output\n";
			foreach my $TEMP_CHRO (sort keys %{$ashash}) {
				for(my $cluster_index=0;$cluster_index<@{$ashash->{$TEMP_CHRO}};$cluster_index++) {
					print ">$TEMP_CHRO\_cluster$cluster_index\n";
					foreach my $gene (@{$ashash->{$TEMP_CHRO}->[$cluster_index]}) {
						foreach my $temp_line (@{$gene}) {
							print join("\t",@{$temp_line}),"\n";
						}
					}
				}
			}
		}

	}

}

sub load_sorted_sam	{
	my ($SAM_FILE,$OUTFILE,$CLUSTER_GENE,$distribution,$CLUSTER_AS,$directional,$promoter_length,$downstream_length)=@_;
	open (IN,$SAM_FILE) || die $!;
	my @gene_reads=();
	my @intergenic_reads=();
	my $current_chr=(sort keys(%{$CLUSTER_GENE}))[0];
	my %reads_distribution=();
	my $cluster_index=0;
	mkdir("$current_dir/$current_chr");
	chdir("$current_dir/$current_chr");
	my $cluster_start=min(map {$_->[0]->[0]->[3]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
	my $cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
	$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
	$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
	my $tag_chromosome_end=0;
	while (<IN>) {
		chomp;
		next if (/^\@|^\#/) ;
		my @line=split(/\t+/,$_);
		next if ($line[3]==0 or $line[2]=~/Zv/i) ;
		my $read_strand = $line[1] & (0x0010) ? "-" : "+";
#		next if ($line[5]!~/^[\dMNISP=XHD]+$/) ;
		my $read_end;
		if ($line[5]=~/^(\d+)M$/) {
			$read_end=$line[3]+$1-1;
		}
		elsif($line[5]=~/^(\d+)M(\d+)N(\d+)M$/){
			$read_end=$line[3]+$1+$2+$3-1;
		}
		else{
			next;
		}
		if ($line[2] ne $current_chr) {
			if ($tag_chromosome_end==0) {
				&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$CLUSTER_AS->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png") if (scalar(@gene_reads)>10) ;
			}

			mkdir("$current_dir/$line[2]");
			chdir("$current_dir/$line[2]");
			print "Done SAM $current_chr\n";
			$current_chr=$line[2];
			@gene_reads=();
			@intergenic_reads=();
			$cluster_index=0;
#			$cluster_end=(sort {$a->[0]->[0]->[4]<=>$b->[0]->[0]->[4]} map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]})[-1]->[0]->[0]->[4];
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$line[2]}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
		}
		if ($line[3]>$cluster_end ) {
			if ($cluster_index+1>$#{$CLUSTER_GENE->{$current_chr}}) {
				if ($tag_chromosome_end==0) {
					&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$CLUSTER_AS->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png") if (scalar(@gene_reads)>10) ;
				}
				$distribution->{"intergenic"}++;
				$tag_chromosome_end=1;
				next;
			}else{
				$tag_chromosome_end=0;
			}
			&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$CLUSTER_AS->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png") if (scalar(@gene_reads)>10) ;
			@gene_reads=();
			$cluster_index++;
			
			$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
			$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
			$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
			$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			while ($line[3]>$cluster_end) {
				$cluster_index++;
				if ($cluster_index>$#{$CLUSTER_GENE->{$current_chr}}){
					$cluster_index--;
					$tag_chromosome_end=1;
					last;
				}
				$cluster_end=max(map {$_->[0]->[0]->[4]} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]});
				$cluster_end_strand= (map {$_->[0]->[0]->[6]} grep {$_->[0]->[0]->[4]==$cluster_end} @{$CLUSTER_GENE->{$current_chr}->[$cluster_index]})[0];
				$cluster_start= $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[6] eq "+" ? $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$promoter_length : $CLUSTER_GENE->{$current_chr}->[$cluster_index]->[0]->[0]->[0]->[3]-$downstream_length;
				$cluster_end= $cluster_end_strand eq "-" ? $cluster_end+$promoter_length : $cluster_end+$downstream_length;
			}
		}
		$current_chr=$line[2];
		if ($line[3]<$cluster_start) {
			$distribution->{"intergenic"}++;
			push @intergenic_reads,[($line[3],$read_strand,$line[5])];
#			draw (\@intergenic_reads,);
		}else{
			push @gene_reads,[($line[3],$read_strand,$line[5],$line[9])];
		}

		if (eof(IN)) {
			if ($tag_chromosome_end==0) {
				&draw(\@gene_reads,$CLUSTER_GENE->{$current_chr}->[$cluster_index],$CLUSTER_AS->{$current_chr}->[$cluster_index],$cluster_start,$cluster_end,"$current_chr\_cluster$cluster_index.png") if (scalar(@gene_reads)>10) ;
			}
		}
	}
	close(IN);
}

sub draw {
	my ($reads_ref,$cluster,$ascluster,$Cluster_start,$Cluster_end,$BASE_NAME)=@_;
	open (OUT,">$BASE_NAME") or die $!;
	my $left=300;
	my $top=300;
	my $scale=1;
	my $scale_height=5000;
	my $scale_width=8000;
	my $Width=$scale_width/$scale+2*$left;

	my $cluster_Length=abs($Cluster_end-$Cluster_start);
	my $ratio=$scale_width/$cluster_Length;

	my $RNA_line_height=0.02;
	my $AS_line_height=0.01;
	my $read_height=0.002;
	my (@junction,@no_junction);
	&differentiate_junction($reads_ref,\@junction,\@no_junction);
	my (@tidy_junction,@tidy_no_junction);
	my $junction_reads_count=&tidy_tag_line(\@junction,\@tidy_junction,$merge_dentical_reads);
	my $nojuction_reads_count=&tidy_tag_line(\@no_junction,\@tidy_no_junction,$merge_dentical_reads);
	my $as_count=sum(map {$#{$_}+1} @{$ascluster});
	$as_count = 0 if not defined($as_count);
	
	my $RNA_count=sum(map {$#{$_}} @{$cluster});

	my $Height=($RNA_count*$RNA_line_height+(scalar(@tidy_junction)+scalar(@tidy_no_junction)+10)*$read_height+0.06)*$scale_height/$scale+2*$top;
	my $img = GD::Simple->new($Width,$Height);

	$img->font(gdSmallFont);
	$img->bgcolor('blue');
	$img->fgcolor('blue');
	$img->fontsize(40);
	my $RNA_index=0;
	my $temp_end=10e10;
	foreach my $temp_gene (@{$cluster}) {
#		foreach my $temp_RNA (@{$temp_gene}[1..$#{$temp_gene}]) {
		foreach my $temp_RNA (@{$temp_gene}) {
#			next if ($temp_RNA->[0]->[2]=~/^mRNA|miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/) ;
			next if ($temp_RNA->[0]->[2]=~/^gene|pseudogene|transposable_element_gene$/) ;
			$RNA_index-- if ($temp_RNA->[0]->[3] > $temp_end+100*$scale/$ratio) ;
			$img->penSize(5,5);
			######draw RNA information###########
			$img->fgcolor('blue');
			$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.009)*$scale_height/$scale);
			$img->angle(0);
			$img->font('Times:italic');
#			print join("\t",@{$temp_RNA->[0]}),"\n";
			$img->string(join(" ",@{$temp_RNA->[0]}[0,3,4,6,-1]));
			######draw mRNA|ncRNA length#########
			$img->fgcolor('blue');
			$img->line($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale,$left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);

			print join("\t",@{$temp_RNA->[0]}),"\n";
			######draw each CDS|three_prime_UTR|five_prime_UTR|ncExon########
			foreach my $temp_exon (@{$temp_RNA}[1..$#{$temp_RNA}]) {
				print join("\t",@{$temp_exon}),"\n";
				$img->penSize(1,1);
				if ($temp_gene->[0]->[0]->[2]=~/^mRNA$/) {
					next if ($temp_exon->[2]=~/^exon|start_codon|stop_codon$/) ;
				}
				if ($temp_exon->[2]=~/^three_prime_UTR|five_prime_UTR$/) {
					$img->bgcolor('fuchsia');
					$img->fgcolor('fuchsia');
				}else{
					$img->bgcolor('aqua');
					$img->fgcolor('aqua');
				}
				my $exon_x1=($temp_exon->[3]-$Cluster_start)*$ratio/$scale;
				my $exon_x2=($temp_exon->[4]-$Cluster_start)*$ratio/$scale;
				$img->rectangle($left+$exon_x1,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale,$left+$exon_x2,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}
			#####mark the transcript direction#####
			if ($temp_RNA->[0]->[6] eq "+") {
				$img->fgcolor('red');
				$img->penSize(3,3);
				$img->moveTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale);
				$img->angle(45);
				$img->lineTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale+50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->moveTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale+50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->angle(135);
				$img->lineTo($left+($temp_RNA->[0]->[4]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}else{
				$img->fgcolor('blue');
				$img->penSize(3,3);
				$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.01)*$scale_height/$scale);
				$img->angle(135);
				$img->lineTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale-50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->moveTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale-50,$top+($RNA_line_height*$RNA_index+0.015)*$scale_height/$scale);
				$img->angle(45);
				$img->lineTo($left+($temp_RNA->[0]->[3]-$Cluster_start)*$ratio/$scale,$top+($RNA_line_height*$RNA_index+0.02)*$scale_height/$scale);
			}
			$RNA_index++;
			$temp_end = $temp_RNA->[-1]->[4];
		}
	}

	if (@tidy_junction){
		foreach my $line (@tidy_junction) {
			map {print join("\t",@{$_}),"\n";} @{$line};
		}
	}
	if (@tidy_no_junction) {
		foreach my $line (@tidy_no_junction) {
			map {print join("\t",@{$_}),"\n";} @{$line};
		}
	}
	#############draw_RNA_frame##################
	$img->bgcolor(undef);
	$img->fgcolor('black');
	$img->penSize(1,1);
	$img->rectangle($left,$top,$left+$scale_width/$scale,$top+$RNA_index*$RNA_line_height*$scale_height/$scale);
	#############draw x-axis#####################
	$img->penSize(5,5);
	$img->fgcolor('black');
	my $locate=$RNA_index*$RNA_line_height+0.02;
	my $count = 0;
	for (my $i=$Cluster_start;$i<=$Cluster_end ;$i+=int($cluster_Length*0.01)) {
		if($count%10==0){$count++;next;}
		$img->moveTo($left+($i-$Cluster_start)*$ratio/$scale,$top+$locate*$scale_height/$scale);
		$img->lineTo($left+($i-$Cluster_start)*$ratio/$scale,$top+($locate-0.005)*$scale_height/$scale);
		$count++;
	}
	for (my $i=$Cluster_start;$i<=$Cluster_end ;$i+=int($cluster_Length*0.1)) {
		$img->moveTo($left+($i-$Cluster_start)*$ratio/$scale,$top+$locate*$scale_height/$scale);
		$img->lineTo($left+($i-$Cluster_start)*$ratio/$scale,$top+($locate-0.01)*$scale_height/$scale);
		$img->moveTo($left+($i-$Cluster_start)*$ratio/$scale-10,$top+($locate+0.015)*$scale_height/$scale);
		$img->angle(0);
#		$img->fontsize(40);
		$img->font('Times:italic');
		$img->string($i);
	}
	$img->moveTo($left,$top+$locate*$scale_height/$scale);
	$img->lineTo($left+$cluster_Length*$ratio/$scale,$top+$locate*$scale_height/$scale);

	############draw AS ######################
	$locate+=0.029;
	$img->moveTo($left,$top+$locate*$scale_height/$scale);
	$img->string("Total as: $as_count");
	$locate+=0.001;
	my $hmax=0;
	$hmax=($as_count+2)*$AS_line_height;
	$img->penSize(1,1);
	$img->angle(0);
	$img->bgcolor(undef);
	$img->fgcolor('black');
	$img->rectangle($left,$top+$locate*$scale_height/$scale,$left+$scale_width/$scale,$top+($locate+$hmax)*$scale_height/$scale);
	my $AS_index=0;
	my $as_temp_end=10e10;
	foreach my $temp_as_content (@{$ascluster}) {
		foreach my $temp_as (@{$temp_as_content}) {
			# $AS_index-- if ($temp_as->[1] > $as_temp_end+100*$scale/$ratio) ;
			$img->penSize(5,5);
			######draw RNA information###########
			$img->fgcolor('blue');
			$img->moveTo($left+($temp_as->[1]-$Cluster_start)*$ratio/$scale,$top+($AS_line_height*$AS_index+0.009+$locate)*$scale_height/$scale);
			$img->angle(0);
			$img->font('Times:italic');
	#			print join("\t",@{$temp_RNA->[0]}),"\n";
			$img->string($temp_as->[4]." ".$temp_as->[8]);
			######draw mRNA|ncRNA length#########
			$img->fgcolor('blue');
			$img->line($left+($temp_as->[1]-$Cluster_start)*$ratio/$scale,$top+($AS_line_height*$AS_index+0.015+$locate)*$scale_height/$scale,$left+($temp_as->[2]-$Cluster_start)*$ratio/$scale,$top+($AS_line_height*$AS_index+0.015+$locate)*$scale_height/$scale);

			$img->penSize(1,1);
			$img->fgcolor('black');
			$img->line($left+($temp_as->[1]-$Cluster_start)*$ratio/$scale,$top,$left+($temp_as->[1]-$Cluster_start)*$ratio/$scale,$scale_height/$scale);
			$img->line($left+($temp_as->[2]-$Cluster_start)*$ratio/$scale,$top,$left+($temp_as->[2]-$Cluster_start)*$ratio/$scale,$scale_height/$scale);

			# print join("\t",@{$temp_as}),"\n";

			$AS_index++;
			# $as_temp_end = $temp_as->[-1]->[4];
		}
	}
	############draw reads ######################
	$locate+=$AS_index*$AS_line_height+0.03;
	$img->moveTo($left,$top+$locate*$scale_height/$scale);
	$junction_reads_count="Total junction reads: $junction_reads_count";
	$img->string($junction_reads_count);
	$locate+=0.001;
	$locate+=&draw_hits($read_height,$left,$top,$scale_height,$scale_width,$scale,\@tidy_junction,$locate,$ratio,$Cluster_start,$img)+0.009;
	$img->angle(0);
	$img->moveTo($left,$top+$locate*$scale_height/$scale);
	$nojuction_reads_count=" Total exonic reads:$nojuction_reads_count";
	$img->string($nojuction_reads_count);
	$locate+=0.001;
	$locate+=&draw_hits($read_height,$left,$top,$scale_height,$scale_width,$scale,\@tidy_no_junction,$locate,$ratio,$Cluster_start,$img);
	binmode OUT;
	print OUT $img->png;
	close OUT;
#	print "Done draw $BASE_NAME\n";
}

sub differentiate_junction {
	my ($array_ref,$return_ref1,$return_ref2)=@_;
	for (my $i=0;$i<@{$array_ref} ;$i++) {
		if ($array_ref->[$i]->[2]=~/N/) {
			push @{$return_ref1},$array_ref->[$i];
		}else{
			push @{$return_ref2},$array_ref->[$i];
		}
	}
}

sub draw_hits {
	my ($Read_height,$Left,$Top,$Scale_height,$Scale_width,$Scale,$array_ref,$Locate,$cluster_length_ratio,$Cluster_start,$IMG)=@_;
	my $max=0;
	my @tidy_tag=@{$array_ref};
	$max=(scalar(@tidy_tag)+5)*$Read_height;
	$IMG->penSize(1,1);
	$IMG->angle(0);
	$IMG->bgcolor(undef);
	$IMG->fgcolor('black');
	$IMG->rectangle($Left,$Top+$Locate*$Scale_height/$Scale,$Left+$Scale_width/$Scale,$Top+($Locate+$max)*$Scale_height/$Scale);
	$IMG->penSize(5,5);
	for (my $i=0;$i<@tidy_tag ;$i++) {
		foreach my $read (@{$tidy_tag[$i]}) {
			my $START=$read->[0];
			my $tag_strand = $read->[1];
			my $matched_size=$read->[2];
			my @tag=split(/[M|N]/,$matched_size);
			my (@Reads_start,@Matched_length);
			&reads_split_0($START,$matched_size,\@Reads_start,\@Matched_length);
			for(my $j=0;$j<@Reads_start;$j++) {
				$IMG->moveTo($Left+($Reads_start[$j]-$Cluster_start)*$cluster_length_ratio/$Scale,$Top+($Read_height*($i+1)+$Locate)*$Scale_height/$Scale);
				if ($j%2==1) {
					$IMG->penSize(1,1);
					$IMG->fgcolor('green');
				}else{
					my @rgb=();
					if ($tag_strand eq "+") {
						@rgb=(255,0,0);
						if (defined($read->[3]) && $read->[3]=~/^\d+$/){
							my @start=(210,180,140);
							my @end=(94,38,18);
							@rgb=&color_gradient_ten($read->[3],1,100,\@start,\@rgb);
						}
					}else{
						@rgb=(0,0,255);
						if (defined($read->[3]) && $read->[3]=~/^\d+$/ ){
							my @start=(176,224,230);
							my @end=(11,23,70);
							@rgb=&color_gradient_ten($read->[3],1,100,\@start,\@rgb);
						}
					}
					my $index = $IMG->translate_color(@rgb) ;
					$IMG->fgcolor($index);
#					$IMG->fgcolor($tag_strand eq "+" ? 'red' : 'blue');
					$IMG->penSize(5,5);
				}
				$IMG->line($Matched_length[$j]*$cluster_length_ratio/$Scale);
#				if (defined($read->[3])){		#draw identical reads number
#					$IMG->fgcolor('black');
#					$IMG->font(gdSmallFont);
#					$IMG->fontsize(10);
#					$IMG->moveTo($Left+($Reads_start[$j]-$Cluster_start+0.4*$Matched_length[$j])*$cluster_length_ratio/$Scale,$Top+($Read_height*($i+1)+$Locate)*$Scale_height/$Scale);
#					$IMG->string($read->[3]);
#				}
			}
		}
	}
	##draw y-axis
	$IMG->fgcolor('black');
	$IMG->penSize(5,5);
	for (my $k=$Locate;$k<=$Locate+$max ;$k+=5*$Read_height) {
		$IMG->moveTo($Left-75,$Top+$k*$Scale_height/$Scale);
		$IMG->line(15);
		$IMG->moveTo($Left-50,$Top+($k+0.003)*$Scale_height/$Scale);
		$IMG->fontsize(40);
		$IMG->font('Times:italic');
		$IMG->string(int(($k-$Locate)/$Read_height+0.5));
	}
	$IMG->moveTo($Left-75,$Top+$Locate*$Scale_height/$Scale);
	$IMG->angle(90);
	$IMG->line($max*$Scale_height/$Scale);
	return $max ;
}

sub color_gradient_ten () {#datanow,data1,data2,color1,color2#######颜色随值变化的函数
	my @svg_x=@_;
	my @out_color;
	if ($svg_x[0] >=$svg_x[2]) {
		@out_color=@{$svg_x[4]};
	}
	elsif ($svg_x[0] <=$svg_x[1]) {
		@out_color=@{$svg_x[3]};
	}
	else {
		my $new_red=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[0]-$svg_x[3]->[0])+$svg_x[3]->[0]);
		my $new_gre=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[1]-$svg_x[3]->[1])+$svg_x[3]->[1]);
		my $new_blu=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($svg_x[4]->[2]-$svg_x[3]->[2])+$svg_x[3]->[2]);
		@out_color=($new_red,$new_gre,$new_blu);
	}
	return @out_color;
}

sub reads_split_0 {
	my ($start,$string,$reads_start,$matched_length)=@_;
	@{$matched_length}=split(/[M|N]/,$string);
	undef(@{$reads_start});
	push @{$reads_start},$start;
	push @{$reads_start},map {$start+sum(@{$matched_length}[0..$_])} (0..$#{$matched_length}-1) if (scalar(@{$matched_length}>1));
}

sub tidy_tag_line1{
	my ($hash_ref,$arr_ref,$Merge_identical_reads)=@_;
	my $Reads_count=scalar(@{$hash_ref});
	if ($Merge_identical_reads) {
	#################combine the identical reads ##################
		my %hash_ref_tmp;
		foreach my $s (@{$hash_ref}) {
			$hash_ref_tmp{join("\t",@{$s})}++;
		}
		@{$hash_ref}=();
		foreach my $t (sort {(split/\s+/,$a)[0]<=>(split/\s+/,$b)[0]} keys %hash_ref_tmp) {
			push @{$hash_ref},[(split(/\s+/,$t),$hash_ref_tmp{$t})];
		}
	}
	else{
		@{$hash_ref}=sort {$a->[0]<=>$b->[0]} @{$hash_ref};
	}
	################################################################
	my $line=0;
	my $end=0;
	while (@{$hash_ref}) {
		$end = 0;
		for (my $k=0;$k<@{$hash_ref} ;$k++) {
			if ($hash_ref->[$k]->[0]>$end) {#############一位加一，相当于相邻两线之间可以加1
				my @tmp_size=split(/[M|N]/,$hash_ref->[$k]->[2]);
				$end = $hash_ref->[$k]->[0] + sum(@tmp_size) - 1;
				push @{$arr_ref->[$line]},splice(@{$hash_ref},$k,1);
				$k--;
			}
		}
		$line++;
	}
	return $Reads_count;
}

sub tidy_tag_line{
	my ($hash_ref,$arr_ref,$Merge_identical_reads)=@_;
	my $Reads_count=scalar(@{$hash_ref});
	if ($Merge_identical_reads) {
	#################combine the identical reads ##################
		my %hash_ref_tmp;
		foreach my $s (@{$hash_ref}) {
			$hash_ref_tmp{join("\t",@{$s})}++;
		}
		@{$hash_ref}=();
		foreach my $t (sort {(split/\s+/,$a)[1] cmp (split/\s+/,$b)[1] or (split/\s+/,$a)[0]<=>(split/\s+/,$b)[0]} keys %hash_ref_tmp) {
			push @{$hash_ref},[(split(/\s+/,$t),$hash_ref_tmp{$t})];
		}
	}
	else{
		@{$hash_ref}=sort {$a->[1] cmp $b->[1] or $a->[0]<=>$b->[0]} @{$hash_ref};
	}
	################################################################
	my $line=0;
	my $end=0;
	while (grep {$_->[1] eq "+"} @{$hash_ref}) {
		$end = 0;
		for (my $k=0;$k<@{$hash_ref} ;$k++) {
			last if ($hash_ref->[$k]->[1] eq "-") ;
			if ($hash_ref->[$k]->[0]>$end) {#############一位加一，相当于相邻两线之间可以加1
				my @tmp_size=split(/[M|N]/,$hash_ref->[$k]->[2]);
				$end = $hash_ref->[$k]->[0] + sum(@tmp_size) - 1;
				push @{$arr_ref->[$line]},splice(@{$hash_ref},$k,1);
				$k--;
			}
		}
		$line++;
	}
	while (@{$hash_ref}) {
		$end = 0;
		for (my $k=0;$k<@{$hash_ref} ;$k++) {;
			if ($hash_ref->[$k]->[0]>$end) {#############一位加一，相当于相邻两线之间可以加1
				my @tmp_size=split(/[M|N]/,$hash_ref->[$k]->[2]);
				$end = $hash_ref->[$k]->[0] + sum(@tmp_size) - 1;
				push @{$arr_ref->[$line]},splice(@{$hash_ref},$k,1);
				$k--;
			}
		}
		$line++;
	}
	return $Reads_count;
}

sub foreach_gene{
	my ($CHR,$FG,$LG,$GENE,$DISTRIBUTION,$GENE_READS,$promoter_length)=@_;
	my ($promoter_start,$promoter_end);
	for (my $i=$FG;$i<=$LG ;$i++) {
		foreach my $temp_RNA (@{$GENE->{$CHR}->[$i]}[1..$#{$GENE->{$CHR}->[$i]}]) {
			my @Gene_reads=@{$GENE_READS};
			my ($gene_read_count,$exon_count)=(0,0);
			foreach my $reads_ref (@{$GENE_READS}) {
				my (@reads_start,@reads_center,@matched_length,@reads_count);
				&reads_split($reads_ref->[0],$reads_ref->[-1],\@reads_start,\@reads_center,\@matched_length,\@reads_count);
				if ($temp_RNA->[0]->[6] eq "+") {
					$promoter_length=$temp_RNA->[0]->[3] if ($temp_RNA->[0]->[3]-$promoter_length<0);
					($promoter_start,$promoter_end)=($temp_RNA->[0]->[3]-$promoter_length,$temp_RNA->[0]->[3]);
					foreach my $index (0..$#reads_center) {
						if ($reads_center[$index]<$promoter_start) {
							$DISTRIBUTION->{"intergenic"}+=$reads_count[$index];
						}
						elsif ( $reads_center[$index]>=$promoter_start && $reads_center[$index]<$promoter_end ) {
							$DISTRIBUTION->{"promoter"}+=$reads_count[$index];
						}
						elsif($reads_center[$index]>$temp_RNA->[0]->[4]){
							$DISTRIBUTION->{"intergenic"}+=$reads_count[$index];
						}
						elsif($reads_center[$index]>=$temp_RNA->[0]->[3] && $reads_center[$index]<=$temp_RNA->[0]->[4]){
							$gene_read_count+=$reads_count[$index];
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
									next if ($temp_line->[2]=~/^exon|start_codon|stop_codon$/) ;
								}
								if ($reads_center[$index]>=$temp_line->[3] && $reads_center[$index]<=$temp_line->[4]) {
									$exon_count+=$reads_count[$index];
									$DISTRIBUTION->{$temp_line->[2]}+=$reads_count[$index];
								}
							}
						}
					}
				}
				else{
					($promoter_start,$promoter_end)=($temp_RNA->[0]->[4],$temp_RNA->[0]->[4]+$promoter_length);
					foreach my $index (0..$#reads_center) {
						if ($reads_center[$index]<$temp_RNA->[0]->[3]) {
							$DISTRIBUTION->{"intergenic"}+=$reads_count[$index];
						}
						elsif ( $reads_center[$index]>$promoter_start && $reads_center[$index]<=$promoter_end ) {
							$DISTRIBUTION->{"promoter"}+=$reads_count[$index];
						}
						elsif($reads_center[$index]>$promoter_end){
							$DISTRIBUTION->{"intergenic"}+=$reads_count[$index];
						}
						elsif($reads_center[$index]>=$temp_RNA->[0]->[3] && $reads_center[$index]<=$temp_RNA->[0]->[4]){
							$gene_read_count+=$reads_count[$index];
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
									next if ($temp_line->[2]=~/^exon|start_codon|stop_codon$/) ;
								}
								if ($reads_center[$index]>=$temp_line->[3] && $reads_center[$index]<=$temp_line->[4]) {
									$exon_count+=$reads_count[$index];
									$DISTRIBUTION->{$temp_line->[2]}+=$reads_count[$index];
								}
							}
						}
					}
				}
			}
			$DISTRIBUTION->{"intron"}+=$gene_read_count-$exon_count;
		}
	}
}

sub foreach_gene_depth{
	my ($CHR,$RS,$RE,$FG,$LG,$GENE,$GENE_READS,$FILE_HANDEL,$promoter_length)=@_;
	my ($promoter_start,$promoter_end);
	my %Depth=();
	foreach my $reads_ref (@{$GENE_READS}) {
		map {$Depth{$_}++} ($reads_ref->[0]..$reads_ref->[1]);
	}
	for (my $i=$FG;$i<=$LG ;$i++) {
		foreach my $temp_RNA (@{$GENE->{$CHR}->[$i]}[1..$#{$GENE->{$CHR}->[$i]}]) {
			my %sum=("intergenic"=>0,"promoter"=>0,"intron"=>0,"exon"=>0,"CDS"=>0,"five_prime_UTR"=>0,"three_prime_UTR"=>0);
			my %region_length=("promoter"=>0,"intron"=>0,"exon"=>0,"CDS"=>0,"five_prime_UTR"=>0,"three_prime_UTR"=>0);
			my %depth=%Depth;
			if ($temp_RNA->[0]->[6] eq "+") {
				$promoter_length=$temp_RNA->[0]->[3] if ($temp_RNA->[0]->[3]-$promoter_length<0);
				($promoter_start,$promoter_end)=($temp_RNA->[0]->[3]-$promoter_length,$temp_RNA->[0]->[3]);
				$region_length{"promoter"}+=$promoter_length;
				foreach my $pos (sort {$a<=>$b} keys %depth) {
					if ($pos<= $promoter_start) {
						$sum{"intergenic"}+=delete($depth{$pos})if (exists($depth{$pos})) ;
					}
					elsif($pos>$promoter_start && $pos<$promoter_end){
						$sum{"promoter"}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
					}
					elsif($pos>$temp_RNA->[0]->[4]){
						$sum{"intergenic"}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
					}
					elsif($pos>=$promoter_end && $pos<=$temp_RNA->[0]->[4]){
						if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								next if ($temp_line->[2]=~/^exon|start_codon|stop_codon$/) ;
								if ($pos>=$temp_line->[3] && $pos<=$temp_line->[4]) {
									$sum{$temp_line->[2]}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
								}
								$region_length{$temp_line->[2]}+=abs($temp_line->[4]-$temp_line->[3])+1;
							}
						}
						elsif($temp_RNA->[0]->[2]=~/^miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/){
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								if ($pos>=$temp_line->[3] && $pos<=$temp_line->[4]) {
									$sum{$temp_line->[2]}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
								}
								$region_length{$temp_line->[2]}+=abs($temp_line->[4]-$temp_line->[3])+1;
							}
						}
					}
				}
			}
			else{
				($promoter_start,$promoter_end)=($temp_RNA->[0]->[4],$temp_RNA->[0]->[4]+$promoter_length);
				$region_length{"promoter"}+=$promoter_length;
				foreach my $pos (sort {$a<=>$b} keys %depth) {
					if ($pos< $temp_RNA->[0]->[3]) {
						$sum{"intergenic"}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
					}
					elsif($pos>=$temp_RNA->[0]->[3] && $pos<=$temp_RNA->[0]->[4]){
						if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								next if ($temp_line->[2]=~/^exon|start_codon|stop_codon$/) ;
								$sum{$temp_line->[2]}+=delete($depth{$pos}) if ($pos>=$temp_line->[3] && $pos<=$temp_line->[4]) ;
								$region_length{$temp_line->[2]}+=abs($temp_line->[4]-$temp_line->[3])+1;
							}
						}
						elsif($temp_RNA->[0]->[2]=~/^miRNA|mRNA_TE_gene|ncRNA|pseudogenic_transcript|rRNA|snoRNA|snRNA|tRNA|transcript$/){
							foreach my $temp_line (@{$temp_RNA}[1..$#{$temp_RNA}]) {
								$sum{$temp_line->[2]}+=delete($depth{$pos}) if ($pos>=$temp_line->[3] && $pos<=$temp_line->[4]) ;
								$region_length{$temp_line->[2]}+=abs($temp_line->[4]-$temp_line->[3])+1;
							}
						}
					}
					elsif($pos>$promoter_start && $pos<$promoter_end){
						$sum{"promoter"}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
					}
					elsif($pos>=$promoter_end){
						$sum{"intergenic"}+=delete($depth{$pos}) if (exists($depth{$pos})) ;
					}
				}
			}
			$sum{"intron"}=sum(values(%depth)) if (%depth);
			$region_length{"intron"}+=abs($temp_RNA->[0]->[4]-$temp_RNA->[0]->[3])+1-$region_length{"CDS"}-$region_length{"three_prime_UTR"}-$region_length{"five_prime_UTR"}-$region_length{"exon"};
			my %average_depth=();
			delete($sum{"intergenic"});
			foreach my $temp_region (keys %sum) {
#				next if ($temp_region eq "intergenic") ;
				if ($region_length{$temp_region}==0) {
					$average_depth{$temp_region}=0;
				}else{
					$average_depth{$temp_region}=sprintf("%10.5f",$sum{$temp_region}/$region_length{$temp_region});
				}
			}
			print $FILE_HANDEL join("\t",map {$sum{$_}} sort keys %sum),"\n";
			print $FILE_HANDEL join("\t",map {$region_length{$_}} sort keys %region_length),"\n";
			print $FILE_HANDEL join("\t",map {$average_depth{$_}} sort keys %average_depth),"\n";
		}
	}
}

sub reads_split{
	my ($start,$string,$reads_start,$reads_center,$matched_length,$reads_count)=@_;
	my @array=split(/[M|N]/,$string);
#	print "@array\n";
	my @matched_index=grep {$_%2==0} (0..$#array);
	@{$matched_length}=map {$array[$_]} @matched_index;
	my $read_length=sum(@{$matched_length});
	@{$reads_count}=map {$_/$read_length} @{$matched_length};
	undef(@{$reads_start});
	push @{$reads_start},$start;
	push @{$reads_start},map {$start+sum(@array[0..$_])} grep {$_%2==1} (0..$#array) if (scalar(@array>1));
	@{$reads_center}=map {$reads_start->[$_]-1+int(($matched_length->[$_]+1)/2)} (0..$#{$reads_start});
}

#sub get_intron{
#	my ($RNA,$Intron)=@_;
#	undef(@{$Intron});
#	my $i=0;
#	foreach my $content (@{$RNA}[1..$#{$RNA}]) {
#		$i++;
#		if ($temp_RNA->[0]->[2]=~/^mRNA$/) {
#			next if ($temp_line->[2]=~/^exon|start_codon|stop_codon$/) ;
#		}
#		
#		push @{$Intron},[($content->[],$content->[],$content->[6])];
#	}
#}

sub overlap_indentify {
	my ($start1,$end1,$start2,$end2)=@_;
	my %hash=();
	map {$hash{$_}++} ($start1..$end1);
	map {$hash{$_}++} ($start2..$end2);
	my @overlap=grep {$hash{$_}==2} keys(%hash);
	return scalar(@overlap);
}

#open (OUT,">$outfile") || die "Can't creat $outfile\n";		#open a file as writing
#close(OUT);

sub cat
#function:quit redundance
#input:($para,@array), para is the merge length
#output:(@array),
#for example (0,1,3,4,7,5,8)->(1,3,4,8) (1,1,3,4,7,5,8)->(1,8)
{
	my($merge,@input) = @_;
	my $i = 0;
	my @output = ();
	my %hash = ();
	my $each = 0;
	my $begin = "";
	my $end = 0;
	for ($i=0;$i<@input;$i+=2){
		my $Qb = $input[$i];
		my $Qe = $input[$i+1];
		if($Qb > $Qe) { next; }
		if(defined($hash{$Qb}))	{
			if($hash{$Qb} < $Qe) {
				$hash{$Qb} = $Qe;
			}
		}
		else {
			$hash{$Qb} = $Qe;
		}
		$Qb = 0;
	}
	foreach $each (sort {$a <=> $b} keys %hash){
		if($begin eq ""){
			$begin = $each;
			$end = $hash{$each};
		}else{
			if($hash{$each} > $end){
				if($each > $end + $merge){
					push(@output,$begin);
					push(@output,$end);
					$begin = $each;
					$end = $hash{$each};
				}
				else {
					$end = $hash{$each};
				}
			}
		}
	}
	if(keys %hash > 0){
		push(@output,$begin);
		push(@output,$end);
	}
	%hash = ();
	return(@output);
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