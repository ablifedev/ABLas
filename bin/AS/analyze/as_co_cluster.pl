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
GetOptions( \%opts, "cluster=s","sj=s","ir=s","name=s","o=s", "h" );

if (   !defined( $opts{cluster} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{ir} )
	|| !defined( $opts{name} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-cluster       as cluster file          ex:cluster.txt
		
		-sj            sample sj list           ex:file1,file2,file3;

		-ir            sample ir list           ex:file1,file2,file3;

		-name          sample name list         ex:file1_name,file2_name,file3_name;

		-o             out file                 must be given;

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################
my $current_dir = `pwd`;chomp($current_dir);

# my @as_group  = split(/,/,$opts{as});
my $cluster_file = $opts{cluster};
my @sj_group  = split(/,/,$opts{sj});
my @ir_group  = split(/,/,$opts{ir});
my @name_group  = split(/,/,$opts{name});
my $out_file = $opts{o};

print $out_file,"\n";

open OUT, ">$out_file" || die;
print OUT "#Cluster_name\t";
print OUT "Chr\tStrand\tss1\tss2\tType";
foreach my $asname(@name_group){
	print OUT "\t$asname\:Total\t$asname\:Sj\t$asname\:otherSj\t$asname\:Ratio";
}
print OUT "\n";


my %junction    = ();    #junctions hash
&load_all_sj();

my %ir    = ();    #ir hash
&load_all_ir();

&load_AS();


sub load_AS {
	open AS, $cluster_file || die;
	my %total=();
	my %as=();
	my $flag=0;
	my $end="";
	my $cluster="";
	while (<AS>) {              #把所有的junction添加进hash
		# >cluster_mxe_chr10:uc001lae.4
		# chr10	+	114917828;114920450	114920378;114925314	model
		# chr10	+	114917828;114918476	114918426;114925314	MXE
		chomp;
		next if ( $_ =~ /^#/ );
		my @line = split(/\t/);
		if(/^>cluster_(\w+)_/){
			if($flag==0){
				$flag=1;
				# print OUT $_,"\n";
				$cluster = $_;
				$cluster =~s/^>//;
				$end=$1;
				next;
			}
			#TODO output
			foreach my $splice (keys %as){
				print OUT $cluster,"\t",$splice;
				for(my $i=0;$i<=$#name_group;$i++){
					my $ratio = 0;
					$ratio = sprintf("%.2f",$as{$splice}{$name_group[$i]} / $total{$name_group[$i]}) if $total{$name_group[$i]} > 0;
					my $other = $total{$name_group[$i]} - $as{$splice}{$name_group[$i]};
					print OUT "\t$total{$name_group[$i]}\t$as{$splice}{$name_group[$i]}\t$other\t$ratio";
				}
				print OUT "\n";
			}
			
			%as=();
			%total=();
			$end=$1;
			# print OUT $_,"\n";
			$cluster = $_;
			$cluster =~s/^>//;
		}elsif($line[4] eq "IntronR"){
			for(my $i=0;$i<=$#name_group;$i++){
				$ir{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"} = 0 if not defined($ir{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"});
				$as{$_}{$name_group[$i]}=$ir{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"};
				$total{$name_group[$i]}+=$ir{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"};
			}
		}
		elsif($end ne "mxe"){
			for(my $i=0;$i<=$#name_group;$i++){
				$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"} = 0 if not defined($junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"});
				$as{$_}{$name_group[$i]}=$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"};
				$total{$name_group[$i]}+=$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$line[2]:$line[3]"};
			}
		}elsif($end eq "mxe"){
			for(my $i=0;$i<=$#name_group;$i++){
				my ($s1,$s2)=split(/;/,$line[2]);
				my ($e1,$e2)=split(/;/,$line[3]);
				$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s1:$e1"} = 0 if not defined($junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s1:$e1"});
				$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s2:$e2"} = 0 if not defined($junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s2:$e2"});
				my $avgsj = ($junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s1:$e1"}+$junction{$name_group[$i]}{$line[0]}{$line[1]}->{"$s2:$e2"})/2;
				$total{$name_group[$i]}+=$avgsj;
				$as{$_}{$name_group[$i]}=$avgsj;
			}
		}


	}
	close AS;
}



###子过程-----load_all_sj($sj_file)
sub load_all_sj {
	for(my $i=0;$i<=$#sj_group;$i++){
		open SJ, $sj_group[$i] || die;
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
			$start = $left + $leftsize;
			$end   = $right - $rightsize + 1;
			if(defined($junction{$name_group[$i]}{$chr}{$strand}->{"$start:$end"})){
				$junction{$name_group[$i]}{$chr}{$strand}->{"$start:$end"} += $readsNum;
			}else{
				$junction{$name_group[$i]}{$chr}{$strand}->{"$start:$end"} = $readsNum;
			}
		}
		close SJ;
		#调试信息
		print "done reading SJ ...........", "\n\n";
	}
}


###子过程-----load_all_sj($sj_file)
sub load_all_ir {
	for(my $i=0;$i<=$#ir_group;$i++){
		open IR, $ir_group[$i] || die;
		while (<IR>) {              #把所有的junction添加进hash
			#chr    ss1     ss2     strand  IntronR ratio   intron_ave_boudaryreads intron_sj       isoform
			#chr22   36587845        36587847        -       IntronR 0.166666666666667       2       10      uc003aox.3
			chomp;
			next if ( $_ =~ /^#/ );
			my @line=split(/\t/);
			if(defined($ir{$name_group[$i]}{$line[0]}{$line[3]}->{"$line[1]:$line[2]"})){
				$ir{$name_group[$i]}{$line[0]}{$line[3]}->{"$line[1]:$line[2]"} += $line[8];
			}else{
				$ir{$name_group[$i]}{$line[0]}{$line[3]}->{"$line[1]:$line[2]"} = $line[8];
			}
		}
		close IR;
		#调试信息
		print "done reading IR ...........", "\n\n";
	}
}

sub sumsj {
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	return sum(@num);

}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/;/,$sjlist);
	my $avg = (sum(@num)) / (scalar(@num)) ;
	return sprintf("%.1f",$avg);
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