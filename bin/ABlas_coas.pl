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
GetOptions( \%opts, "i=s","name=s","p=s", "o=s", "h" );

if (   !defined( $opts{i} )
	|| !defined( $opts{name} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-i          input dir              summary directory
		
		-name       sample name list       ex:file1_name,file2_name,file3_name;

		-o          out file               must be given;

		-p          postfix                default is ".miso_summary"

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

my $input_dir = $opts{i};
$input_dir = &AbsolutePath("dir",$input_dir);

my $postfix='.miso_summary';
$postfix=$opts{p} if defined $opts{p};

my @file_group;
my @name_group  = split(/,/,$opts{name});
my $out_file = $opts{o};

map{push @file_group,"$input_dir/$_$postfix"} @name_group;
map{print "$input_dir/$_$postfix\n"} @name_group;
if(scalar(@file_group) != scalar(@name_group)){
	print "as file nums must be name as name nums.\n";
	exit;
}

open OUT, ">$out_file" || die;
print OUT "#Chromosome\tIsoform\tType\tStart\tEnd\tStrand\tFeature";
foreach my $asname(@name_group){
	print OUT "\t$asname\:Alt_sj_reads\t$asname\:Model_sj_reads\t$asname\:Alt_avg_sj_reads\t$asname\:Model_avg_sj_reads\t$asname\:Alt_ratio";
}
print OUT "\n";

my %as = ();

for(my $i=0;$i<=$#file_group;$i++){
	open AS, $file_group[$i] || die;
	while (<AS>) {              #把所有的junction添加进hash
		#NW_001867363.1	rna14966,rna14967	SE	3943355	3940584	5	-	5,22	flanks "3943355^,3940584-"; structure "0,1-2^"; splice_chain ",3941587-3941431^"; splice_position "3941587,3941431";	0.27
		chomp;
		next if ( $_ =~ /^#/i );
		my @line = split(/\t/);
		my $key = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line[8]";
		my $alt_avg=&avgsj($line[5]);
		$alt_avg=sprintf("%.2f",$alt_avg);
		my $modle_avg=&avgsj($line[7]);
		$modle_avg=sprintf("%.2f",$modle_avg);
		my $value = "$line[5]\t$line[7]\t$alt_avg\t$modle_avg\t$line[9]";
		$as{$key}->{$name_group[$i]} = $value;
	}
	close AS;
}



foreach my $key (sort keys %as){
	print OUT $key;
	for(my $i=0;$i<=$#name_group;$i++){
		if(not defined($as{$key}->{$name_group[$i]})){
			$as{$key}->{$name_group[$i]} = "-\t-\t-\t-\t-";
		}
		print OUT "\t$as{$key}->{$name_group[$i]}";
	}
	print OUT "\n";
}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	my $avg = (sum(@num)) / (scalar(@num)) ;
	return $avg;
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