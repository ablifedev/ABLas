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
GetOptions( \%opts, "as=s","sj=s","o=s", "h" );

if (   !defined( $opts{as} )
	|| !defined( $opts{sj} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-as        as file        must be given;

		-sj        all sj                     must be given;

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

my $as_file = $opts{as};
my $sj_file = $opts{sj};
my $out_file = $opts{o};
open OUT, ">$out_file" || die;

my %junction    = ();    #junctions hash
&load_all_sj($sj_file);

my %gene=();

&load_AS($as_file);


#子过程-----load_AS
sub load_AS {
	# Chr01   21094   23082   -       15355   23082   A3SS    0.5     50      50      PAC:26820096
	# Chr01   76405   76694   -       76405,76599     76485,76694     ES      0.5     50      50      PAC:26821253
	my ($asfile) = @_;
	my $count = 0;
	open AS, $asfile    || die;
	while (<AS>) {
		chomp;
		print OUT $_,"\n" if ( $_ =~ /^\#/ );
		next if ( $_ =~ /^\#/ );
		my @line= split /\t/;
		if(not defined($gene{$line[7]}{$line[6]})){
			$gene{$line[7]}{$line[6]}=0;
		}else{
			$gene{$line[7]}{$line[6]}++;
		}
		print OUT ">",$_,"\t",$line[6],"_",$gene{$line[7]}{$line[6]},"\n";
		my $chr = $line[0];
		my $type = $line[6];
		my $strand = $line[3];

		if($type eq "ES"){
			$junction{$chr}{$strand}->{"$line[1]:$line[2]"} = 0 if not defined($junction{$chr}{$strand}->{"$line[1]:$line[2]"});
			my $alt_sj=$junction{$chr}{$strand}->{"$line[1]:$line[2]"};
			print OUT "$chr\t$strand\t$line[1]\t$line[2]\t",$junction{$chr}{$strand}->{"$line[1]:$line[2]"},"\talt\n";
			my ($start1,$start2)=split(/;/,$line[4]);
			my ($end1,$end2)=split(/;/,$line[5]);
			$junction{$chr}{$strand}->{"$start1:$end1"} = 0 if not defined($junction{$chr}{$strand}->{"$start1:$end1"});
			$junction{$chr}{$strand}->{"$start2:$end2"} = 0 if not defined($junction{$chr}{$strand}->{"$start2:$end2"});
			my $model_sj1=$junction{$chr}{$strand}->{"$start1:$end1"};
			my $model_sj2=$junction{$chr}{$strand}->{"$start2:$end2"};
			print OUT "$chr\t$strand\t$start1\t$end1\t",$junction{$chr}{$strand}->{"$start1:$end1"},"\tmodel\n";
			print OUT "$chr\t$strand\t$start2\t$end2\t",$junction{$chr}{$strand}->{"$start2:$end2"},"\tmodel\n";
			
		}
		# Chr01   446749,448113   448113,448250   -       446749  448250  cassetteExon    0.5     100     100     PAC:26821290
		elsif($type eq "cassetteExon"){
			$junction{$chr}{$strand}->{"$line[4]:$line[5]"} = 0 if not defined($junction{$chr}{$strand}->{"$line[4]:$line[5]"});
			my $model_sj=$junction{$chr}{$strand}->{"$line[4]:$line[5]"};
			print OUT "$chr\t$strand\t$line[4]\t$line[5]\t",$junction{$chr}{$strand}->{"$line[4]:$line[5]"},"\tmodel\n";
			my ($start1,$start2)=split(/;/,$line[1]);
			my ($end1,$end2)=split(/;/,$line[2]);
			$junction{$chr}{$strand}->{"$start1:$end1"} = 0 if not defined($junction{$chr}{$strand}->{"$start1:$end1"});
			$junction{$chr}{$strand}->{"$start2:$end2"} = 0 if not defined($junction{$chr}{$strand}->{"$start2:$end2"});
			my $alt_sj1=$junction{$chr}{$strand}->{"$start1:$end1"};
			my $alt_sj2=$junction{$chr}{$strand}->{"$start2:$end2"};
			print OUT "$chr\t$strand\t$start1\t$end1\t",$junction{$chr}{$strand}->{"$start1:$end1"},"\talt\n";
			print OUT "$chr\t$strand\t$start2\t$end2\t",$junction{$chr}{$strand}->{"$start2:$end2"},"\talt\n";
		}
		# Chr01   351716,350505   352170,351664   -       351992,350505   352170,351928   MXE     0.5     100     100     PAC:26823778
		elsif($type eq "MXE"){
			my ($start1,$start2)=split(/;/,$line[1]);
			my ($end1,$end2)=split(/;/,$line[2]);
			$junction{$chr}{$strand}->{"$start1:$end1"} = 0 if not defined($junction{$chr}{$strand}->{"$start1:$end1"});
			$junction{$chr}{$strand}->{"$start2:$end2"} = 0 if not defined($junction{$chr}{$strand}->{"$start2:$end2"});
			my $alt_sj1=$junction{$chr}{$strand}->{"$start1:$end1"};
			my $alt_sj2=$junction{$chr}{$strand}->{"$start2:$end2"};
			print OUT "$chr\t$strand\t$start1\t$end1\t",$junction{$chr}{$strand}->{"$start1:$end1"},"\talt\n";
			print OUT "$chr\t$strand\t$start2\t$end2\t",$junction{$chr}{$strand}->{"$start2:$end2"},"\talt\n";

			($start1,$start2)=split(/;/,$line[4]);
			($end1,$end2)=split(/;/,$line[5]);
			$junction{$chr}{$strand}->{"$start1:$end1"} = 0 if not defined($junction{$chr}{$strand}->{"$start1:$end1"});
			$junction{$chr}{$strand}->{"$start2:$end2"} = 0 if not defined($junction{$chr}{$strand}->{"$start2:$end2"});
			my $model_sj1=$junction{$chr}{$strand}->{"$start1:$end1"};
			my $model_sj2=$junction{$chr}{$strand}->{"$start2:$end2"};
			print OUT "$chr\t$strand\t$start1\t$end1\t",$junction{$chr}{$strand}->{"$start1:$end1"},"\tmodel\n";
			print OUT "$chr\t$strand\t$start2\t$end2\t",$junction{$chr}{$strand}->{"$start2:$end2"},"\tmodel\n";
		}
		else{
			$junction{$chr}{$strand}->{"$line[1]:$line[2]"} = 0 if not defined($junction{$chr}{$strand}->{"$line[1]:$line[2]"});
			my $alt_sj=$junction{$chr}{$strand}->{"$line[1]:$line[2]"};
			print OUT "$chr\t$strand\t$line[1]\t$line[2]\t",$junction{$chr}{$strand}->{"$line[1]:$line[2]"},"\talt\n";
			$junction{$chr}{$strand}->{"$line[4]:$line[5]"} = 0 if not defined($junction{$chr}{$strand}->{"$line[4]:$line[5]"});
			my $model_sj=$junction{$chr}{$strand}->{"$line[4]:$line[5]"};
			print OUT "$chr\t$strand\t$line[4]\t$line[5]\t",$junction{$chr}{$strand}->{"$line[4]:$line[5]"},"\tmodel\n";
		}

		$count++;

	}
	close AS;
	#调试信息
	print "done reading AS ...........", "\n\n";
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
		$start = $left + $leftsize;
		$end   = $right - $rightsize + 1;
		if(defined($junction{$chr}{$strand}->{"$start:$end"})){
			$junction{$chr}{$strand}->{"$start:$end"} += $readsNum;
		}else{
			$junction{$chr}{$strand}->{"$start:$end"} = $readsNum;
		}
	}
	close SJ;
	#调试信息
	print "done reading SJ ...........", "\n\n";
}

sub sumsj {
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
	return sum(@num);

}

sub avgsj{
	my ($sjlist) = @_;
	my @num = split(/,/,$sjlist);
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