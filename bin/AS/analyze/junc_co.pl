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
GetOptions( \%opts, "sj=s","name=s","o=s", "h" );

if (   !defined( $opts{sj} )
    || !defined( $opts{name} )
    || !defined( $opts{o} )
    || defined( $opts{h} ) )
{
    print <<"Usage End.";

        Version:$ver

    Usage:perl $0

        -sj       sample sj list           ex:file1,file2,file3;

        -name       sample name list       ex:file1_name,file2_name,file3_name;same order with sj

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
my $current_dir = `pwd`;chomp($current_dir);

my @sj_group  = split(/,/,$opts{sj});
my @name_group  = split(/,/,$opts{name});
my $out_file = $opts{o};

open OUT, ">$out_file" || die;
print OUT "#Chr\tStart\tEnd\tStrand\t",join("\t",@name_group),"\n";



my %junction    = ();    #junctions hash
&load_all_sj();

&co_Junc();


#子过程-----co_Junc
sub co_Junc {

    foreach my $key (sort keys %junction){
        print OUT $key;
        for(my $i=0;$i<=$#name_group;$i++){
            if(not defined($junction{$key}->{$name_group[$i]})){
                $junction{$key}->{$name_group[$i]} = 0;
            }
            print OUT "\t$junction{$key}->{$name_group[$i]}";
        }
        print OUT "\n";
    }

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
            my $key = "$chr\t$start\t$end\t$strand";
            $junction{$key}->{$name_group[$i]}=$readsNum;
        }
        close SJ;
    }
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