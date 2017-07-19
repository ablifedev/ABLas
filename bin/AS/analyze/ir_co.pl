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
GetOptions( \%opts, "as=s", "exp=s", "name=s", "l=s", "o=s", "h" );

if (   !defined( $opts{as} )
    || !defined( $opts{exp} )
    || !defined( $opts{name} )
    || !defined( $opts{o} )
    || defined( $opts{h} ) )
{
    print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-as          sample as list             ex:file1,file2,file3;

		-exp       sample exp list           ex:file1,file2,file3;

		-name       sample name list       ex:file1_name,file2_name,file3_name;

		-l        label,known|novel       default is "novel"

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
my $current_dir = `pwd`;
chomp($current_dir);
my $label = $opts{l} || "novel";
$opts{as} =~ s/,$//;
$opts{exp} =~ s/,$//;
$opts{name} =~ s/,$//;

my @as_group   = split( /,/, $opts{as} );
my @exp_group  = split( /,/, $opts{exp} );
my @name_group = split( /,/, $opts{name} );
my $out_file   = $opts{o};

open OUT, ">$out_file" || die;
if ( $label ne "novel" ) {
    print OUT "#ID\tChr\tAss1\tAss2\tStrand\tType\tIsoform\tKorN";
    foreach my $asname (@name_group) {
        print OUT
"\t$asname\:Ratio\t$asname\:intron_ave_base\t$asname\:ends_exon_ave_base\t$asname\:isfiltered";
    }
    print OUT "\n";
}

my %as     = ();
my %filter = ();
&load_AS();

my %exphash = ();    #exphash
&load_all_exp();

&co_AS();

sub load_AS {
    for ( my $i = 0 ; $i <= $#as_group ; $i++ ) {
        open AS, $as_group[$i] || die;
        while (<AS>) {    #把所有的junction添加进hash
             #Chr01   8317    9654    -       IntronR 0.42    1.31    3.10    PAC:26823939
            chomp;
            next if ( $_ =~ /^#/ );
            my @line = split(/\t/);
            my $key =
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[-1]\t$label";
            $as{$key}{"as"} = 1;
            $filter{$key}{ $name_group[$i] } = 1;
        }
        close AS;
    }
}

#子过程-----load_AS
sub co_AS {
    my $count=1;
    foreach my $key ( sort keys %as ) {
        my $id = $label."_IR_".$count;
        $count++;
        print OUT $id,"\t",$key;
        for ( my $i = 0 ; $i <= $#name_group ; $i++ ) {
            $exphash{$key}->{ $name_group[$i] } = "-\t-\t-\t-"
              if not defined( $exphash{$key}->{ $name_group[$i] } );
            print OUT "\t$exphash{$key}->{$name_group[$i]}";
        }
        print OUT "\n";
    }

}

###子过程-----load_all_sj($sj_file)
sub load_all_exp {
    for ( my $i = 0 ; $i <= $#exp_group ; $i++ ) {
        open EXP, $exp_group[$i] || die;
        while (<EXP>) {

   #Chr01   75223   75307   +       -       0.09    18.53   215.34  PAC:26821253
            chomp;
            next if ( $_ =~ /^#/i );
            my @line = split(/\t/);
            my $key =
"$line[0]\t$line[1]\t$line[2]\t$line[3]\tIntronR\t$line[-1]\t$label";
            my $filterflag =
              defined( $filter{$key}{ $name_group[$i] } ) ? "y" : "n";
            my $value = "";
            next if not defined( $line[6] );
            next if not defined( $line[7] );

            if ( $line[6] + $line[7] == 0 ) {
                $value = "-\t-\t-\t-";
            }
            else {
                $value = "$line[5]\t$line[6]\t$line[7]\t$filterflag";
            }
            $exphash{$key}->{ $name_group[$i] } = $value;
        }
        close EXP;

        #调试信息
        print "done reading EXP ...........", "\n\n";
    }
}

sub sumsj {
    my ($sjlist) = @_;
    my @num = split( /,/, $sjlist );
    return sum(@num);

}

sub avgsj {
    my ($sjlist) = @_;
    my @num = split( /,/, $sjlist );
    my $avg = ( sum(@num) ) / ( scalar(@num) );
    return sprintf( "%.1f", $avg );
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
