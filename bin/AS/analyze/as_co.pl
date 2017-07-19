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
GetOptions( \%opts, "as=s", "sj=s", "name=s", "l=s", "o=s", "h" );

if (   !defined( $opts{as} )
    || !defined( $opts{sj} )
    || !defined( $opts{name} )
    || !defined( $opts{o} )
    || defined( $opts{h} ) )
{
    print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-as          sample as list             ex:file1,file2,file3;

		-sj       sample sj list           ex:file1,file2,file3;

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
$opts{sj} =~ s/,$//;
$opts{name} =~ s/,$//;
my @as_group   = split( /,/, $opts{as} );
my @sj_group   = split( /,/, $opts{sj} );
my @name_group = split( /,/, $opts{name} );
my $out_file   = $opts{o};

open OUT, ">$out_file" || die;
if ( $label ne "novel" ) {
    print OUT "#ID\tChr\tAss1\tAss2\tStrand\tMss1\tMss2\tType\tIsoform\tKorN";
    foreach my $asname (@name_group) {
        print OUT
"\t$asname\:Ratio\t$asname\:Alt_average_sj\t$asname\:Model_average_sj\t$asname\:Alt_sj\t$asname\:Model_sj\t$asname\:isfiltered";
    }
    print OUT "\n";
}

my %junction = ();    #junctions hash
&load_all_sj();

my %as     = ();
my %filter = ();
&load_AS();
&co_AS();

sub load_AS {
    for ( my $i = 0 ; $i <= $#as_group ; $i++ ) {
        open AS, $as_group[$i] || die;
        while (<AS>) {    #把所有的junction添加进hash
             #Chr01   3020557 3021014 -       3020557 3021000 A5SS    0.16    12      61      12      61      PAC:26820620
            chomp;
            next if ( $_ =~ /^#/ );
            my @line = split(/\t/);
            my $key =
"$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[-1]\t$label";
            $as{$key}{"as"} = 1;
            $filter{$key}{ $name_group[$i] } = 1;
        }
        close AS;
    }
}

#子过程-----load_AS
sub co_AS {

    foreach my $key ( sort keys %as ) {
        for ( my $i = 0 ; $i <= $#name_group ; $i++ ) {
            my $filterflag =
              defined( $filter{$key}{ $name_group[$i] } ) ? "y" : "n";
            my @line   = split( /\t/, $key );
            my $chr    = $line[0];
            my $type   = $line[6];
            my $strand = $line[3];

            if ( $type eq "ES" ) {
                $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[1]:$line[2]"} = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$line[1]:$line[2]"} );
                my $alt_sj = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[1]:$line[2]"};
                my ( $start1, $start2 ) = split( /;/, $line[4] );
                my ( $end1,   $end2 )   = split( /;/, $line[5] );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start1:$end1"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start1:$end1"} );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start2:$end2"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start2:$end2"} );
                my $model_sj1 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start1:$end1"};
                my $model_sj2 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start2:$end2"};
                my $model_sj = "$model_sj1;$model_sj2";

                my $alt_avg_sj   = $alt_sj;
                my $model_avg_sj = &avgsj($model_sj);

                if ( $alt_avg_sj + $model_avg_sj == 0 ) {
                    $as{$key}->{ $name_group[$i] } = "-\t-\t-\t-\t-\t-";
                    next;
                }
                my $alt_ratio = sprintf( "%.2f",
                    $alt_avg_sj / ( $alt_avg_sj + $model_avg_sj ) );

                $as{$key}->{ $name_group[$i] } =
                    $alt_ratio . "\t"
                  . $alt_avg_sj . "\t"
                  . $model_avg_sj . "\t"
                  . $alt_sj . "\t"
                  . $model_sj . "\t"
                  . $filterflag;
            }

# Chr01   446749,448113   448113,448250   -       446749  448250  cassetteExon    0.5     100     100     PAC:26821290
            elsif ( $type eq "cassetteExon" ) {
                $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[4]:$line[5]"} = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$line[4]:$line[5]"} );
                my $model_sj = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[4]:$line[5]"};
                my ( $start1, $start2 ) = split( /;/, $line[1] );
                my ( $end1,   $end2 )   = split( /;/, $line[2] );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start1:$end1"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start1:$end1"} );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start2:$end2"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start2:$end2"} );
                my $alt_sj1 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start1:$end1"};
                my $alt_sj2 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start2:$end2"};
                my $alt_sj = "$alt_sj1;$alt_sj2";

                my $model_avg_sj = $model_sj;
                my $alt_avg_sj   = &avgsj($alt_sj);

                if ( $alt_avg_sj + $model_avg_sj == 0 ) {
                    $as{$key}->{ $name_group[$i] } = "-\t-\t-\t-\t-\t-";
                    next;
                }
                my $alt_ratio = sprintf( "%.2f",
                    $alt_avg_sj / ( $alt_avg_sj + $model_avg_sj ) );

                $as{$key}->{ $name_group[$i] } =
                    $alt_ratio . "\t"
                  . $alt_avg_sj . "\t"
                  . $model_avg_sj . "\t"
                  . $alt_sj . "\t"
                  . $model_sj . "\t"
                  . $filterflag;
            }

# Chr01   351716,350505   352170,351664   -       351992,350505   352170,351928   MXE     0.5     100     100     PAC:26823778
            elsif ( $type eq "MXE" ) {
                my ( $start1, $start2 ) = split( /;/, $line[1] );
                my ( $end1,   $end2 )   = split( /;/, $line[2] );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start1:$end1"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start1:$end1"} );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start2:$end2"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start2:$end2"} );
                my $alt_sj1 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start1:$end1"};
                my $alt_sj2 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start2:$end2"};
                my $alt_sj = "$alt_sj1;$alt_sj2";

                ( $start1, $start2 ) = split( /;/, $line[4] );
                ( $end1,   $end2 )   = split( /;/, $line[5] );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start1:$end1"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start1:$end1"} );
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start2:$end2"}
                  = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$start2:$end2"} );
                my $model_sj1 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start1:$end1"};
                my $model_sj2 = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$start2:$end2"};
                my $model_sj = "$model_sj1;$model_sj2";

                my $model_avg_sj = &avgsj($model_sj);
                my $alt_avg_sj   = &avgsj($alt_sj);

                if ( $alt_avg_sj + $model_avg_sj == 0 ) {
                    $as{$key}->{ $name_group[$i] } = "-\t-\t-\t-\t-\t-";
                    next;
                }
                my $alt_ratio = sprintf( "%.2f",
                    $alt_avg_sj / ( $alt_avg_sj + $model_avg_sj ) );

                $as{$key}->{ $name_group[$i] } =
                    $alt_ratio . "\t"
                  . $alt_avg_sj . "\t"
                  . $model_avg_sj . "\t"
                  . $alt_sj . "\t"
                  . $model_sj . "\t"
                  . $filterflag;
            }
            else {
                $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[1]:$line[2]"} = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$line[1]:$line[2]"} );
                my $alt_sj = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[1]:$line[2]"};
                $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[4]:$line[5]"} = 0
                  if not defined( $junction{ $name_group[$i] }{$chr}{$strand}
                      ->{"$line[4]:$line[5]"} );
                my $model_sj = $junction{ $name_group[$i] }{$chr}{$strand}
                  ->{"$line[4]:$line[5]"};

                my $alt_avg_sj   = $alt_sj;
                my $model_avg_sj = $model_sj;

                if ( $alt_avg_sj + $model_avg_sj == 0 ) {
                    $as{$key}->{ $name_group[$i] } = "-\t-\t-\t-\t-\t-";
                    next;
                }
                my $alt_ratio = sprintf( "%.2f",
                    $alt_avg_sj / ( $alt_avg_sj + $model_avg_sj ) );

                $as{$key}->{ $name_group[$i] } =
                    $alt_ratio . "\t"
                  . $alt_avg_sj . "\t"
                  . $model_avg_sj . "\t"
                  . $alt_sj . "\t"
                  . $model_sj . "\t"
                  . $filterflag;
            }

        }
    }

    my $count=1;
    foreach my $key ( sort keys %as ) {
        my $id = $label."_NIR_".$count;
        $count++;
        print OUT $id,"\t",$key;
        for ( my $i = 0 ; $i <= $#name_group ; $i++ ) {
            if ( not defined( $as{$key}->{ $name_group[$i] } ) ) {
                $as{$key}->{ $name_group[$i] } = "-\t-\t-\t-\t-\t-";
            }
            print OUT "\t$as{$key}->{$name_group[$i]}";
        }
        print OUT "\n";
    }

}

###子过程-----load_all_sj($sj_file)
sub load_all_sj {
    for ( my $i = 0 ; $i <= $#sj_group ; $i++ ) {
        open SJ, $sj_group[$i] || die;
        while (<SJ>) {    #把所有的junction添加进hash
                #默认为读取通过tophat获取的junctions.bed文件
             #chr20   199843  204678  JUNC00000001    24      -       199843  204678  255,0,0 2       65,70   0,4765
             #chr20   275722  276146  JUNC00000015    3       +       275722  276146  255,0,0 2       64,69   0,355
            chomp;
            next if ( $_ =~ /^track/i );
            my $start = 0;
            my $end   = 0;
            my (
                $chr,      $left,       $right,      $feature,
                $readsNum, $strand,     $thickStart, $thickEnd,
                $itemRgb,  $blockCount, $blockSizes, $blockStarts
            ) = split /\s+/;
            my ( $leftsize, $rightsize ) =
              ( split /\,/, $blockSizes )[ 0 .. 1 ];
            $start = $left + $leftsize;
            $end   = $right - $rightsize + 1;

            if (
                defined(
                    $junction{ $name_group[$i] }{$chr}{$strand}->{"$start:$end"}
                )
              )
            {
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start:$end"} +=
                  $readsNum;
            }
            else {
                $junction{ $name_group[$i] }{$chr}{$strand}->{"$start:$end"} =
                  $readsNum;
            }
        }
        close SJ;

        #调试信息
        print "done reading SJ ...........", "\n\n";
    }
}

sub sumsj {
    my ($sjlist) = @_;
    my @num = split( /,/, $sjlist );
    return sum(@num);

}

sub avgsj {
    my ($sjlist) = @_;
    my @num = split( /;/, $sjlist );
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
