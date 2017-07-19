#!/usr/bin/perl -w
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;
use POSIX ":sys_wait_h";
use Bio::DB::Sam;

my %opts;
GetOptions( \%opts, "bam=s", "fa=s", "chrlen=s", "anno=s", "od=s", "h" );

if (   !defined( $opts{bam} )
    || !defined( $opts{chrlen} )
    || !defined( $opts{fa} )
    || !defined( $opts{anno} )
    || !defined( $opts{od} )
    || defined( $opts{h} ) )
{
    print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-bam         splice junction        must be given;

		-chrlen      chromsome lenght       must be given;

		-fa          genome fa file         must be given;

		-anno        anno file              must be given;

		-od          out directory          must be given;

	Usage End.

    exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

my $bam_file    = $opts{bam};
my $fa_file     = $opts{fa};
my $chrlen_file = $opts{chrlen};
my $anno_file   = $opts{anno};
my $outdir      = $opts{od};

my %chrlen = ();
open CHRLEN, $chrlen_file || die;
while (<CHRLEN>) {

    #chr1	201139347
    chomp;
    my ( $chr, $len ) = split /\s+/;
    $chrlen{$chr} = $len;

}
close CHRLEN;

my $outfile = $outdir . "/" . "EXP_exon";
open EXON, ">$outfile" || die;
print EXON "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
$outfile = $outdir . "/" . "EXP_intron";
open INTRON, ">$outfile" || die;
print INTRON "#Chromesome\tFeature\tStrand\tStart\tEnd\tGene\tBaseNumber\n";
close EXON;
close INTRON;

## == number of proc ==
my $num_proc = 0;

## == number of collected ==
my $num_collect = 0;

my $collect;

## == get the child signal ==
$SIG{CHLD} = sub { $num_proc-- };

my $i = 0;
my @pids;
foreach my $chr ( sort keys %chrlen ) {
    $i++;
    ## == fork a new process ==
    my $pid = fork();

    if ( !defined($pid) ) {
        print "Error in fork: $!";
        exit 1;
    }

    if ( $pid == 0 ) {
        ## == child proc ==
        print "begin $chr\n";
        push @pids, $pid;
        my $outfile = $outdir . "/" . "EXP_exon";
        open EXON, ">$outfile\_$chr" || die;
        $outfile = $outdir . "/" . "EXP_intron";
        open INTRON, ">$outfile\_$chr" || die;
        my $bam = Bio::DB::Sam->new(
            -fasta => "$fa_file",
            -bam   => "$bam_file"
        );

        open ANNO, $anno_file || die;

        #NC_008467.1     mRNA       +       36663   37331   POPTRDRAFT_750402

        while (<ANNO>) {
            chomp;
            my ( $thischr, $feature, $strand, $start, $end, $gene ) = split;
            next if $thischr ne $chr;
            next if $feature ne "exon" && $feature ne "intron";

            my %base_p = ();
            my %base_n = ();
            &readbam( $bam, $thischr, $start, $end, \%base_p, \%base_n );
            my $base_num =
              &sumbase( $thischr, $start, $end, $strand, \%base_p, \%base_n );

            print EXON
              "$thischr\t$feature\t$strand\t$start\t$end\t$gene\t$base_num\n"
              if $feature eq "exon";
            print INTRON
              "$thischr\t$feature\t$strand\t$start\t$end\t$gene\t$base_num\n"
              if $feature eq "intron";
        }

        close ANNO;
        close EXON;
        close INTRON;

        print "done $chr\n";
        exit 0;
    }

    $num_proc++;
    ## == if need to collect zombies ==
    if ( ( $i - $num_proc - $num_collect ) > 0 ) {
        while ( ( $collect = waitpid( -1, WNOHANG ) ) > 0 ) {
            $num_collect++;
            print $num_collect, "\n";
        }
    }

    do {
        sleep(1);
    } until ( $num_proc < 25 );

}

## 等待所有子进程全部完成
while ( ( $i - $num_collect ) > 0 ) {
    sleep(1);
    while ( ( $collect = waitpid( -1, WNOHANG ) ) > 0 ) {
        $num_collect++;
        print $num_collect, "\n";
    }
}

`cat $outdir/EXP_exon_* >> $outdir/EXP_exon && rm -rf $outdir/EXP_exon_*`;
`cat $outdir/EXP_intron_* >> $outdir/EXP_intron && rm -rf $outdir/EXP_intron_*`;
print "all done\n";

sub readbam {
    my ( $bam, $chr, $start, $end, $base_p, $base_n ) = @_;
    my $left_block   = 0;
    my $middle_block = 0;
    my $right_block  = 0;
    my $len          = 0;
    my @alignments   = $bam->get_features_by_location(
        -seq_id => $chr,
        -start  => $start,
        -end    => $end
    );
    for my $a (@alignments) {
        my $seqid  = $a->seq_id;
        my $start  = $a->start;
        my $end    = $a->end;
        my $strand = $a->strand;
        my $cigar  = $a->cigar_str;
        next if ( $cigar !~ /^\d+M$/ && $cigar !~ /^\d+M\d+N\d+M/ );

        if ( $strand == 1 ) {    #正向
            if ( $cigar =~ m/^(\d+)M$/ ) {
                $len = $1;
                ##记录base数
                for ( my $i = $start ; $i < $start + $len ; $i++ ) {
                    $base_p->{$chr}->{$i} = 0
                      if not defined( $base_p->{$chr}->{$i} );
                    $base_p->{$chr}->{$i}++;
                }
            }
            elsif ( $cigar =~ /(\d+)M(\d+)N(\d+)M/ ) {
                $left_block   = $1;
                $middle_block = $2;
                $right_block  = $3;
                for ( my $i = $start ; $i < $start + $left_block ; $i++ ) {
                    $base_p->{$chr}->{$i} = 0
                      if not defined( $base_p->{$chr}->{$i} );
                    $base_p->{$chr}->{$i}++;
                }
                for (
                    my $i = $start + $left_block + $middle_block ;
                    $i < $start + $left_block + $middle_block + $right_block ;
                    $i++
                  )
                {
                    $base_p->{$chr}->{$i} = 0
                      if not defined( $base_p->{$chr}->{$i} );
                    $base_p->{$chr}->{$i}++;
                }

            }
        }
        else {    #负向
            if ( $cigar =~ m/^(\d+)M$/ ) {
                $len = $1;
                ##记录base数
                for ( my $i = $start ; $i < $start + $len ; $i++ ) {
                    $base_n->{$chr}->{$i} = 0
                      if not defined( $base_n->{$chr}->{$i} );
                    $base_n->{$chr}->{$i}++;
                }
            }
            elsif ( $cigar =~ /(\d+)M(\d+)N(\d+)M/ ) {
                $left_block   = $1;
                $middle_block = $2;
                $right_block  = $3;
                for ( my $i = $start ; $i < $start + $left_block ; $i++ ) {
                    $base_n->{$chr}->{$i} = 0
                      if not defined( $base_n->{$chr}->{$i} );
                    $base_n->{$chr}->{$i}++;
                }
                for (
                    my $i = $start + $left_block + $middle_block ;
                    $i < $start + $left_block + $middle_block + $right_block ;
                    $i++
                  )
                {
                    $base_n->{$chr}->{$i} = 0
                      if not defined( $base_n->{$chr}->{$i} );
                    $base_n->{$chr}->{$i}++;
                }
            }
        }
    }
}

#####方法：取base总和
sub sumbase {
    my ( $chr, $start, $end, $strand, $base_p, $base_n ) = @_;
    my $sum_value = 0;
    if ( $strand eq "+" ) {
        for ( my $j = $start ; $j <= $end ; $j++ ) {
            $base_p->{$chr}->{$j} = 0 if not defined( $base_p->{$chr}->{$j} );
            $sum_value += $base_p->{$chr}->{$j};
        }
    }
    if ( $strand eq "-" ) {
        for ( my $j = $start ; $j <= $end ; $j++ ) {
            $base_n->{$chr}->{$j} = 0 if not defined( $base_n->{$chr}->{$j} );
            $sum_value += $base_n->{$chr}->{$j};
        }
    }

    return $sum_value;
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
