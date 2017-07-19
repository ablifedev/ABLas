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
GetOptions( \%opts, "bam=s", "fa=s", "chrlen=s", "gff=s", "mod=s", "o=s", "h" );

if (   !defined( $opts{bam} )
    || !defined( $opts{chrlen} )
    || !defined( $opts{fa} )
    || !defined( $opts{gff} )
    || !defined( $opts{mod} )
    || !defined( $opts{o} )
    || defined( $opts{h} ) )
{
    print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

		-bam         splice junction        must be given;

		-chrlen      chromsome lenght       must be given;

		-fa          genome fa file         must be given;

		-gff         gff file               must be given;

		-mod         gene model file        must be given;

		-o           out file               must be given;

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
my $gff_file    = $opts{gff};
my $mod_file    = $opts{mod};
my $out_file    = $opts{o};

my %chrlen = ();
open CHRLEN, $chrlen_file || die;
while (<CHRLEN>) {

    #chr1	201139347
    chomp;
    my ( $chr, $len ) = split /\s+/;
    $chrlen{$chr} = $len;

}
close CHRLEN;

#将gene model保存进哈希以便查询
my %gene_model = ();
open MOD, $mod_file || die;
while (<MOD>) {
    chomp;
    next if ( $_ =~ /chrom/ );
    my @line = split;
    $gene_model{ $line[1] }{$line[0]}=1;
}
close MOD;

open OUT, ">$out_file" || die;
print OUT
"#chr\tss1\tss2\tstrand\tIntronR\tintron_ave_base/ends_exon_ave_base\tintron_ave_base\tends_exon_ave_base\tforward_isoform\n";
close OUT;

open LOG, ">$out_file\.log" || die;
print LOG
"#chr\tss1\tss2\tstrand\tIntronR\tintron_ave_base/ends_exon_ave_base\tintron_ave_base\tends_exon_ave_base\tforward_isoform\n";
close LOG;

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

        my $bam = Bio::DB::Sam->new(
            -fasta => "$fa_file",
            -bam   => "$bam_file"
        );

        ###遍历gff，按基因检测AS事件
        my (
            $gene,          $isoform,         $forward_chr,
            $forward_gene,  $forward_isoform, $gene_count,
            $forward_start, $forward_end,     $exon_region_value,
            $intron_region_value
        ) = ( "", "", "", "", "", 1, 0, 0, 1, 0 );
        my $forward_strand = "";
        my %exon           = ();
        my @exon_hash      = ();    #基因的exon区域
        my @intron_hash    = ();    #基因的intron区域
        my %sj_start       = ();    #基因中所有起点的sj信息
        my %sj_end         = ();    #基因中所有终点的sj信息
        my %gene_region    = ();
        my $gene_start     = 0;
        my $gene_end       = 0;
        my $gene_strand    = "";
        my %as_count = ( "IntronR" => 0 );
        my $exon_number = 0;

#chr1    hg18_knownGene  mRNA    1310954 1323585 .       -       .       ID=CCNL2.1.5;Parent=CCNL2.1
#chr1    hg18_knownGene  exon    1323476 1323585 0.000000        -       .       ID=exon:CCNL2.1.5:1;Parent=CCNL2.1.5
#chr1    hg18_knownGene  intron  1320781 1323475 0.000000        -       .       ID=intron:CCNL2.1.5:2;Parent=CCNL2.1.5
#chr1    hg18_knownGene  exon    1320637 1320780 0.000000        -       .       ID=exon:CCNL2.1.5:2;Parent=CCNL2.1.5
        open OUT, ">$out_file\_$chr" || die;
		open LOG, ">$out_file\.log_$chr" || die;

        open GFF, $gff_file || die;
        while (<GFF>) {
            chomp;
            next if ( $_ =~ /^\#/ );
            my (
                $thischr, $source, $feature, $start, $end,
                $other1,  $strand, $other2,  $info
            ) = split /\t/;
            next if $thischr ne $chr;
            next
              if ( $feature =~ /chromosome/
                || $feature =~ /gene/
                || $feature =~ /protein/ );
            if ( $feature =~ /^mRNA|transcript$/ ) {
                $info =~ m/ID=([\w\:\.\-]+);.*Parent=([\w\:\.\-]+)/;
                $isoform = "$1";
                $gene    = "$2";

            }
            if ( $feature =~ /^exon$/ ) {
                $info =~ m/Parent=([\w\:\.\-]+)/;
                $isoform = "$1";
            }
            next if (  not defined( $gene_model{$gene}{$isoform} ) );

            if (   $isoform ne $forward_isoform
                && $forward_isoform ne ""
                && $exon_number > 1 ) {
                ##处理基因方法，检测AS事件
                &handling_gene(
                    $bam,            \%exon,        \@exon_hash,
                    \@intron_hash,   \%gene_region, $forward_isoform,
                    $forward_strand, $forward_chr,  $gene_start,
                    $gene_end,       \*OUT, \*LOG
                );
                %exon        = ();
                %gene_region = ();
                %sj_start    = ();
                %sj_end      = ();
                @exon_hash   = ();
                @intron_hash = ();
                $gene_count++;
                $exon_region_value   = 1;
                $intron_region_value = 0;
                $exon_number = 0;
            }
            if ( $feature =~ /^mRNA|transcript$/ ) {
                $gene_start  = $start;
                $gene_end    = $end;
                $gene_strand = $strand;
                %exon        = ();
                %gene_region = ();
                %sj_start    = ();
                %sj_end      = ();
                @exon_hash   = ();
                @intron_hash = ();
                $gene_count++;
                $exon_region_value   = 1;
                $intron_region_value = 0;
                $exon_number = 0;
            }
            if ( $feature =~ /^exon$/ ) {
                $exon_number++;
                push @exon_hash, "$start:$end";
                if ( $exon_region_value > 1 ) {    #给intron区域赋值
                    my $temp_start = 0;
                    my $temp_end   = 0;
                    if ( $strand eq "+" ) {
                        $temp_start            = $forward_end;
                        $temp_end              = $start;
                        $sj_start{$temp_start} = $temp_end;
                        $sj_end{$temp_end}     = $temp_start;
                    }
                    if ( $strand eq "-" ) {
                        $temp_start          = $end;
                        $temp_end            = $forward_start;
                        $sj_start{$temp_end} = $temp_start;
                        $sj_end{$temp_start} = $temp_end;
                    }
                    push @intron_hash, "$temp_start:$temp_end";

                    $gene_region{$intron_region_value} =
                      "$temp_start:$temp_end";
                }
                $gene_region{$exon_region_value} = "$start:$end";

                $exon_region_value += 1;
                $intron_region_value -= 1;
                $forward_end     = $end;
                $forward_start   = $start;
                $forward_isoform = $isoform;
                $forward_strand  = $strand;
                $exon{$start} =
                  $end;    #将exon起始终止位点记录进入hash %exon
            }
            $forward_gene = $gene;
            $forward_chr  = $thischr;
        }
        close GFF;

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

print "all done\n";

`cat $out_file\_* >> $out_file && rm -rf $out_file\_*`;
`cat $out_file\.log\_* >> $out_file\.log && rm -rf $out_file\.log\_*`;

# ###
# print "Done detect IntronR...\n";
# ###
# ###END detect IntronR

sub handling_gene {
    my (
        $bam,         $exon,            $exon_hash, $intron_hash,
        $gene_region, $forward_isoform, $strand,    $chr,
        $gene_start,  $gene_end,        $outfile_handle, $log_handle
    ) = @_;
    my %base_p = ();
    my %base_n = ();

    # print $forward_isoform, "\n";
    &readbam( $bam, $chr, $gene_start, $gene_end, \%base_p, \%base_n );
    my $sum_intron_len  = 0;
    my $sum_intron_base = 0;
    my $ave_base_all    = 0;
    my $intron_num      = 0;

    foreach my $index ( sort { $a <=> $b } keys %{$gene_region} ) {
        next if $index > 0;
        my ( $intron_start, $intron_end ) = split /\:/, $gene_region->{$index};
        $sum_intron_len += $intron_end - $intron_start + 1;
        $sum_intron_base +=
          &sumbase( $chr, $intron_start, $intron_end, $strand, \%base_p,
            \%base_n );
        $intron_num++;
    }
    $ave_base_all = $sum_intron_base / $sum_intron_len if $sum_intron_len != 0;

    ##
    # print "\n",$chr,"\t",$strand,"\n";

  LABEL1: foreach my $index ( sort { $a <=> $b } keys %{$gene_region} ) {
        next if $index > 0;
        my $ends_exon_len      = 0;
        my $ends_exon_base     = 0;
        my $ends_exon_ave_base = 0;

        my ( $intron_start, $intron_end ) = split /\:/, $gene_region->{$index};

        # print $intron_start,":",$intron_end,"\n";

        ##如果有其他as事件则跳过。
# next if(defined($as{$chr}{$intron_start}) || defined($as{$chr}{$intron_end}));

        my $intron_len = $intron_end - $intron_start - 1;

        # print "========",$intron_len,"\n";
        next if $intron_len <= 0;
        my $intron_base = &sumbase(
            $chr,
            $intron_start + 1,
            $intron_end - 1,
            $strand, \%base_p, \%base_n
        );
        my $intron_ave_base = 0;
        $intron_ave_base = $intron_base / $intron_len if $intron_len > 0;

        ##筛选条件3：intron区base平均depth要大于左右两边exon区base平均depth的20%.
        my ( $left_exon_start, $left_exon_end ) = split /\:/,
          $gene_region->{ 0 - $index };
        $ends_exon_len += $left_exon_end - $left_exon_start + 1;
        $ends_exon_base +=
          &sumbase( $chr, $left_exon_start, $left_exon_end, $strand, \%base_p,
            \%base_n );

# print "exon1: ",$left_exon_start,"\t",$left_exon_end,"\t",$ends_exon_base,"\n";
        my ( $right_exon_start, $right_exon_end ) = split /\:/,
          $gene_region->{ 1 - $index };
        $ends_exon_len += $right_exon_end - $right_exon_start + 1;
        $ends_exon_base +=
          &sumbase( $chr, $right_exon_start, $right_exon_end, $strand,
            \%base_p, \%base_n );
        $ends_exon_ave_base = $ends_exon_base / $ends_exon_len;

# print "exon2: ",$right_exon_start,"\t",$right_exon_end,"\t",$ends_exon_base,"\n";
# print "exon3: ",$ends_exon_len,"\t",$ends_exon_base,"\t",$ends_exon_ave_base,"\n";
        my $ratio = 1;
        $ratio = $intron_ave_base / $ends_exon_ave_base
          if $ends_exon_ave_base != 0;
        $ratio              = sprintf( "%.2f", $ratio );
        $ends_exon_ave_base = sprintf( "%.2f", $ends_exon_ave_base );
        $intron_ave_base    = sprintf( "%.2f", $intron_ave_base );
        print $log_handle
"$chr\t$intron_start\t$intron_end\t$strand\t-\t$ratio\t$intron_ave_base\t$ends_exon_ave_base\t$forward_isoform\n";
        next LABEL1 if $intron_ave_base < $ends_exon_ave_base * 0.2;

        # print "$ends_exon_ave_base\n";

        # print $intron_base,"\t",$intron_ave_base,"\n";
        next LABEL1 if $intron_base < 150;

# next LABEL1 if $intron_ave_base < $ave_base_all * 2;		##筛选条件1：intron区base平均depth要大于总intron区base平均depth.
        next LABEL1 if $intron_ave_base < $ave_base_all * 2 && $intron_num >= 2;

# print $ave_base_all,"\t",$intron_ave_base,"\t",$intron_base,"\t",$intron_len,"\n";

        ##筛选条件2：左右两边border reads至少有一边大于1。
        my $dd = 0;
        for ( my $i = $intron_start ; $i <= $intron_start + 4 ; $i++ ) {

            # next LABEL1 if &sumbase($chr, $i, $i)<1;
            if ( &sumbase( $chr, $i, $i, $strand, \%base_p, \%base_n ) < 1 ) {
                $dd++;
                last;
            }
        }
        for ( my $i = $intron_end - 4 ; $i <= $intron_end ; $i++ ) {

            # next LABEL1 if &sumbase($chr, $i, $i)<1;
            if ( &sumbase( $chr, $i, $i, $strand, \%base_p, \%base_n ) < 1 ) {
                $dd++;
                last;
            }
        }
        next LABEL1 if $dd == 2;

        # print "have border reads\n";

        ##输出intron
#chr11   AS      ES      2923452 2929065 0.66    -       .       ID=ES_2;gene=NAP1L4.1;mRNA=NAP1L4.1.1;
# my $idnum = 0;
# $as_count{"IntronR"}++;
# $idnum = $as_count{"IntronR"};
# my $asid = "IntronR\_$idnum";

# print OUT "$chr\tAS\tIntronR\t$intron_start\t$intron_end\t\.\t$strand\t\.\tID=$asid\;gene=$forward_gene\;mRNA=$forward_isoform\;\n";
# print OUT "$chr\tAS\texon\t$left_exon_start\t$left_exon_end\t\.\t$strand\t\.\tParent=$asid\n";
# print OUT "$chr\tAS\texon\t$right_exon_start\t$right_exon_end\t\.\t$strand\t\.\tParent=$asid\n";

        print $outfile_handle
"$chr\t$intron_start\t$intron_end\t$strand\tIntronR\t$ratio\t$intron_ave_base\t$ends_exon_ave_base\t$forward_isoform\n";
    }
}

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
