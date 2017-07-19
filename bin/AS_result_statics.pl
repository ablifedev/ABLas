#!/usr/bin/perl -w
# 
# Copyright (c)   AB_Life 2011
# Writer:         dongchen <dongchen@ablife.cc>
# Program Date:   2011.06.28
# Modifier:       dongchen <dongchen@ablife.cc>
# Last Modified:  2011.06.28
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your programmeyou must write the detailed timediscriptionsparameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions(\%opts,"indir=s","knir=s","kir=s","od=s");

if(!defined($opts{indir}) || !defined($opts{od}) )
{
	print <<"	Usage End.";

	Description:This programme is used for 
		
		Version: $ver

	Usage:perl $0

		-indir           in directory          must be given

		-knir            nir file              must be given

		-kir             ir file               must be given

		-od              report result dir     default is ./

	Usage End.

	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################

#Write your program
my $indir = $opts{indir};
my $nir_file = $opts{knir};
my $ir_file = $opts{kir};
my $outdir = $opts{od} || "./";

chdir($outdir);

my %Config = ();
&Load_config_file($nir_file,$ir_file,\%Config);

my @AS_type = ("3pMXE","5pMXE","A3SS","A3SS&ES","A5SS","A5SS&ES","cassetteExon","ES","MXE","IntronR");
my %Result = ();
my @Event = ();
my $total = 0;

&Exon_exp($indir,\%Result);

&SJ_stat($indir,\%Result);

&AS_stat($indir,\%Result);

&RAS_stat($indir,\%Result);


my $outfile = "Exon_detection_results_in_RNA-seq_data.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "#Sample\tDetected_exon\tAnnotatedExon\tRatio\n";
foreach my $sam (sort keys %{$Result{"exon_exp"}}) {
	next if($sam eq "total" || $sam eq "TotalNum");
	print OUT $sam,"\t",$Result{"exon_exp"}{$sam},"\t",$Result{"exon_exp"}{"TotalNum"},"\t",sprintf("%.2f",$Result{"exon_exp"}{$sam}/$Result{"exon_exp"}{"TotalNum"}*100),"%\n";
}
print OUT "Total\t",$Result{"exon_exp"}{"total"},"\t",$Result{"exon_exp"}{"TotalNum"},"\t",sprintf("%.2f",$Result{"exon_exp"}{"total"}/$Result{"exon_exp"}{"TotalNum"}*100),"%\n";
close OUT;

$outfile = "Splicing_Junction_analysis_of_samples_from_RNA-seq_data.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "#Sample\tAll_Junction\tNovel_Junction\tKnown_Junction\n";
foreach my $sam (sort keys %{$Result{"SJ"}}) {
	next if($sam eq "Total");
	print OUT $sam,"\t",$Result{"SJ"}{$sam}{"all"},"\t",$Result{"SJ"}{$sam}{"novel"},"\t",$Result{"SJ"}{$sam}{"known"},"\n";
}
print OUT "Total\t",$Result{"SJ"}{"Total"}{"all"},"\t",$Result{"SJ"}{"Total"}{"novel"},"\t",$Result{"SJ"}{"Total"}{"known"},"\n";
close OUT;

$outfile = "All_known_AS_events_detected_from_all_samples.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "#Sample\tknownAS(KAS)\tAllAS\tKAS%\n";
my $all_known = 0;
my $all_all = 0;
foreach my $sam (sort keys %{$Result{"AS"}}) {
	
	if ($sam eq "All") {
		foreach my $event (sort @AS_type) {
			$all_known += $Result{"AS"}{$sam}{"All_known"}{$event};
			$all_all += $Result{"AS"}{$sam}{"All_total"}{$event};
		}
		next;
	}
	print OUT $sam,"\t";
	my $sum_k = 0;
	my $sum_a = 0;
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"known_result"}{$event}) {
			$sum_k += $Result{"AS"}{$sam}{"known_result"}{$event};
			$sum_a += $Result{"AS"}{$sam}{"all_result"}{$event};
		} else {
			$sum_a += $Result{"AS"}{$sam}{"all_result"}{$event};
			next;
		}
	}
	print OUT $sum_k,"\t",$sum_a,"\t",sprintf("%.2f",$sum_k/$sum_a*100),"%\n";
}
print OUT "Total\t",$all_known,"\t",$all_all,"\t",sprintf("%.2f",$all_known/$all_all*100),"%\n";
close OUT;


$outfile = "All_known_AS_events_detected_by_Type.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "#Type\tAnnotated\tDetected\tRatio\n";
my $sum_k_bytype = 0;
my $sum_k_d_bytype = 0;
foreach my $event(sort @AS_type) {
	#print $event,"\n";
	print OUT $event,"\t",$Config{"AS_type"}{$event},"\t",$Result{"AS"}{"All"}{"All_known"}{$event},"\t",sprintf("%.2f",$Result{"AS"}{"All"}{"All_known"}{$event}/$Config{"AS_type"}{$event}*100),"%\n" if (defined($Config{"AS_type"}{$event}) && $Config{"AS_type"}{$event}>0);
	$sum_k_bytype += $Config{"AS_type"}{$event};
	$sum_k_d_bytype += $Result{"AS"}{"All"}{"All_known"}{$event};
}
print OUT "Total\t",$sum_k_bytype,"\t",$sum_k_d_bytype,"\t",sprintf("%.2f",$sum_k_d_bytype/$sum_k_bytype*100),"%\n";
close OUT;


$outfile = "Novel_AS_events_detected_from_all_samples.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "#Sample\tnovel AS(NAS)\tAllAS\tNAS%\n";
my $all_novel = 0;
foreach my $sam (sort keys %{$Result{"AS"}}) {
	if ($sam eq "All") {
		foreach my $event (sort @AS_type) {
			if(exists $Result{"AS"}{$sam}{"All_novel"}{$event}) {
				$all_novel += $Result{"AS"}{$sam}{"All_novel"}{$event};
			}
		}
		next;
	}
	print OUT $sam,"\t";
	my $sum_k = 0;
	my $sum_a = 0;
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"novel_result"}{$event}) {
			$sum_k += $Result{"AS"}{$sam}{"novel_result"}{$event};
			$sum_a += $Result{"AS"}{$sam}{"all_result"}{$event};
		} else {
			$sum_a += $Result{"AS"}{$sam}{"all_result"}{$event};
			next;
		}
	}
	print OUT $sum_k,"\t",$sum_a,"\t",sprintf("%.2f",$sum_k/$sum_a*100),"%\n";
}
print OUT "Total\t",$all_novel,"\t",$all_all,"\t",sprintf("%.2f",$all_novel/$all_all*100),"%\n";
close OUT;


$outfile = "Classification_of_all_the_detected_AS_events.xls";
open(OUT,">$outfile") || die "$!\n";
#print OUT "Sample\t",join("\t", sort keys %{$Result{"AS"}{"All"}{"All_total"}}),"\n";
print OUT "#Sample\t",join("\t",sort @AS_type),"\tTotal\tDetected junction\tAS/100SJ\n";
foreach my $sam (sort keys %{$Result{"AS"}}) {
	next if ($sam eq "All") ;
	print OUT $sam;
	my $total = 0;
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"all_result"}{$event}) {
			print OUT "\t",$Result{"AS"}{$sam}{"all_result"}{$event};
			$total += $Result{"AS"}{$sam}{"all_result"}{$event};
		} else {
			print OUT "\t",0;
		}
	}
	# print "all:\t",$sam,"\t",$total,"\n";
	print OUT "\t",$total,"\t",$Result{"SJ"}{$sam}{"all"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{$sam}{"all"});
	print OUT "\n";
}
print OUT "Total";
$total = 0;
foreach my $event (sort @AS_type) {
	if(exists $Result{"AS"}{"All"}{"All_total"}{$event}) {
		print OUT "\t",$Result{"AS"}{"All"}{"All_total"}{$event};
		$total += $Result{"AS"}{"All"}{"All_total"}{$event};
	} else {
		print OUT "\t0";
	}
}
print OUT "\t",$total,"\t",$Result{"SJ"}{"Total"}{"all"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{"Total"}{"all"});
close OUT;


$outfile = "Classification_of_all_the_detected_known_AS_events.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "Sample\t",join("\t",sort @AS_type),"\tTotal\tDetected junction\tAS/100SJ\n";
#print OUT "Sample\t",join("\t", sort keys %{$Result{"AS"}{"All"}{"All_total"}}),"\n";
foreach my $sam (sort keys %{$Result{"AS"}}) {
	next if ($sam eq "All") ;
	print OUT $sam;
	my $total = 0;
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"known_result"}{$event}) {
			print OUT "\t",$Result{"AS"}{$sam}{"known_result"}{$event};
			$total += $Result{"AS"}{$sam}{"known_result"}{$event};
		} else {
			print OUT "\t",0;
		}
	}
	print OUT "\t",$total,"\t",$Result{"SJ"}{$sam}{"known"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{$sam}{"known"});
	print OUT "\n";
}
print OUT "Total";
$total = 0;
foreach my $event (sort @AS_type) {
	if(exists $Result{"AS"}{"All"}{"All_known"}{$event}) {
		print OUT "\t",$Result{"AS"}{"All"}{"All_known"}{$event};
		$total += $Result{"AS"}{"All"}{"All_known"}{$event};
	} else {
		print OUT "\t0";
	}
}
print OUT "\t",$total,"\t",$Result{"SJ"}{"Total"}{"known"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{"Total"}{"known"});
close OUT;


$outfile = "Classification_of_all_the_detected_novel_AS_events.xls";
open(OUT,">$outfile") || die "$!\n";
# print OUT "Classification of all the detected novel AS events\n";
print OUT "Sample\t",join("\t",sort @AS_type),"\tTotal\tDetected junction\tAS/100SJ\n";
#print OUT "Sample\t",join("\t", sort keys %{$Result{"AS"}{"All"}{"All_total"}}),"\n";
foreach my $sam (sort keys %{$Result{"AS"}}) {
	next if ($sam eq "All") ;
	print OUT $sam;
	my $total = 0;
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"novel_result"}{$event}) {
			print OUT "\t",$Result{"AS"}{$sam}{"novel_result"}{$event};
			$total += $Result{"AS"}{$sam}{"novel_result"}{$event};
		} else {
			print OUT "\t",0;
		}
	}
	print OUT "\t",$total,"\t",$Result{"SJ"}{$sam}{"novel"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{$sam}{"novel"});
	print OUT "\n";
}
print OUT "Total";
$total = 0;
foreach my $event (sort @AS_type) {
	if(exists $Result{"AS"}{"All"}{"All_novel"}{$event}) {
		print OUT "\t",$Result{"AS"}{"All"}{"All_novel"}{$event};
		$total += $Result{"AS"}{"All"}{"All_novel"}{$event};
	} else {
		print OUT "\t0";
	}
}
print OUT "\t",$total,"\t",$Result{"SJ"}{"Total"}{"novel"},"\t",sprintf("%.2f",$total*100/$Result{"SJ"}{"Total"}{"novel"});
close OUT;


$outfile = "Distribution_of_each_class_of_All_AS_events.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "Sample\t",join("\t",@AS_type),"\tTotal\n";
#print OUT "Sample\t",join("\t", sort keys %{$Result{"AS"}{"All"}{"All_total"}}),"\tTotal\n";
foreach my $sam (sort keys %{$Result{"AS"}}) {
	next if ($sam eq "All") ;
	@Event = ();
	$total = 0;
	print OUT $sam;
	# foreach my $events (sort keys %{$Result{"AS"}{$sam}{"all_result"}}) {
	# 	push @Event,$Result{"AS"}{$sam}{"all_result"}{$events};
	# 	$total += $Result{"AS"}{$sam}{"all_result"}{$events};
	# }
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"all_result"}{$event}) {
			push @Event,$Result{"AS"}{$sam}{"all_result"}{$event};
			$total += $Result{"AS"}{$sam}{"all_result"}{$event};
		}
	}
	foreach my $event (sort @AS_type) {
		if(exists $Result{"AS"}{$sam}{"all_result"}{$event}) {
			print OUT "\t",sprintf("%.2f",$Result{"AS"}{$sam}{"all_result"}{$event}/$total*100),"%";
		} else {
			print OUT "\t0%";
		}
	}
	#foreach my $events (@Event) {
		#print OUT "\t",sprintf("%.2f",$events/$total*100),"%";
	#}
	print OUT "\t",sprintf("%.2f",$total/$total*100),"%";
	print OUT "\n";
}
print OUT "Total";
@Event = ();
$total = 0;
# foreach my $events (sort keys %{$Result{"AS"}{"All"}{"All_total"}}) {
# 	push @Event,$Result{"AS"}{"All"}{"All_total"}{$events};
# 	$total += $Result{"AS"}{"All"}{"All_total"}{$events};
# }
foreach my $event (sort @AS_type) {
	if(exists $Result{"AS"}{"All"}{"All_total"}{$event}) {
		push @Event,$Result{"AS"}{"All"}{"All_total"}{$event};
		$total += $Result{"AS"}{"All"}{"All_total"}{$event};
	}
}
foreach my $event (sort @AS_type) {
	if(exists $Result{"AS"}{"All"}{"All_total"}{$event}) {
		print OUT "\t",sprintf("%.2f",$Result{"AS"}{"All"}{"All_total"}{$event}/$total*100),"%";
	} else {
		print OUT "\t0%";
	}
}
#foreach my $events (@Event) {
	#print OUT "\t",sprintf("%.2f",$events/$total*100),"%";
#}
print OUT "\t",sprintf("%.2f",$total/$total*100),"%\n";
close OUT;


$outfile = "Classification_of_all_RASE_events.xls";
open(OUT,">$outfile") || die "$!\n";
print OUT "Type\t";
foreach my $grou (sort keys %{$Result{"RAS"}}) {
	print OUT $grou,"\t\t";
}
print OUT "\n";
#print OUT join("\t\t\t",sort keys ),"\n";
print OUT "\t",join("\t",("Up\tDown") x scalar(keys %{$Result{"RAS"}})),"\n";
my %allrase = ();
foreach my $group (sort keys %{$Result{"RAS"}}) {
	$allrase{$group}{"up"} = 0;
	$allrase{$group}{"down"} = 0;
}
foreach my $type (sort @AS_type) {
	print OUT $type;
	foreach my $group (sort keys %{$Result{"RAS"}}) {
		if (exists $Result{"RAS"}{$group}{"all"}{$type}) {
			$Result{"RAS"}{$group}{"all"}{$type}{"up"} = 0 if not defined($Result{"RAS"}{$group}{"all"}{$type}{"up"});
			$Result{"RAS"}{$group}{"all"}{$type}{"down"} = 0 if not defined($Result{"RAS"}{$group}{"all"}{$type}{"down"});
			# $Result{"RAS"}{$group}{"all"}{$type}{"total"} = 0 if not defined($Result{"RAS"}{$group}{"all"}{$type}{"total"});
			$allrase{$group}{"up"} += $Result{"RAS"}{$group}{"all"}{$type}{"up"};
			$allrase{$group}{"down"} += $Result{"RAS"}{$group}{"all"}{$type}{"down"};
			print OUT "\t",$Result{"RAS"}{$group}{"all"}{$type}{"up"},"\t",$Result{"RAS"}{$group}{"all"}{$type}{"down"};
		} else {
			print OUT "\t0\t0";
		}
	}
	print OUT "\n";
}
print OUT "Total";
foreach my $group (sort keys %{$Result{"RAS"}}) {
	print OUT "\t",$allrase{$group}{"up"},"\t",$allrase{$group}{"down"};
}
print OUT "\n";
close OUT;

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

# subroutines

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Load_config_file {
	my ($nir,$ir,$conf_ref) = @_;
	open(IN,"<$nir") || die "$!\n";
	while (<IN>) {
		chomp;
		next if(/^#|^\+/);
		my @type = split("\t");
		$conf_ref->{"AS_type"}{$type[0]} = $type[1];
	}
	close(IN);
	open(IN,"<$ir") || die "$!\n";
	while (<IN>) {
		chomp;
		next if(/^#|^\+/);
		my @type = split("\t");
		$conf_ref->{"AS_type"}{$type[0]} = $type[1];
	}
	close(IN);
}

sub Exon_exp {
	my ($dir,$result_ref) = @_;
	my @file = glob "$dir/as_stat/EXP/*base";
	#print join("\n",@file),"\n";
	foreach my $f (@file) {
		my $basename = basename($f);
		$basename =~ s/\_EXP\_exon\_base//g;
		$result_ref->{"exon_exp"}{"TotalNum"} = 0;
		open(IN,"<$f") || die "$!\n";
		while (<IN>) {
			chomp;
			next if($_!~/^\d+/);
			next if(/^BaseNumber/);
			$result_ref->{"exon_exp"}{"TotalNum"}++;
			my @line=split;
			my $total = sum(@line);
			$result_ref->{"exon_exp"}{$basename}++ if ($total>0);
		}
		close(IN);
	}
}

sub SJ_stat {
	my ($dir,$result_ref) = @_;
	# print $dir,"\n";
	my @file = glob "$dir/as_nas/*/SJ/sjNum";
	# print join("\n",@file),"\n";
	foreach my $f (@file) {
		next if ($f=~/Stat/);
		my @dir = split(/\//,$f);
		my $sample = $dir[-3];
		$sample=~s/_NAS$//;
		#print $sample,"\n";
		#exit;
		open(IN,"<$f") || die "$!\n";
		while (<IN>) {
			chomp;
			my @line = split(/\s+/);
			if(/all/) {
				$result_ref->{"SJ"}{$sample}{"all"} = $line[1];
			} elsif (/novel/) {
				$result_ref->{"SJ"}{$sample}{"novel"} = $line[1];
			} else {
				$result_ref->{"SJ"}{$sample}{"known"} = $line[1];
			}
		}
		close(IN);
	}
	my $all_sj = `wc -l $dir/as_stat/SJ/all_sj.bed`;chomp($all_sj);
	my $all_known_sj = `wc -l $dir/as_stat/SJ/all_sj_known.bed`;chomp($all_known_sj);
	my $all_novel_sj = `wc -l $dir/as_stat/SJ/all_sj_novel.bed`;chomp($all_novel_sj);
	my @num = split(/\s+/,$all_sj);
	$result_ref->{"SJ"}{"Total"}{"all"} = $num[0];
	@num = split(/\s+/,$all_known_sj);
	$result_ref->{"SJ"}{"Total"}{"known"} = $num[0];
	@num = split(/\s+/,$all_novel_sj);
	$result_ref->{"SJ"}{"Total"}{"novel"} = $num[0];
}

sub AS_stat {
	my ($dir,$result_ref) = @_;
	my @file = glob "$dir/as_nas/*/AS/AS_events*result";
	my $type = "";
	foreach my $f (@file) {
		my @dir = split(/\//,$f);
		my $sample = $dir[-3];
		$sample=~s/_NAS$//;
		if($f=~/known_result/) {
			$type = "known_result";
		} elsif ($f=~/novel_result/) {
			$type = "novel_result";
		} else {
			$type = "all_result";
		}
		open(IN,"<$f") || die "$!\n";
		while (<IN>) {
			chomp;
			my @line = split(/\s+/);
			next if ($line[-1] eq "Type");
			$result_ref->{"AS"}{$sample}{$type}{$line[-1]} = $line[-2];
		}
		close(IN);
	}
	@file = glob "$dir/as_stat/AS/*result";
	foreach my $f (@file) {
		if($f=~/all\_available\_AS/) {
			$type = "All_total";
		} elsif ($f=~/all\_available\_known/) {
			$type = "All_known";
		} else {
			$type = "All_novel";
		}
		open(IN,"<$f") || die "$!\n";
		while (<IN>) {
			chomp;
			my @line = split(/\s+/);
			next if ($line[-2] eq "\+Type");
			$result_ref->{"AS"}{"All"}{$type}{$line[-2]} = $line[-1];
		}
		close(IN);
	}
}

sub RAS_stat {
	my ($dir,$result_ref) = @_;
	my @file = glob "$dir/as_ras/*_vs_*/*result";
	my $group = "";
	foreach my $f (@file) {
		my @dir = split(/\//,$f);
		$group = $dir[-2];
		open(IN,"<$f") || die "$!\n";
		while (<IN>) {
			chomp;
			next if(/Type/);
			my @line = split(/\s+/);
			if ($dir[-1] =~ /all.*down\.result/) {
				$result_ref->{"RAS"}{$group}{"all"}{$line[0]}{"down"} = $line[1];
			} elsif ($dir[-1] =~ /all.*up\.result/) {
				$result_ref->{"RAS"}{$group}{"all"}{$line[0]}{"up"} = $line[1];
			} elsif ($dir[-1] =~ /all.*diff\.result/) {
				$result_ref->{"RAS"}{$group}{"all"}{$line[0]}{"total"} = $line[1];
			} elsif ($dir[-1] =~ /known.*down\.result/) {
				$result_ref->{"RAS"}{$group}{"known"}{$line[0]}{"down"} = $line[1];
			} elsif ($dir[-1] =~ /known.*up\.result/) {
				$result_ref->{"RAS"}{$group}{"known"}{$line[0]}{"up"} = $line[1];
			} elsif ($dir[-1] =~ /known.*diff\.result/) {
				$result_ref->{"RAS"}{$group}{"known"}{$line[0]}{"total"} = $line[1];
			} elsif ($dir[-1] =~ /novel.*down\.result/) {
				$result_ref->{"RAS"}{$group}{"novel"}{$line[0]}{"down"} = $line[1];
			} elsif ($dir[-1] =~ /novel.*up\.result/) {
				$result_ref->{"RAS"}{$group}{"novel"}{$line[0]}{"up"} = $line[1];
			} elsif ($dir[-1] =~ /novel.*diff\.result/) {
				$result_ref->{"RAS"}{$group}{"novel"}{$line[0]}{"total"} = $line[1];
			} else {
				print $f,"\n";
				next;
			}
		}
		close(IN);
	}
}

