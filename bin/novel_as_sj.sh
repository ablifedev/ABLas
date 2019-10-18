#!/bin/bash
perl "$1" -gff novel_as_SE_uniq -t SE -sj "$2" -o novel_as_SE_uniq_sj 
perl "$1" -gff novel_as_EE_uniq -t EE -sj "$2" -o novel_as_EE_uniq_sj 
perl "$1" -gff novel_as_II_uniq -t II -sj "$2" -o novel_as_II_uniq_sj 
perl "$1" -gff novel_as_IR_uniq -t IR -sj "$2" -o novel_as_IR_uniq_sj 
perl "$1" -gff novel_as_MXE_uniq -t MXE -sj "$2" -o novel_as_MXE_uniq_sj 
perl "$1" -gff novel_as_A3SS_uniq -t A3SS -sj "$2" -o novel_as_A3SS_uniq_sj 
perl "$1" -gff novel_as_A5SS_uniq -t A5SS -sj "$2" -o novel_as_A5SS_uniq_sj 
perl "$1" -gff novel_as_AFE_uniq -t AFE -sj "$2" -o novel_as_AFE_uniq_sj 
perl "$1" -gff novel_as_ALE_uniq -t ALE -sj "$2" -o novel_as_ALE_uniq_sj 
perl "$1" -gff novel_as_CE_uniq -t CE -sj "$2" -o novel_as_CE_uniq_sj 
perl "$1" -gff raw_introR_events_uniq -t IR -sj "$2" -o raw_introR_events_uniq_sj
more raw_introR_events_uniq_sj | perl -ne 'BEGIN{open IN,"../../genome_known_as/genome_known_as_IR_uniq";%knownir=();while(<IN>){chomp;@line=split(/\t/);@dd=split(/;/,$line[8]);$key=$line[0].":".$dd[-1];$knownir{$key}=1;}}chomp;@line=split(/\t/);@dd=split(/;/,$line[8]);$key=$line[0].":".$dd[-2];if($knownir{$key}!=1){print $_,"\n";}' > raw_introR_events_uniq_sj_N
more raw_introR_events_uniq_sj | perl -ne 'BEGIN{open IN,"../../genome_known_as/genome_known_as_IR_uniq";%knownir=();while(<IN>){chomp;@line=split(/\t/);@dd=split(/;/,$line[8]);$key=$line[0].":".$dd[-1];$knownir{$key}=1;}}chomp;@line=split(/\t/);@dd=split(/;/,$line[8]);$key=$line[0].":".$dd[-2];if($knownir{$key}==1){print $_,"\n";}' > raw_introR_events_uniq_sj_K
