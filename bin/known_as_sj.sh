#!/bin/bash
perl "$1" -gff genome_known_as_SE_uniq -t SE -sj "$2" -o genome_known_as_SE_uniq_sj 
perl "$1" -gff genome_known_as_EE_uniq -t EE -sj "$2" -o genome_known_as_EE_uniq_sj 
perl "$1" -gff genome_known_as_II_uniq -t II -sj "$2" -o genome_known_as_II_uniq_sj 
perl "$1" -gff genome_known_as_IR_uniq -t IR -sj "$2" -o genome_known_as_IR_uniq_sj 
perl "$1" -gff genome_known_as_MXE_uniq -t MXE -sj "$2" -o genome_known_as_MXE_uniq_sj 
perl "$1" -gff genome_known_as_A3SS_uniq -t A3SS -sj "$2" -o genome_known_as_A3SS_uniq_sj 
perl "$1" -gff genome_known_as_A5SS_uniq -t A5SS -sj "$2" -o genome_known_as_A5SS_uniq_sj 
perl "$1" -gff genome_known_as_AFE_uniq -t AFE -sj "$2" -o genome_known_as_AFE_uniq_sj 
perl "$1" -gff genome_known_as_ALE_uniq -t ALE -sj "$2" -o genome_known_as_ALE_uniq_sj 
