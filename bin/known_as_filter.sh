more genome_known_as_SE_uniq_sj | perl -ne 'use List::Util qw(max);chomp;@line=split;@tmp=split(",",$line[5]);if(max(@tmp)>0 && $line[7]>0){print $_,"\n";}' > genome_known_as_SE_uniq_sj_filter 
more genome_known_as_MXE_uniq_sj | perl -ne 'use List::Util qw(max);chomp;@line=split;@tmp=split(",",$line[5]);@tmp2=split(",",$line[7]);if(max(@tmp)>0 && max(@tmp2)>0){print $_,"\n";}' > genome_known_as_MXE_uniq_sj_filter 
more genome_known_as_II_uniq_sj | perl -ne 'use List::Util qw(max);chomp;@line=split;@tmp=split(",",$line[5]);@tmp2=split(",",$line[7]);if(max(@tmp)>0 && max(@tmp2)>0){print $_,"\n";}' > genome_known_as_II_uniq_sj_filter 
more genome_known_as_EE_uniq_sj | perl -ne 'use List::Util qw(max);chomp;@line=split;@tmp=split(",",$line[5]);@tmp2=split(",",$line[7]);if(max(@tmp)>0 && max(@tmp2)>0){print $_,"\n";}' > genome_known_as_EE_uniq_sj_filter 
more genome_known_as_IR_uniq_sj | perl -ne 'use List::Util qw(max);chomp;@line=split;if($line[7]>0){print $_,"\n";}' > genome_known_as_IR_uniq_sj_filter 
more genome_known_as_ALE_uniq_sj | perl -ne 'chomp;@line=split;if($line[5]>0 && $line[7]>0){print $_,"\n";}' > genome_known_as_ALE_uniq_sj_filter 
more genome_known_as_AFE_uniq_sj | perl -ne 'chomp;@line=split;if($line[5]>0 && $line[7]>0){print $_,"\n";}' > genome_known_as_AFE_uniq_sj_filter 
more genome_known_as_A3SS_uniq_sj | perl -ne 'chomp;@line=split;if($line[5]>0 && $line[7]>0){print $_,"\n";}' > genome_known_as_A3SS_uniq_sj_filter 
more genome_known_as_A5SS_uniq_sj | perl -ne 'chomp;@line=split;if($line[5]>0 && $line[7]>0){print $_,"\n";}' > genome_known_as_A5SS_uniq_sj_filter