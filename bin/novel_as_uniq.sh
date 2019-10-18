less -S novel_as_IR | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_IR_uniq 
less -S novel_as_SE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_SE_uniq 
less -S novel_as_II | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_II_uniq 
less -S novel_as_EE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_EE_uniq 
less -S novel_as_MXE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_MXE_uniq 
less -S novel_as_A3SS | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);next if($line[-1]!~/\"1\-,2\-\"/ && $line[-1]!~/\"2\-,1\-\"/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_A3SS_uniq 
less -S novel_as_A5SS | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);next if($line[-1]!~/\"1\^,2\^\"/ && $line[-1]!~/\"2\^,1\^\"/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_A5SS_uniq 
less -S novel_as_AFE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_AFE_uniq 
less -S novel_as_ALE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_ALE_uniq 
less -S novel_as_CE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > novel_as_CE_uniq 
less -S raw_introR_events | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[6]\t$line2[-1]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n > raw_introR_events_uniq
