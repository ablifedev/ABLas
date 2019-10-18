cat "$1"/genome_known_as/genome_known_as_SE_uniq_sj_filter "$1"/"$2"/AS/novel_as_SE_uniq_sj_filter "$1"/"$3"/AS/novel_as_SE_uniq_sj_filter > SE &&
less -S SE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > SE_uniq &&
perl "$4" -gff SE_uniq -t SE -sj "$5" -o "$2"_SE_uniq_sj &&
perl "$4" -gff SE_uniq -t SE -sj "$6" -o "$3"_SE_uniq_sj &&
less -S "$2"_SE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$line[7]\n";' > "$2"_SE_uniq_sj_avg &&
less -S "$3"_SE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$line[7]\n";' > "$3"_SE_uniq_sj_avg &&
less -S "$3"_SE_uniq_sj_avg|awk '{print $6"\t"$7}' > SE_"$3"_tmp &&
paste "$2"_SE_uniq_sj_avg SE_"$3"_tmp > SE_uniq_sj_avg &&
R --slave <"$7" --args SE_uniq_sj_avg SE_uniq_sj_avg_fisher_test ./ &&
less -S SE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > SE_uniq_sj_avg_fisher_test_final &&
less -S SE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > SE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/"$2"/AS/novel_as_CE_uniq_sj_filter "$1"/"$3"/AS/novel_as_CE_uniq_sj_filter > CE &&
less -S CE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > CE_uniq &&
perl "$4" -gff CE_uniq -t CE -sj "$5" -o "$2"_CE_uniq_sj &&
perl "$4" -gff CE_uniq -t CE -sj "$6" -o "$3"_CE_uniq_sj &&
less -S "$2"_CE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$line[7]\n";' > "$2"_CE_uniq_sj_avg &&
less -S "$3"_CE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$line[7]\n";' > "$3"_CE_uniq_sj_avg &&
less -S "$3"_CE_uniq_sj_avg|awk '{print $6"\t"$7}' > CE_"$3"_tmp &&
paste "$2"_CE_uniq_sj_avg CE_"$3"_tmp > CE_uniq_sj_avg &&
R --slave <"$7" --args CE_uniq_sj_avg CE_uniq_sj_avg_fisher_test ./ &&
less -S CE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[5]/($line[5]+$line[6]);$s=$line[7]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > CE_uniq_sj_avg_fisher_test_final &&
less -S CE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > CE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_A3SS_uniq_sj_filter "$1"/"$2"/AS/novel_as_A3SS_uniq_sj_filter "$1"/"$3"/AS/novel_as_A3SS_uniq_sj_filter > A3SS &&
less -S A3SS | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > A3SS_uniq &&
perl "$4" -gff A3SS_uniq -t A3SS -sj "$5" -o "$2"_A3SS_uniq_sj &&
perl "$4" -gff A3SS_uniq -t A3SS -sj "$6" -o "$3"_A3SS_uniq_sj &&
less -S "$2"_A3SS_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_A3SS_uniq_sj_avg &&
less -S "$3"_A3SS_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_A3SS_uniq_sj_avg &&
less -S "$3"_A3SS_uniq_sj_avg|awk '{print $6"\t"$7}' > A3SS_"$3"_tmp &&
paste "$2"_A3SS_uniq_sj_avg A3SS_"$3"_tmp > A3SS_uniq_sj_avg &&
R --slave <"$7" --args A3SS_uniq_sj_avg A3SS_uniq_sj_avg_fisher_test ./ &&
less -S A3SS_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > A3SS_uniq_sj_avg_fisher_test_final &&
less -S A3SS_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > A3SS_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_A5SS_uniq_sj_filter "$1"/"$2"/AS/novel_as_A5SS_uniq_sj_filter "$1"/"$3"/AS/novel_as_A5SS_uniq_sj_filter > A5SS &&
less -S A5SS | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > A5SS_uniq &&
perl "$4" -gff A5SS_uniq -t A5SS -sj "$5" -o "$2"_A5SS_uniq_sj &&
perl "$4" -gff A5SS_uniq -t A5SS -sj "$6" -o "$3"_A5SS_uniq_sj &&
less -S "$2"_A5SS_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_A5SS_uniq_sj_avg &&
less -S "$3"_A5SS_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_A5SS_uniq_sj_avg &&
less -S "$3"_A5SS_uniq_sj_avg|awk '{print $6"\t"$7}' > A5SS_"$3"_tmp &&
paste "$2"_A5SS_uniq_sj_avg A5SS_"$3"_tmp > A5SS_uniq_sj_avg &&
R --slave <"$7" --args A5SS_uniq_sj_avg A5SS_uniq_sj_avg_fisher_test ./ &&
less -S A5SS_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > A5SS_uniq_sj_avg_fisher_test_final &&
less -S A5SS_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > A5SS_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_AFE_uniq_sj_filter "$1"/"$2"/AS/novel_as_AFE_uniq_sj_filter "$1"/"$3"/AS/novel_as_AFE_uniq_sj_filter > AFE &&
less -S AFE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > AFE_uniq &&
perl "$4" -gff AFE_uniq -t AFE -sj "$5" -o "$2"_AFE_uniq_sj &&
perl "$4" -gff AFE_uniq -t AFE -sj "$6" -o "$3"_AFE_uniq_sj &&
less -S "$2"_AFE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_AFE_uniq_sj_avg &&
less -S "$3"_AFE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_AFE_uniq_sj_avg &&
less -S "$3"_AFE_uniq_sj_avg|awk '{print $6"\t"$7}' > AFE_"$3"_tmp &&
paste "$2"_AFE_uniq_sj_avg AFE_"$3"_tmp > AFE_uniq_sj_avg &&
R --slave <"$7" --args AFE_uniq_sj_avg AFE_uniq_sj_avg_fisher_test ./ &&
less -S AFE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > AFE_uniq_sj_avg_fisher_test_final &&
less -S AFE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > AFE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_ALE_uniq_sj_filter "$1"/"$2"/AS/novel_as_ALE_uniq_sj_filter "$1"/"$3"/AS/novel_as_ALE_uniq_sj_filter > ALE &&
less -S ALE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > ALE_uniq &&
perl "$4" -gff ALE_uniq -t ALE -sj "$5" -o "$2"_ALE_uniq_sj &&
perl "$4" -gff ALE_uniq -t ALE -sj "$6" -o "$3"_ALE_uniq_sj &&
less -S "$2"_ALE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_ALE_uniq_sj_avg &&
less -S "$3"_ALE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_ALE_uniq_sj_avg &&
less -S "$3"_ALE_uniq_sj_avg|awk '{print $6"\t"$7}' > ALE_"$3"_tmp &&
paste "$2"_ALE_uniq_sj_avg ALE_"$3"_tmp > ALE_uniq_sj_avg &&
R --slave <"$7" --args ALE_uniq_sj_avg ALE_uniq_sj_avg_fisher_test ./ &&
less -S ALE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > ALE_uniq_sj_avg_fisher_test_final &&
less -S ALE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > ALE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_II_uniq_sj_filter "$1"/"$2"/AS/novel_as_II_uniq_sj_filter "$1"/"$3"/AS/novel_as_II_uniq_sj_filter > II &&
less -S II | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > II_uniq &&
perl "$4" -gff II_uniq -t II -sj "$5" -o "$2"_II_uniq_sj &&
perl "$4" -gff II_uniq -t II -sj "$6" -o "$3"_II_uniq_sj &&
less -S "$2"_II_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_II_uniq_sj_avg &&
less -S "$3"_II_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_II_uniq_sj_avg &&
less -S "$3"_II_uniq_sj_avg|awk '{print $6"\t"$7}' > II_"$3"_tmp &&
paste "$2"_II_uniq_sj_avg II_"$3"_tmp > II_uniq_sj_avg &&
R --slave <"$7" --args II_uniq_sj_avg II_uniq_sj_avg_fisher_test ./ &&
less -S II_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > II_uniq_sj_avg_fisher_test_final &&
less -S II_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > II_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_MXE_uniq_sj_filter "$1"/"$2"/AS/novel_as_MXE_uniq_sj_filter "$1"/"$3"/AS/novel_as_MXE_uniq_sj_filter > MXE &&
less -S MXE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > MXE_uniq &&
perl "$4" -gff MXE_uniq -t MXE -sj "$5" -o "$2"_MXE_uniq_sj &&
perl "$4" -gff MXE_uniq -t MXE -sj "$6" -o "$3"_MXE_uniq_sj &&
less -S "$2"_MXE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;@dd2=split(/,/,$line[7]);$num2=$#dd2+1;$avg2=sum(@dd2)/$num2;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$avg2\n";' > "$2"_MXE_uniq_sj_avg &&
less -S "$3"_MXE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;@dd2=split(/,/,$line[7]);$num2=$#dd2+1;$avg2=sum(@dd2)/$num2;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$avg2\n";' > "$3"_MXE_uniq_sj_avg &&
less -S "$3"_MXE_uniq_sj_avg|awk '{print $6"\t"$7}' > MXE_"$3"_tmp &&
paste "$2"_MXE_uniq_sj_avg MXE_"$3"_tmp > MXE_uniq_sj_avg &&
R --slave <"$7" --args MXE_uniq_sj_avg MXE_uniq_sj_avg_fisher_test ./ &&
less -S MXE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > MXE_uniq_sj_avg_fisher_test_final &&
less -S MXE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > MXE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_EE_uniq_sj_filter "$1"/"$2"/AS/novel_as_EE_uniq_sj_filter "$1"/"$3"/AS/novel_as_EE_uniq_sj_filter > EE &&
less -S EE | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[-1]);$key="$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;print $_,"\n";' > EE_uniq &&
perl "$4" -gff EE_uniq -t EE -sj "$5" -o "$2"_EE_uniq_sj &&
perl "$4" -gff EE_uniq -t EE -sj "$6" -o "$3"_EE_uniq_sj &&
less -S "$2"_EE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;@dd2=split(/,/,$line[7]);$num2=$#dd2+1;$avg2=sum(@dd2)/$num2;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$avg2\n";' > "$2"_EE_uniq_sj_avg &&
less -S "$3"_EE_uniq_sj | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);@dd=split(/,/,$line[5]);$num=$#dd+1;$avg=sum(@dd)/$num;@dd2=split(/,/,$line[7]);$num2=$#dd2+1;$avg2=sum(@dd2)/$num2;print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$avg\t$avg2\n";' > "$3"_EE_uniq_sj_avg &&
less -S "$3"_EE_uniq_sj_avg|awk '{print $6"\t"$7}' > EE_"$3"_tmp &&
paste "$2"_EE_uniq_sj_avg EE_"$3"_tmp > EE_uniq_sj_avg &&
R --slave <"$7" --args EE_uniq_sj_avg EE_uniq_sj_avg_fisher_test ./ &&
less -S EE_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > EE_uniq_sj_avg_fisher_test_final &&
less -S EE_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > EE_uniq_sj_avg_fisher_test_final_filter

cat "$1"/genome_known_as/genome_known_as_IR_uniq_sj_filter "$1"/"$2"/AS/novel_as_IR_uniq_sj2_filter "$1"/"$2"/AS/raw_introR_events_uniq_sj2_filter "$1"/"$2"/AS/raw_introR_events_uniq_sj_K "$1"/"$3"/AS/novel_as_IR_uniq_sj2_filter "$1"/"$3"/AS/raw_introR_events_uniq_sj2_filter "$1"/"$3"/AS/raw_introR_events_uniq_sj_K> IR &&
less -S IR | perl -ne 'BEGIN{%as=();}chomp;@line=split(/\t/);@line2=split(";",$line[8]);$key="$line[0]\t$line[2]\t$line[6]\t$line2[-2]";if(defined($as{$key}->{"gene"})){$as{$key}->{"gene"}.=",$line[1]";}else{$as{$key}->{"gene"}="$line[0]\t$line[1]";$as{$key}->{"list"}=$_;}END{foreach $thiskey (sort keys %as){$tmp=$as{$thiskey}->{"list"};$tmp=~s/^\S+\s+\S+\s+//;print $as{$thiskey}->{"gene"},"\t",$tmp,"\n";}}' | sort -k 1,1 -k 4,4n | perl -ne 'chomp;$_=~s/\sID\s\S+//;$_=~s/\s\d+$//;print $_,"\n";' > IR_uniq &&
perl "$4" -gff IR_uniq -t IR -sj "$5" -o "$2"_IR_uniq_sj &&
perl "$4" -gff IR_uniq -t IR -sj "$6" -o "$3"_IR_uniq_sj &&
perl5.14 "$8" -gff "$2"_IR_uniq_sj -fa "$11" -bam "$9" -o "$2"_IR_uniq_sj2 &&
perl5.14 "$8" -gff "$3"_IR_uniq_sj -fa "${11}" -bam "${10}" -o "$3"_IR_uniq_sj2 &&
less -S "$2"_IR_uniq_sj2 | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$2"_IR_uniq_sj_avg &&
less -S "$3"_IR_uniq_sj2 | perl -ne 'use List::Util qw(sum);chomp;@line=split(/\t/);print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[7]\n";' > "$3"_IR_uniq_sj_avg &&
less -S "$3"_IR_uniq_sj_avg|awk '{print $6"\t"$7}' > IR_"$3"_tmp &&
paste "$2"_IR_uniq_sj_avg IR_"$3"_tmp > IR_uniq_sj_avg &&
R --slave <"$7" --args IR_uniq_sj_avg IR_uniq_sj_avg_fisher_test ./ &&
less -S IR_uniq_sj_avg_fisher_test | perl -ne 'chomp;@line=split;next if($line[5]+$line[6]==0 || $line[7]+$line[8]==0);$f=$line[6]/($line[5]+$line[6]);$s=$line[8]/($line[7]+$line[8]);$ratio=$s-$f;print "$_\t$ratio\n";' > IR_uniq_sj_avg_fisher_test_final &&
less -S IR_uniq_sj_avg_fisher_test_final|awk '$10<0.05 && ($11>0.2||$11<-0.2)'|less -S > IR_uniq_sj_avg_fisher_test_final_filter