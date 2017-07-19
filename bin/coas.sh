less -S "$1" | perl -ne '
BEGIN{%deg1=();
%deg2=();
open IN,"'"$1"'";
while(<IN>){chomp;
@line=split;
$key="$line[0]:$line[1]:$line[2]:$line[3]:$line[4]";
$deg1{$key}=$_;
}open IN2,"'"$2"'";
while(<IN2>){chomp;
@line=split;
$key="$line[0]:$line[1]:$line[2]:$line[3]:$line[4]";
$deg2{$key}=$_;
}
$file1="'"$1"'";
$file2="'"$2"'";
$file1=~s/\S+\///;
$file2=~s/\S+\///;
open OUT1,">'"$3"'-$file1";
open OUT2,">'"$3"'-$file2";
}chomp;
@line=split(/\s+/);
$key="$line[0]:$line[1]:$line[2]:$line[3]:$line[4]";
if(defined($deg2{$key})){print OUT1 $_,"\n";
}END{foreach $deg_gene (keys %deg2){if(defined($deg1{$deg_gene})){print OUT2 $deg2{$deg_gene},"\n";
}}}'