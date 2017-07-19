less -S "$1" | perl -ne '
BEGIN{%deg1=();
%deg2=();
open IN,"'"$1"'";
while(<IN>){chomp;
@line=split;
$key="$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]";
$deg1{$key}=join("\t",@line[8..19]);
}open IN2,"'"$2"'";
while(<IN2>){chomp;
@line=split;
if($_=~/^X\.Chr/){
$head=join("\t",@line[8..19]);
next;
}
$key="$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]";
$deg2{$key}=join("\t",@line[8..19]);
}
open OUT1,">'"$3"'";
}chomp;
if($_=~/^X\.Chr/){print OUT1 $_,"\t",$head,"\n";next;}
@line=split(/\t/);
$key="$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]";
if(defined($deg2{$key})){print OUT1 $_,"\t",$deg2{$key},"\n";
}'