#!/usr/bin/perl -w

# Copyright (c)   AB_Life 2013
# Writer:         chengchao
# Program Date:   2013.

use strict;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#Before writing your programme£¬you must write the detailed time¡¢discriptions¡¢parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

if (scalar(@ARGV)<1) {
	print"Usage:$0 fafile chrlenfile\n";
	exit;
}

my $fa = $ARGV[0];
my $chrlenfile = $ARGV[1];

`less -S $fa | perl -ne 'chomp;if(\$_!~/^>/){\$n+=length(\$_);} if((\$_=~/^>/)&&(\$n!=0)){print "\$chr\t\$n\n";\$n=0;} if(\$_=~/^>/){\$chr = \$_;\$chr=~s/>//;} END{print "\$chr\\t\$n\\n";}' > $chrlenfile`;

print "Done get chromosome length from fasta file.\n";

