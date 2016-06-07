#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: annotateTree.pl
#
#        USAGE: ./annotateTree.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 05/27/2016 05:24:43 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 3){print "\n$0 <tree> <meta> <c,f>\n\n"; exit;}

writeScript($ARGV[0],$ARGV[1]);
`Rscript tree.r`;


sub writeScript{

open OUT, ">tree.r" or die;

print OUT '
library(ape)
';
print OUT "tree.raw<-read.tree(\"$ARGV[0]\")\n";
print OUT "meta<-read.table(\"$ARGV[1]\")\n";
print OUT '
tree.raw2<-drop.tip(tree.raw,setdiff(meta$V1,tree.raw$tip.label))
tree.raw3<-drop.tip(tree.raw2,setdiff(tree.raw$tip.label,meta$V1))
meta<-meta[match(tree.raw3$tip.label,meta$V1),]
if (length(is.na(meta$V2))>0){
	meta$V2[which(is.na(meta$V2))]<-"NA"
}
tree<-tree.raw3
x11()
';
print OUT "plot(tree,show.tip.label=F,type=\"$ARGV[2]\")";
print OUT '
library(RColorBrewer)
cols.uniq<-colorRampPalette(brewer.pal(11,"Set3"))(length(unique(meta$V2)))
meta.uniq<-unique(meta$V2)
symb.uniq<-rep(c(21,22),(round(length(meta.uniq)/2)+1))
symb<-symb.uniq[match(meta$V2,meta.uniq)]
cols.uniq[match("NA",meta.uniq)]<-"black"

cols<-cols.uniq[match(meta$V2,meta.uniq)]
tiplabels(pch=symb,tip=match(meta$V1,tree$tip.label),bg=cols)
temp<-locator(1)
legend(temp$x,temp$y,fill=cols.uniq,legend=meta.uniq)

locator(1)
';

}


