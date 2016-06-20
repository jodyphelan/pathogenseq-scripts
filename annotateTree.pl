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
use Getopt::Long qw(GetOptions);

my $treeFile;
my $annFile;
my $type = "u";
my $out = "FALSE";
my $root = "FALSE";
my $png = "FALSE";

my $usage =  '
annotateTree.pl -t <tree> -a <annotation>

	--tree|-t	Tree file
	--ann|-a	Annotation text file containing two columns - first with samples sames and sedcond with metadata
	--type|-y	Type of tree [c,u,p,f]
	--out|-o	String containing name of sample to not plot
	--root|-r	String contianing name of sample to root tree at
	--png|-p	Output to png file

';


GetOptions(
    'tree|t=s' => \$treeFile,
    'ann|a=s' => \$annFile,
    'type|y=s' => \$type,
	'out|o=s' => \$out,
	'root|r=s' => \$root,
	'png|p=s' => \$png,
	
) or die $usage;

if (!$treeFile or !$annFile){
	print $usage; exit;
}


writeScript($treeFile,$annFile,$type,$out,$png);
`Rscript tree.r`;


sub writeScript{

my ($treeFile,$annFile,$type,$out,$png) = @_;

open OUT, ">tree.r" or die;

print OUT '
library(ape)
';
print OUT "tree.raw<-read.tree(\"$treeFile\")\n";
print OUT "meta<-read.table(\"$annFile\")\n";

print OUT '
tree.raw2<-drop.tip(tree.raw,setdiff(meta$V1,tree.raw$tip.label))
tree.raw3<-drop.tip(tree.raw2,setdiff(tree.raw$tip.label,meta$V1))
meta<-meta[match(tree.raw3$tip.label,meta$V1),]
if (length(which(is.na(meta$V2)))>0){
	meta$V2[which(is.na(meta$V2))]<-"NA"
}
';
print OUT "if (\"$out\"==\"FALSE\"){\n";
print OUT '
tree.raw4<-tree.raw3
} else {
';
print OUT "\ttree.raw4<-drop.tip(tree.raw3,\"$out\")\n}\n";
print OUT "if (\"$root\"==\"FALSE\"){\n";
print OUT '
tree<-tree.raw4
} else {
';
print OUT "\ttree<-root(tree.raw4,\"$root\")\n}\n";

print OUT '
x11()
';
print OUT "plot(tree,show.tip.label=F,type=\"$type\")";
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
dev.off()
';

if ($png ne "FALSE"){

	print OUT "png(\"$png\",width=1280,height=1024)\n";
	print OUT "plot(tree,show.tip.label=F,type=\"$type\")";
	print OUT '
tiplabels(pch=symb,tip=match(meta$V1,tree$tip.label),bg=cols)
legend(temp$x,temp$y,fill=cols.uniq,legend=meta.uniq)
';
	print OUT "dev.off()\n";
}

}

