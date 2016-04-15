#!/usr/bin/perl 
#===============================================================================
#
#	  FILE:  snpMatrix2pca.pl
#
#	 USAGE:  ./snpMatrix2pca.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#	  BUGS:  ---
#	 NOTES:  ---
#	AUTHOR:  Jody Phelan (JP), jody.phelan@lshtm.ac.uk
#      COMPANY:  LSHTM
#      VERSION:  1.0
#      CREATED:  04/11/2016 05:17:58 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (scalar @ARGV != 3){ print "\nsnpMatrix2pca.pl <snp_mat.bin> <annotations> <threads>  \n\n"; exit;}

open OUT, ">pca.r" or die;
print OUT 'library("data.table")
library(amap)
';
print OUT "meta<-read.table(\"$ARGV[1]\", sep = \"\\t\")\n";
print OUT "temp <- fread(\"$ARGV[0]\",header=T, sep=\"\\t\")\n";
print OUT 'temp<-as.data.frame(temp)
snps<-temp[,4:dim(temp)[2]]
samples<-colnames(snps)
if(sum(is.na(match(samples,meta$V1)))>0){warning("Differing number of samples in datasets");}
snps<-snps[,match(intersect(samples,meta$V1),samples)]
meta<-meta[match(intersect(meta$V1,samples),meta$V1),]
';
print OUT "dists <- as.matrix(Dist(t(snps),method=\"manhattan\", nbproc=$ARGV[2]))\n";
print OUT 'results.pca<- cmdscale(dists,k=10,eig=T)
vars <- round(results.pca$eig/sum(results.pca$eig)*100,1)
write.table(dists,"samples.dists")
';

print OUT '

meta<-meta[match(rownames(results.pca$points),meta$V1),]
temp<-names(table(meta$V2))
meta$col<-rep("NA",dim(meta)[1])
for (l in temp){meta$col[which(meta$V2==l)]<-rainbow(length(temp))[match(l,temp)]}
pcaPlot<-function(x,y){plot(results.pca$points[,x], results.pca$points[,y], col=meta$col, pch=20, xlab=paste("PC",x," (",vars[x],"%)",sep=""),ylab=paste("PC",y," (",vars[y],"%)", sep=""))}
x11()

for (i in 1:9){
	pcaPlot(i,i+1)
	loc<-locator(1)
	legend(loc$x,loc$y,fill=rainbow(length(temp)),legend=temp)
	png(paste("PC",i,"vPC",i+1,".png",sep=""))
	pcaPlot(i,i+1)
	legend(loc$x,loc$y,fill=rainbow(length(temp)),legend=temp)
	dev.off()
locator(1)
}


';
close(OUT);

my $a = `Rscript pca.r`;
