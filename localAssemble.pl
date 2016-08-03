#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd 'abs_path';
use Cwd;



my $base_dir = `pwd`;
chomp $base_dir;
my $sample;
my $chr;
my $start;
my $end;
my $minKmer = 19;
my $maxKmer = 79;
my $extraLength = 2000;
my $ref;

my $usage =  '
localAssemble.pl
	
	--sample|-s	Sample
	--chr|-c	Chromosome
	--begin|-b	Start position
	--end|-e	End Position
	--ref|-r	Reference
	--base|-d	Base directory
	--minKmer|n	Minimum Kmer			
	--maxKmer|x	Maximum Kmer
	--length|l	Extra length to add to region

';

GetOptions(
	'sample|s=s'	=>\$sample,
	'chr|c=s'	=>\$chr,
	'begin|b=s' => \$start,
	'end|e=s' => \$end,
	'ref|r=s' => \$ref,
	'minKmer|n=s'	=> \$minKmer,
	'maxKmer|x=s'	=> \$maxKmer,
	'length|l=s'	=> \$extraLength,
	'base|d=s'	=>\$base_dir,
) or die $usage;

if (!$sample or !$chr or !$start or !$end or !$ref){print $usage;exit;}

$ref = abs_path($ref);
$base_dir = abs_path($base_dir);
my $script_dir = abs_path($0);
$script_dir =~ s/localAssemble.pl//;
my $samtools = "$script_dir/samtools-1.3.1/samtools";

if (!-e "$base_dir/bam/$sample.bam"){print "$base_dir/bam/$sample.bam not found\n";exit;}

my $regionStart = $start - $extraLength;
my $regionEnd = $end + $extraLength;

if (-d "tempIndel"){`rm -r tempIndel`;}
`mkdir tempIndel`;
chdir("tempIndel");
my $LA = fullLocalAssembly($base_dir,$sample,$chr,$regionStart,$regionEnd);
if ($LA eq "FAIL"){	print "FAIL\n"; exit;}
my $res = bwaAlign($ref,$start,$end);
sam2graph($chr,$regionStart,$regionEnd,$start,$end,$base_dir);





sub fullLocalAssembly{

if (-e "contigs.fa"){
    `rm -r *`;
    print "Revmoing old run\n";
}

my ($base_dir,$sample,$chr,$start,$end) = @_;
`sambamba view -F "not (unmapped or mate_is_unmapped) and mapping_quality >=30" -o filt.bam $base_dir/bam/$sample.bam $chr:$start-$end -f bam`;


`velveth test $minKmer,$maxKmer,2 -shortPaired -bam filt.bam`;

	my %nNodes;
    my %assembly;
    my %expCov;
    my %covCut;
    for (my $i=$minKmer; $i<$maxKmer; $i=$i+2){
        my $res = `velvetg test_$i -cov_cutoff auto -exp_cov auto -clean yes| tail -3 | tr '\n' ' '`;
        $res =~ m/Estimated Coverage = ([\d\.]+).+Estimated Coverage cutoff = ([\d\.]+).+has (\d+) nodes .+n50 of (\d+)/;
        if ($res =~ /EMPTY/){
            next;
        }
        my ($exp_cov,$cov_cut,$nNode,$n50) = ($1,$2,$3,$4);
		$nNodes{$nNode} = $i;
        $assembly{$nNode}{$n50} = $i;
        $expCov{$i} = $exp_cov;
        $covCut{$i} = $cov_cut;
        print "$i\t$nNode\t$n50\t$exp_cov\t$cov_cut\n";
    }

	my $bestnNode = (sort {$a<=>$b} keys %nNodes)[0];
	print "Optimum number of contigs = $bestnNode\n";
    my $bestn50 = (sort {$b<=>$a} keys %{$assembly{$bestnNode}})[0];
	print "Optimum n50 $bestn50\n";
    my $best = $assembly{$bestnNode}{$bestn50};

    `mv test_$best k$best`;
    `rm -r test*`;
    `ln -s k$best/contigs.fa .`;

    my $exp_cov = $expCov{$best};
    my $cov_cut = $covCut{$best};



}


sub bwaAlign{
##### start #####
my ($ref,$bp1,$bp2) = @_;
`bwa bwasw $ref contigs.fa 2>> err | $samtools view -S - > bwaContigAln.sam 2>> err`;


my ($c1,$c2) = hash_contigs("contigs.fa");
my %hasContigs = %{$c1};
my @arrContigs = @{$c2};

my %hits;
open F, "bwaContigAln.sam" or die;
my $i = 0;
while (<F>){
    chomp;
    my ($contig) = (split /\s+/,$_)[0];
    $hits{$contig} ++;
#   print "Segment found:\t$contig\n";
}
close F;

my @hits;
for( keys %hits){
    if ($hits{$_} >1){
        push @hits,$_;
#       print "Analysing $_ further\n";
    }
}

my $result = checkContigCov($bp1,$bp2);
if ($result eq "ABSENT"){
    return("ABSENT");
} elsif ($result eq "NA"){
    return("NA");
} else {

    if (scalar @hits < 1){
        die"No split contig\n";
    } else {
        open OUT, ">hits.fasta" or die;
        for (@hits){
            print OUT ">$_\n$hasContigs{$_}\n";
        }
    }
    close OUT;
    return("POSSIBLE");
}

##### end
}



sub hash_contigs {
	if (@_ != 1) {
		print "Usage: hash_contigs contigs_file";
		exit;
	}
	my $contigs_file = shift;
	if( $contigs_file =~ /\.gz$/ ){
		open(CONTIGS, $contigs_file, "gunzip -c $contigs_file |" ) or die "Cant open contigs file: $!\n";
	} else {
		open(CONTIGS, $contigs_file) or die "Cant open contigs file: $!\n";
	}
	my %contigs_hash; # hash to store contig names and sequences
		my $contigName;
	my @ids;
	while (<CONTIGS>) {
		if (/^>(\S+)/) {
			$contigName=$1;
			push @ids, $contigName;
		} else {
			chomp;
			$contigs_hash{$contigName} .= $_;
		}
	}
	close(CONTIGS);
	return(\%contigs_hash,\@ids);
}




sub checkContigCov{
##### start ######

my ($bp1,$bp2) = @_;
my %contigs;
open F, "bwaContigAln.sam" or die;
while(<F>){
    chomp;
    my ($contig,$chr,$start,$cigar) = (split /\s+/,$_)[0,2,3,5];
    if ($cigar eq "*"){next;}
    if ($cigar =~ /S/){
        $cigar =~ s/\d+S//g;
    }
    if ($cigar =~ /I/){
        $cigar =~ s/I/M/g;
    }
    if ($cigar =~ /D/){
        $cigar =~ s/D/M/g;
    }

    my @m = split /M/, $cigar;
    my $m;
    for (@m){
        $m += $_;
    }
    my $end = $start+$m;
#   print "$contig\t$chr\t$start\t$m\n";
    for (my $i = $start; $i<=$end; $i++){
        $contigs{$contig}{$i}++;
    }
}
close(F);

foreach my $contig (keys %contigs){
    my $delLen = $bp2-$bp1;
    my $delCover=0;
    for (my $i=$bp1; $i<$bp2; $i++){
        if (exists($contigs{$contig}{$i})){
            $delCover ++;
        }
    }
#   print "\n\n$contig\t$delLen\t$delCover\n";
    if ($delLen eq $delCover){return("ABSENT");}

    my $flank1Start = $bp1 - 50;
    my $flank1End = $bp1;
    my $flank2Start = $bp2;
    my $flank2End = $bp2 + 50;

    my $flank1Cover=0;
    my $flank2Cover=0;

    for (my $i=$flank1Start; $i<$flank1End; $i++){
        if (exists($contigs{$contig}{$i})){
            $flank1Cover ++;
        }
    }
    for (my $i=$flank2Start; $i<$flank2End; $i++){
        if (exists($contigs{$contig}{$i})){
            $flank2Cover ++;
        }
    }
    my $flank1CoverPCT = $flank1Cover/50;
    my $flank2CoverPCT = $flank2Cover/50;
#   print "$contig\t$flank1CoverPCT\t$flank2CoverPCT\n";
    if ($flank1CoverPCT>0.50 && $flank2CoverPCT>0.50){return("POSSIBLE")}
}
return("NA");



########## end ############
}




sub sam2graph{

my ($regionChr,$regionStart,$regionEnd,$bp1,$bp2,$base_dir) = @_;
`echo "$bp1 $bp2" > delly.coords.txt`;
open F, "bwaContigAln.sam" or die;
open OUT, ">samGraph.txt" or die;
$regionStart -= 1000;
$regionEnd += 1000;
while(<F>){
    chomp;
    my ($contig,$chr,$start,$cigar) = (split /\s+/,$_)[0,2,3,5];
    if ($regionChr ne $chr){next;}
    if ($cigar eq "*"){next;}
    if ($cigar =~ /S/){
        $cigar =~ s/\d+S//g;
    }
    if ($cigar =~ /I/){
        $cigar =~ s/I/M/g;
    }
    if ($cigar =~ /D/){
        $cigar =~ s/D/M/g;
    }

    my @m = split /M/, $cigar;
    my $m;
    for (@m){
        $m += $_;
    }
    my $end = $start+$m;
    if ($start<$regionStart or $end>$regionEnd){next;}
    print OUT "$contig\t$start\t$end\n";
}
close(F);

plotResults($base_dir,$sample,$regionStart,$regionEnd);
#### close ####
}




sub plotResults{
#### start #####

my ($base_dir,$sample,$regionStart,$regionEnd) = @_;

`$samtools view filt.bam| cut -f4,9 > insertSizes.txt`;
`tabix $base_dir/coverage/$sample.coverage.gz $chr:$regionStart-$regionEnd | cut -f2,3 > coverage.txt`;

open R, ">sam2graph.r" or die;
print R '

contigs<-read.table("samGraph.txt",stringsAsFactors=F)
xmax<-max(c(contigs$V2,contigs$V3))
xmin<-min(c(contigs$V2,contigs$V3))

ymax<-dim(contigs)[1]
ymin<-0-(2*ymax);
#png("results.png", 800,600)
x11()
plot (xmax,ymax,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch="",yaxt="n")

uniq.contigs<-unique(contigs$V1)

for (i in 1:dim(contigs)[1]){
        height<-match(contigs$V1[i],uniq.contigs)
        rect(contigs$V2[i],height-1,contigs$V3[i],height,col="light grey")

}


yaxis.points<-seq(1:length(uniq.contigs))-0.5
yaxis.names<-substr(uniq.contigs,1,7)
axis(2,yaxis.points,yaxis.names)

segments(xmin,0,xmax,0,lwd=5)


ins<-read.table("insertSizes.txt")
ins$V2<-abs(ins$V2)
if(length(which(ins$V2>10000))>0){
    ins<-ins[-(which(ins$V2>10000)),]
}

par(new=T)
ymax<-max(ins$V2)*3.5

plot (xmax,ymax,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch="",yaxt="n")
points(ins$V1,(ins$V2+max(ins$V2)),pch=20)
axis(2,(seq(0,max(ins$V2),500)+max(ins$V2)),seq(0,max(ins$V2),500))

del<-read.table("delly.coords.txt")
rect(del$V1,ymin,del$V2,ymax,col=paste(substr("#00FF00FF",1,7),40,sep=""),border=F)
segments(xmin,0+max(ins$V2),xmax,0+max(ins$V2),lwd=5)

#res<-read.table("results.txt")
#if (res$V1[1] != "NO"){
#        for (i in 1:dim(res)[1]){
#                rect(res$V5[i],ymin,res$V6[i],ymax,col=paste(substr("#FF00FFFF",1,7),40,sep=""),border=F)
#        }
#}

cov<-read.table("coverage.txt")
par(new=T)
ymax<-max(cov$V2)*4.5

xx<-c(cov$V1,rev(cov$V1))
    yy<-c(cov$V2,rep(0,length(cov$V2)))
plot (xmax,ymax,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch="",yaxt="n")
    polygon(xx,yy,col="light blue",ylim=c(0,ymax),border="NA")

segments(xmin,0,xmax,0,lwd=5)
locator(1)

dev.off()
';
close R;
`Rscript sam2graph.r 2>> err`;

#### end ####
}

