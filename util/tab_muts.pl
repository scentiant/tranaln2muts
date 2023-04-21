#!/usr/bin/perl -w
use strict;

# Michael E. Sparks (michael.sparks2@usda.gov), 20 April 2023

# It is not technically correct to dub indels (or sequencing ambiguities) as
# non-synonymous mutations, though I've done so to simplify the programming.
# The distinctions should be obvious to any biologist, and folks are
# welcome to extend the code to handle them in a more specific manner.

my $DESC = 15; # buffer size for sequence ids
my $CHUNK = 80; # length of sequence reporting per stanza
my $len=-1; # used to record overall length of formatted pep/cds strings

my $PEP = shift;
my $CDS = shift or die "$0 pep_aln cds_aln\n";

# fasta-formatted peps aligned by muscle
# seqs should be merged onto one line
my %pep=();
open(PEP,"<$PEP") or die "$!\n";
while(<PEP>) {
  chomp $_;
  $_=~/^>(\S+)$/;
  my $k=$1;

  $_=<PEP>;
  chomp $_;
  $pep{$k}=$_;
}
close PEP;

# fasta-formatted cdss aligned by tranalign
# seqs should be merged onto one line
my %cds=();
open(CDS,"<$CDS") or die "$!\n";
while(<CDS>) {
  chomp $_;
  $_=~/^>(\S+)$/;
  my $k=$1;

  $_=<CDS>;
  chomp $_;
  $cds{$k}=$_;
}
close CDS;

my %lines=(); # formatted pep/cds lines to break in chunks and report
foreach my $key (keys %pep) {
  # peptide line
  $lines{"${key}p"}="";
  my @seq=split(//,$pep{$key});
  for(my $i=0;$i<=$#seq;++$i) {
    $lines{"${key}p"} .= " ";
    $lines{"${key}p"} .= $seq[$i];
    $lines{"${key}p"} .= "  ";
  }

  if($len > 0) {
    if(length($lines{"${key}p"}) != $len) {
      die "$key len\n";
    }
  }
  else { # set it and forget it - all lines should be of this length
    $len = length($lines{"${key}p"});
  }

  # corresponding cds line
  $lines{"${key}c"}="";
  @seq=split(//,$cds{$key});
  for(my $i=0;$i<=$#seq;$i+=3) {
    $lines{"${key}c"} .= $seq[$i];
    $lines{"${key}c"} .= $seq[$i+1];
    $lines{"${key}c"} .= $seq[$i+2];
    $lines{"${key}c"} .= " ";
  }
  if(length($lines{"${key}c"}) != $len) {
    die "$key len\n";
  }
}

# now identify synonymous and non-synonymous mutations
my @pepdiffs=();
for(my $i=0;$i<$len/4;++$i) {
  @pepdiffs[$i]=0;
  my @ps=();
  foreach my $key (keys %pep) {
    push(@ps,substr($pep{$key},$i,1));
  }
  my $init=$ps[0];
  for(my $j=1;$j<=$#ps;++$j) {
    if($ps[$j] ne $init) { $pepdiffs[$i]=1; }
  }
}

my @cdsdiffs=();
for(my $i=0,my $j=0;$i<$len/4;++$i,$j+=3) { # i indexes CDSs, j indexes nt pos
  @cdsdiffs[$i]=0;
  my @cs=();
  foreach my $key (keys %pep) { # keys of %cds are the same
    push(@cs,substr($cds{$key},$j,3));
  }
  my $init=$cs[0];
  for(my $k=1;$k<=$#cs;++$k) {
    if($cs[$k] ne $init) { $cdsdiffs[$i]=1; }
  }
}

if($#pepdiffs != $#cdsdiffs) {
  die "diff arrs of unequal size? pep: ",$#pepdiffs+1," cds: ",$#cdsdiffs+1,"\n";
}
my $diffline="";
my @ka=();
my @ks=();
for(my $i=0;$i<=$#pepdiffs;++$i) {
  $diffline .= " ";
  if($cdsdiffs[$i]) {
    if($pepdiffs[$i]) {
      $diffline .= "^";
      push(@ka,$i);
    }
    else {
      $diffline .= "=";
      push(@ks,$i);
    }
  }
  elsif($pepdiffs[$i]) {
    die "peps differ but not cds?\n";
  }
  else {
    $diffline .= " ";
  }
  $diffline .= "  ";
}
if(length($diffline) != $len) {
  die "diffline unexpected len (obs: ",length($diffline)," exp: $len)\n";
}


# report stuff
print "position";

# The @order array can be used if you do not care for the default ordering
# specified by sorting keys from %pep (e.g., if you have an ordering in mind
# that takes phylogenetic relatedness into account)
##my @order=("MT753155.1", "MN938851.1", "LdIV1_JGS",
##           "LdIV1_ZY", "LdIV1_NJ", "LdIV1_CT", "KJ629170.1"); # MY PREFERRED ORDERING FOR THIS APPLICATION
my @order=sort keys %pep; # THIS IS A GENERAL WAY OF ORDERING TAXA

for(my $z=0;$z<=$#order;++$z) {
  print "\t",$order[$z];
}
print "\n";

print "Non-synonymous mutations (",$#ka+1," positions):\n";
for(my $i=0;$i<=$#ka;++$i) {
  print $ka[$i]+1;
  for(my $z=0;$z<=$#order;++$z) {
    print "\t",substr($pep{$order[$z]},$ka[$i],1);
  }
  print "\n";
  for(my $z=0;$z<=$#order;++$z) {
    print "\t",substr($cds{$order[$z]},$ka[$i]*3,3);
  }
  print "\n\n";
}

print "\nSynonymous mutations (",$#ks+1," positions):\n";
for(my $i=0;$i<=$#ks;++$i) {
  print $ks[$i]+1;
  for(my $z=0;$z<=$#order;++$z) {
    print "\t",substr($pep{$order[$z]},$ks[$i],1);
  }
  print "\n";
  for(my $z=0;$z<=$#order;++$z) {
    print "\t",substr($cds{$order[$z]},$ks[$i]*3,3);
  }
  print "\n\n";
}

exit 0;
