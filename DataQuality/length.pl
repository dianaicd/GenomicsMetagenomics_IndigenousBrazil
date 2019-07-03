#!/usr/bin/perl

# samtools view $bam |length.pl -o ${bam}_length.txt -type [nuclear|MT|endo]
use Getopt::Long;
my %opts = ();

GetOptions(\%opts, 'o=s', 'type:s');
#GetOptions('rmdamage' => \$damage);
#GetOptions('homozygous' => \$homozygous);

if(scalar(keys(%opts)) < 1){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}
use warnings;
use strict;

my $l = 0;
my @line =();
my %lengths;
my $i = 0;
my $total_bases = 0;
my $average;
my @chromosomes;
#
# if($opts{'type'} eq "MT"){
#   @chromosomes = ("MT");
# }elsif($opts{'type'} eq 'endo'){
#   @chromosomes = (1..22);
#   push @chromosomes, ("X", "Y", "MT");
# }elsif($opts{'type'} eq 'nuclear'){
#   @chromosomes = (1..22);
#   push @chromosomes, ("X", "Y");
# }

# sub is_mito_nuclear_endo{
#   my $chromosomes = shift;
#   my @chromosomes = @{$chromosomes};
#   my %chromosomes = map {$_ => 1} @chromosomes;
#   my $chr = shift;
#   !exists($chromosomes{$chr})
# }

open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

while(<STDIN>){
  chomp($_);
  @line = split("\t", $_);
  # if(is_mito_nuclear_endo(\@chromosomes, $line[2])){
  #   next;
  # }

  $l = scalar(split("", $line[9]));
  $total_bases = $total_bases + $l;
  #print("$l\n");
  if(exists $lengths{$l}){
    $lengths{$l} = $lengths{$l} + 1;
  }else{
    $lengths{$l} = 1;
  }
  $i = $i + 1;
}

if($i >0){
  $average = $total_bases/$i;
}else{
  $average = "no reads of this type";
}


print "$average";

#my @n = sort keys %lengths;

#print "N\tLength\n";

#for my $i (@n){
foreach $i (sort keys %lengths){
  print OUTPUT "$lengths{$i}\t$i\n";
}

close OUTPUT
        or warn $! ? "Error closing : $!"
                   : "Exit status $? ";
