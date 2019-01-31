#!/usr/bin/perl
=head1 NAME
weight_alleles.pl
=head1 DESCRIPTION
This script takes as input a allele frequencies for a reference 
population and genotypes for an ancient sample.
=head1 AUTHOR
Cruz-Davalos Diana Ivette
=head1 VERSION
1.0
=head1 USAGE
perl weight_alleles.pl -refpop test_refpop.freq -ancient test.geno
=head1 INPUT FORMAT
in progress...
# Script to project an ancestral sample to a reference population

=head1 OUTPUT FORMAT

=cut
use warnings;
use Getopt::Long;
my %refpop = ();
my %opts = ();

GetOptions(\%opts, 'refpop=s', 'ancient=s');

if(scalar(keys(%opts)) < 2){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}
if(!open(REFPOP, $opts{'refpop'})){
  die "The reference population file cannot be found.\n"
}
if(!open(ANCIENT, $opts{'ancient'})){
  die "The ancestral population file cannot be found.\n"
}
#open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

# Scan the Reference population allele frequencies
my ($geno, $w_freq, $position, @line, $weight, @popdata, @alleles, $id);
# Read and descard first line (header)
<REFPOP>;
while(<REFPOP>){
    chomp($_);
    @line = split(/\t/, $_);
    $id = join("_", @line[0..1]);
    # switch the alleles if needed
    # report the derived allele first
    if($line[2] eq $line[4]){
        $line[2] = $line[3];
        $line[3] = $line[4];
    }
    # adjust the frequency; report according to the
    # derived allele
    else#($line[3] eq $line[4])
    {
        $line[5] = 1 - $line[5];
    }
    $refpop{$id} = [@line[2..6]];
}

sub weight_allele{
    my ($position, $geno) = @_;
    $weight = "NA";
    @alleles = split(//, $geno);
    @popdata = @{$refpop{$position}};
    if($alleles[0] eq $alleles[1]){
        # Homozygous ancestral?
        if($alleles[0] eq $popdata[2]){
            # weight
            $weight = 0;
        # Homozygous derived?
        }elsif($alleles[0] eq $popdata[0]){
            $weight = 1/$popdata[3];
        # Heterozygous
        }
    }elsif(($alleles[0] eq $popdata[0] & $alleles[1] eq $popdata[1]) | 
        ($alleles[0] eq $popdata[1] & $alleles[1] eq $popdata[0])){
        
        $weight = 1/(2*$popdata[3]);
    }
    #if($weight ne "NA"){ print $weight ; }
    return($weight);

}

$w_freq = "NA";

while(<ANCIENT>){
    chomp($_);
    @line = split(/\t/, $_);
    $position = join("_", @line[0..1]);
    $geno = $line[2];
   # print "$position\n";
    if(exists $refpop{$position}){
        #print "exists!";
        $w_freq = weight_allele($position, $geno);
        if($w_freq ne "NA"){
            print "$w_freq\t@{$refpop{$position}}[3]\n";
        }
    }

}