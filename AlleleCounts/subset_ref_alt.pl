#!/usr/bin/perl
use warnings;
use Getopt::Long;
my %opts = ();

# Get options and open files
GetOptions(\%opts, 'i=s', 'o=s', 'sites:s');
print "input: $opts{'i'}";

open(COUNTS, '<:gzip', $opts{'i'}
    || die "Could not open counts file from $opts{'i'}: $!");

open(SITES, $opts{'sites'}
    || die "Could not open sites file $opts{'sites'}: $!");

open(OUTPUT, '>:gzip', $opts{'o'}
    || die "Could not write to output file $opts{'o'}: $!");

# Allele code by ANGSD
my %marker = (
	"A" => 0,
	"C" => 1,
	"G" => 2,
	"T" => 3
);

my %alleles;
my @sites = ();
my @line;
my $key;
my $i;
my $firstline;
my $nInd;
my @alt_ref;
my @counts;
my $j;
my $k;
my $l;

while(<SITES>){
    chomp($_);
    $_ =~ s/\t/_/ ;

    # Change numbers to bases
    #foreach $key keys %marker{
    #    $_ =~ s/$marker{$key}/$key/;
    #}
    
    @line = split("\t", $_);
    push @sites, $line[0];
    $alleles{$line[0]} = [@line[1..2]];
}

$i = 0;
$firstline = <COUNTS>;
chomp($firstline);
$nInd = scalar(split("\t", $firstline))/4;

print("$nInd\n");
while(<COUNTS>){
    chomp($_);
    @line = split("\t", $_);
    @counts = ();

    # Get alleles
    @alt_ref = @{$alleles{$sites[$i]}}; 
    print "alt_ref @alt_ref\n";
    # Get the columns or reference and alternative alleles
    foreach $j (0..$nInd-1){
        #print "ind$nInd\n";
        $k = $j*4 + $alt_ref[0];
        $l = $j*4 + $alt_ref[1];

        push @counts, $line[$k];
        push @counts, $line[$l];
        print "index $k $l\n";
    }

    print {OUTPUT} join("\t",@counts);
    print {OUTPUT} "\n";
    $i++;
}