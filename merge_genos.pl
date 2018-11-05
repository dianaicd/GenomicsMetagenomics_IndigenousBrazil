#!/usr/bin/perl
=head1 NAME
merge_geno.pl
=head1 DESCRIPTION
This script left merges two beagle genotype likelihoods file.
=head1 AUTHOR
Cruz-Davalos Diana Ivette
=head1 VERSION
1.0
=head1 USAGE
perl merge_genos.pl -g1 beagle1.beagle -g2 beagle2.beagle -id ids.txt
-o output.beagle -homozygous [yes|no]
=head1 INPUT FORMAT
Beagle file as specified in ANGSD website.

=head1 OUTPUT FORMAT
Beagle file as specified in ANGSD website.

=cut


#=============================================================================#
# The action starts here
#=============================================================================#
use Getopt::Long;
#my $homozygous = '';
my %opts = ();
GetOptions(\%opts, 'g1=s', 'g2=s', 'id=s', 'homozygous:s', 'o:s');
#GetOptions('homozygous' => \$homozygous);


if(scalar(keys(%opts)) < 3){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}

my ($line, @genos, $alt, $ref, $length,$x, $k, $j, %id, $ind, $flag,
$nmarkers, $linex);
$flag = 1;
$nmarkers = 0;
my $is_first_line = 1;

if(!open(GENO1, $opts{'g1'})){
  die "The first file cannot be found.\n"
}
if(!open(GENO2, $opts{'g2'})){
  die "The second file cannot be found.\n"
}

if(!open(IDs, $opts{'id'})){
  die "The IDs file cannot be found.\n"
}

open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";
#open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

while (<IDs>){
  chomp($_);
  $id{$_} = 1;
}


while(<GENO2>){
  chomp($_);

  if($opts{'homozygous'} eq "yes"){
    $_  =~ s/\t0.333333\t0.333333\t0.333333/\t0.5\t0\t0.5/g;
  }
  @genos = split(/\t/, $_);
  if($is_first_line){
    $length = scalar(@genos);
    $ind = $length - 3;
    $is_first_line = 0;
  }
  $x = scalar(@genos);
  if($x != $length){
    print "Warning: missing genotypes on position $genos[0] $genos[1]\n.";
    print "Length of vector should be $length but it is $x\n."
  }
  if (exists $id{$genos[0]}){
    $id{$genos[0]} = [@genos[3..$length-1]];
  }else{
    $linex = scalar(@genos);
  }

}
print "nInd: $ind\n";
if($linex >1){
  print "linex: $linex\n";
}else{
  #print "entered here\n";
}
close GENO2
        or warn $! ? "Error closing : $!"
                   : "Exit status $? ";

my $i = 0;
while(<GENO1>){
  chomp($_);
  @genos = split(/\t/, $_);
  $length = scalar(@genos);

  if($flag){
    print OUTPUT $_;
    $k = 3;
    while($k <= $ind){
      $j = ($k + $length)/3 - 2;
      print OUTPUT "\tInd$j\tInd$j\tInd$j";
      $k = $k + 3;
    }
    $flag = 0;
    print OUTPUT "\n";
  }else{
      $i = $i + 1;
      print  OUTPUT $_;
      #
      $length = scalar(@{$id{$genos[0]}});
      if( $length > 1){
          print  OUTPUT "\t";
          print  OUTPUT join("\t", @{$id{$genos[0]}});

      }else{
        $k = 3;
        $nmarkers = $nmarkers + 1;
        #if k starts in 3, then while goes up to $k <= $ind
        while($k <= $ind){
          if($opts{'homozygous'} eq "yes"){
            print OUTPUT "\t0.5\t0\t0.5";
          }else{
            print OUTPUT "\t0.333333\t0.333333\t0.333333";
          }
          $k = $k + 3;
        }
      }
      print OUTPUT "\n";
  }


}

close GENO1
        or warn $! ? "Error closing : $!"
                   : "Exit status $? ";
close OUTPUT
        or warn $! ? "Error closing : $!"
                  : "Exit status $? ";

print "Panel and samples merged on $i sites.\n";
#print $nmarkers;
#print " not present in g2\n";
