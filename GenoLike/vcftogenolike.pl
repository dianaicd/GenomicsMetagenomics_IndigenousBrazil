#!/usr/bin/perl
=head1 NAME
vcftogenolike.pl
=head1 DESCRIPTION
This script takes as input a VCF file and outputs
 a beagle genotype likelihoods file.
=head1 AUTHOR
Cruz-Davalos Diana Ivette
=head1 VERSION
1.0
=head1 USAGE
perl vcftogenolike.pl -i INPUT.vcf -o OUTPUT -rmdamage [yes\no]
 -homozygous [yes\no]
=head1 INPUT FORMAT
VCF file with GT fields. Keep only di-allelic sites.


=head1 OUTPUT FORMAT

=cut

#my $damage = '';
#my $homozygous = '';

use Getopt::Long;
my %opts = ();

GetOptions(\%opts, 'i=s', 'o=s', 'rmdamage:s', 'homozygous:s');
#GetOptions('rmdamage' => \$damage);
#GetOptions('homozygous' => \$homozygous);

if(scalar(keys(%opts)) < 2){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}

#=============================================================================#
# The action starts here
#=============================================================================#

my ($line, @seq, $alt, $ref, $nind,$x, $k, $j, $damaged, $di, $is_first_line);

if(!open(VCF, $opts{'i'})){
  die "The VCF file cannot be found.\n"
}
open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

my %bases = (
	"A" => "1\t0\t0",
	"C" => "0\t1\t0",
	"G" => "0\t0\t1",
	"T" => "0\t0\t0"
);

my %marker = (
	"A" => 0,
	"C" => 1,
	"G" => 2,
	"T" => 3
);

my $nsites = 0;

print OUTPUT "marker\tallele1\tallele2";
my $flag = 1;
$damaged = 0;
$di = 0;
$is_first_line = 1;

while (<VCF>){
	chomp($_);
	$line = $_;
	@seq = split(/\t/,$line);

  # Skip header lines
	if($seq[0] =~ m/#/){
		next;
	}

	$alt = $seq[4];
	$ref = $seq[3];

  # Print sites that are not di-allelic (but not to OUTPUT)
  if($alt eq "\."){
    #print "Position $seq[0] $seq[1] $seq[2] $ref $alt\n";
    $di = $di + 1;
    next;
  }
  # Filter out positions that could be confounded with damage
	if($opts{'rmdamage'} eq "yes"){
		if(($ref eq "C" && $alt eq "T") || ($ref eq "T" && $alt eq "C")){
      $damaged = $damaged + 1;
			next;
		}elsif(($ref eq "G" && $alt eq "A") || ($ref eq "A" && $alt eq "G")){
      $damaged = $damaged + 1;
			next;
		}
	}

  # number of individuals should be number of columns minus 9 first columns
  # minus one as we index starting on 0
  if($is_first_line){
    	$nind = scalar(@seq) - 9;
      print "Total individuals: $nind. \n";
      $is_first_line = 0;
      $nind = $nind - 1;

		  foreach my $k (0..$nind){
			     print OUTPUT "\tInd$k\tInd$k\tInd$k";
		  }print OUTPUT "\n";
	}

	print OUTPUT "$seq[0]_$seq[1]\t$marker{$ref}\t$marker{$alt}";
  $nsites = $nsites + 1;
	foreach my $k (9..$nind+9) {
    print OUTPUT "\t";
    if ($seq[$k] =~ m/0\/1/ || $seq[$k] =~ m/1\/0/){
              if($opts{'homozygous'} eq "yes"){
                $x = int(rand(1));
                if($x){
                  print OUTPUT "1\t0\t0";
                }else{
                  print OUTPUT "0\t0\t1";
                }
              }else{
                print OUTPUT "0\t1\t0";
              }
  	}elsif($seq[$k] =~ m/1\/1/){
			print OUTPUT "0\t0\t1";
		}elsif($seq[$k] =~m/0\/0/){
			print OUTPUT "1\t0\t0";
		}else{
      if($opts{'homozygous'} eq "yes"){
        print OUTPUT "0.5\t0\t0.5";
      }else{
        print OUTPUT "0.333333\t0.333333\t0.333333";
      }

		}

	}
	print OUTPUT "\n";
}

close OUTPUT
        or warn $! ? "Error closing : $!"
                   : "Exit status $? ";

print "Damaged sites: $damaged \nNot di-allelic site: $di\n";
print "Panel with $nsites sites";
