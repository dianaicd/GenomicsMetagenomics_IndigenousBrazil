#!/usr/bin/perl


#=============================================================================#
# The action starts here
#=============================================================================#

# Remove not di-alllelic sites
sub remove_diallelic{

  my ($ref, $alt, $chr, $pos, $id) = @_;

  if($alt eq "\."){
    print "Position $chr $pos $id $ref $alt\n";
    return(1);
  }
}

sub remove_damage{
  my ($option, $ref, $alt) = @_;
  # Filter out positions that could be confounded with damage
  if($option eq "yes"){
    if(($ref eq "C" && $alt eq "T") || ($ref eq "T" && $alt eq "C")){
      return 1;
    }elsif(($ref eq "G" && $alt eq "A") || ($ref eq "A" && $alt eq "G")){
      return 1;
    }
  }
}

sub print_header{
  my ($length_header, $flag) = @_;

  # number of individuals should be number of columns minus 9 first columns
  # minus one as we index starting on 0
  my $nind = $length_header - 9 - 1;
  if($flag){
    foreach my $k (0..$nind){
      print OUTPUT "\tInd$k\tInd$k\tInd$k";
    }
    $flag = 0;
    print OUTPUT "\n";
  }
  return $flag;
}

sub gt_to_genolike{
  my $seq = shift;
  my $homo = shift;
  my @seq = @{$seq};
  my $length_header = scalar(@seq) - 1;
  my $nind = $length_header - 9;

  foreach my $k (9..$length_header) {
    print OUTPUT "\t";
    if ($seq[$k] =~ m/0\/1/ || $seq[$k] =~ m/1\/0/){
              if($homo eq "yes"){
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
      if($homo eq "yes"){
        print OUTPUT "0.5\t0\t0.5";
      }else{
        print OUTPUT "0.333333\t0.333333\t0.333333";
      }

    }

  }
}


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

my ($line, @seq, $alt, $ref, $nind,$x, $k, $j, $l);

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

print OUTPUT "marker\tallele1\tallele2";
my $flag = 1;

while (<VCF>){
	chomp($_);
	$line = $_;
	@seq = split(/\t/,$line);

  # Skip header lines
	if($seq[0] =~ m/#/){
		next;
	}
  $l = scalar(@seq);
  $flag = print_header($l, $flag);

	$alt = $seq[4];
	$ref = $seq[3];
  remove_diallelic($ref, $alt, $seq[0], $seq[1], $seq[2]) && next if m/1/;
  # Filter out positions that could be confounded with damage
  remove_damage($opts{'rmdamage'}, $ref, $alt) && next if m/1/;

	print OUTPUT "$seq[0]_$seq[1]\t$marker{$ref}\t$marker{$alt}";

  gt_to_genolike(\@seq, \$opts{'homozygous'});

	print OUTPUT "\n";
}
