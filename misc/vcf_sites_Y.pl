#!/usr/bin/perl
use Getopt::Long;
my %opts = ();

GetOptions(\%opts, 'i=s', 'o=s');
#GetOptions('rmdamage' => \$damage);
#GetOptions('homozygous' => \$homozygous);

if(scalar(keys(%opts)) < 2){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}
my ($line, @seq, $alt, $ref);

if(!open(VCF, $opts{'i'})){
  die "The VCF file cannot be found.\n"
}
open(OUTPUT, ">$opts{'o'}") || die "Cannot create output destination.\n";

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

    # Matching reference or alternative allele;
    # print position and ref|alt
    if ($seq[9] ~= m/^0:/){
        print OUTPUT "$seq[1]\t$ref\n";
    }elsif($seq[9] ~= m/^1:/){
        print OUTPUT "$seq[1]\t$alt\n";
    }
}