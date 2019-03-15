#!/usr/bin/perl
=head1 NAME

=head1 DESCRIPTION
This script samples alleles from a ped file.
=head1 AUTHOR
Cruz-Davalos Diana Ivette
=head1 VERSION
1.0
=head1 USAGE

=head1 INPUT FORMAT


=head1 OUTPUT FORMAT


=cut


#=============================================================================#
# The action starts here
#=============================================================================#
use Getopt::Long;
#my $homozygous = '';
my %opts = ();
GetOptions(\%opts, 'ped=s', 'out=s');
#GetOptions('homozygous' => \$homozygous);


if(scalar(keys(%opts)) < 2){
   &PrintHelp();
}

sub PrintHelp {
   system "pod2text -c $0 ";
   exit();
}

if(!open(PED, $opts{'ped'})){
  die "The first file cannot be found.\n"
}
open(OUTPUT, ">$opts{'out'}") || die "Cannot create output destination.\n";

my($nSites, @line, $offset, $i, $allele, @new_line);

$offset = 6;

while(<PED>){
    chomp($_);
    @line = split(" ", $_);
    @new_line = @line[0..$offset-1];
    $nSites = (scalar(@line) - $offset)/2;
    foreach $i (0..$nSites){
        $allele = @line[$offset + 2*$i + int(rand(2))];
        push @new_line,$allele;
        push @new_line,$allele;
    }
    
    print OUTPUT "@new_line\n";
}
print $nSites;