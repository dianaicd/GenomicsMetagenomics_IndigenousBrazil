use warnings;
# Hash storing number of mismatches
my %mismatch =(
	"A_C" => 0,
	"A_G" => 0,
	"A_T" => 0,
	"C_T" => 0,
	"C_G" => 0,
	"C_A" => 0,
	"T_C" => 0,
	"T_A" => 0,
	"T_G" => 0,
	"G_A" => 0,
	"G_C" => 0,
	"G_T" => 0,
	"A_N" => 0,
	"C_N" => 0,
	"G_N" => 0,
	"T_N" => 0
);

use Getopt::Long;
my %opts = ();
my ($i, @seq, @reference, $total_bases, $key, $total_poly, $j);
GetOptions(\%opts, 'ref=s', 'align=s');

if(!open(REFERENCE, $opts{'ref'})) {
	die "Cannot open reference fasta file $opts{'ref'}.\n";
}if(!open(MULALIGN, $opts{'align'})){
	die "Cannot open the multiple alignment file $opts{'align'}.\n";
}

while (<REFERENCE>){
	chomp($_);
	@reference = split("\t", $_);
	@reference = split("", $reference[1]);
	$total_bases = scalar(@reference);
}

# To debug
print "Total bases = $total_bases \n";

$total_poly = 0;
while (<MULALIGN>){
	chomp($_);
	@seq = split("\t", $_);
	@seq = split("", $seq[1]);

	for $i (0 .. $total_bases-1) {
		if($seq[$i] eq $reference[$i]){next;}
		$seq[$i] =~ s/-/N/ ;

		$key = $reference[$i];
		$key =~ s/-/N/ ;
		$key = uc $key."_".$seq[$i];
		if(exists $mismatch{$key}){
			$total_poly = $total_poly + 1;
			$mismatch{$key} = $mismatch{$key} + 1;
			if($key eq "C_T" | $key eq "G_A"){
				$j = $i + 1;
				print "MT\t$i\t$j\t";
				$j = $i."_".$j."_".$key;
				print "$j\n";
			}
		}
	}
}


foreach $i (sort keys %mismatch){
	$total_bases = $mismatch{$i}/$total_poly;
	#print "$i\t$mismatch{$i}\t$total_bases\n";
}
