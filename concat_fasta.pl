my $new_seq = "";
my $start = 0;

while(<STDIN>){
	chomp($_);
	if ($_ =~ m/>/ ){
		$new_seq .= "\n";
		if($start){
			print $new_seq ;
		}
		$new_seq = $_."\t" ;
	}else{
		$new_seq .= $_ ;
		$start = 1;
	}

}

print $new_seq."\n";
