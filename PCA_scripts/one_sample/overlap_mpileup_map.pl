#!/usr/bin/perl
use strict;

my $map   = $ARGV[0];
my $mpile = $ARGV[1];
my %map_hash;
my $line;
my @fields;
my $key;

open(MAPFILE, "$map");

#### Making the hash

while(<MAPFILE>)
{
	chomp($_);
	@fields=split(/\t/,$_);
	#print $fields[0];
	$key=$fields[0] ."_".$fields[3];
	#print "$key\n";
	$map_hash{$key}=$_;
}
close(MAPFILE);

#### Scan the hash with mpileup file 
## prints a map file with only the coordinates of the overlapping positions in the mpileup file 
##

my $chr;
open(MPILEUP, $mpile);

while(<MPILEUP>)
{
	@fields=split(/\t/,$_);
	$fields[0] =~ s/chr//g;
	#print $fields[0];
	$key=$fields[0] ."_".$fields[1];
	#print $key;
	print $map_hash{$key};
	print "\n";
}
	

exit;

