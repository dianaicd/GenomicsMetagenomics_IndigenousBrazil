#!/usr/bin/perl
use strict;

my $tped   = $ARGV[0];
my $mpile = $ARGV[1];
my %tped_hash;
my $line;
my @fields;
my $key;

open(TPEDFILE, "$tped");

#### Making the hash

while(<TPEDFILE>)
{
	chomp($_);
	@fields=split(/\t/,$_);
	$key=$fields[0] ."_".$fields[3];      ## Make a hash where the key is the chr and the position coordinates 
	$tped_hash{$key}=$_;		      ## And the value is the line itself 
}
close(TPEDFILE);


#### Scan the hash with mpileup file 

my $chr;
open(MPILEUP, $mpile);

while(<MPILEUP>)
{
	chomp($_);
	@fields=split(/\t/,$_);                   ## The mpileup file is scanned the key is built with the chromosome and position coordinates
	$fields[0] =~ s/chr//g;                   ## The hash built in the first loop is queried and the SNP coordinates (the map file line) 
	$key=$fields[0] ."_".$fields[1];          ## are printed followed by the homozygous allele
	print $tped_hash{$key};
	print "\t",uc($fields[6])," ",uc($fields[6]),"\n";
}
	

exit;

