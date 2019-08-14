#!/usr/bin/perl

# Script to get trailing data
# source: https://unix.stackexchange.com/questions/295702/how-to-get-trailing-data-of-gzip-archive
use strict;
use warnings; 

use IO::Uncompress::Gunzip qw(:all);
use IO::File;

unshift(@ARGV, '-') unless -t STDIN;

my $input_file_name = shift;
my $output_file_name = shift;

if (! defined $input_file_name) {
  die <<END;
Usage:

  $0 ( GZIP_FILE | - ) [OUTPUT_FILE]

  ... | $0 [OUTPUT_FILE]

Extracts the trailing data of a gzip archive.
Outputs to stdout if no OUTPUT_FILE is given.
- as input file file causes it to read from stdin.

Examples:

  $0 archive.tgz trailing.bin

  cat archive.tgz | $0

END
}

my $in = new IO::File "<$input_file_name" or die "Couldn't open gzip file.\n";
gunzip $in => "/dev/null",
  TrailingData => my $trailing;
undef $in;

if (! defined $output_file_name) {
  print $trailing;
} else {
  open(my $fh, ">", $output_file_name) or die "Couldn't open output file.\n";
  print $fh $trailing;
  close $fh;
  print "Output file written.\n";
}