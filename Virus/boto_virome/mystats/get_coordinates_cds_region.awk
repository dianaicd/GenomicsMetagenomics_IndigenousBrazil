#!/bin/awk -f

# Script to obtain the lowest and highest positions of the merged CDS chunk in the parvovirus 
# genome.

# IMPORTANT: the start_pos value was assigned considering the length of parvovirus B19 genome.
BEGIN{start_pos=7000; end_pos=0}

{
      # Ignore the line with spaces at the end of the *.tsv file
      if ($0 ~ /^\s/) {next}

      if ($1 ~ /start/) {next}

      if ($1 < $2) {
            if ($1 < start_pos) {start_pos = $1}
            if ($2 > end_pos) {end_pos = $2}
      } else {
            if ($2 < start_pos) {start_pos = $2}
            if ($1 > end_pos) {end_pos = $1}
      }
}

END {ID=gensub(".tsv", "", "g", FILENAME); 
      printf("%s\t%s\t%s", ID, start_pos, end_pos)}
