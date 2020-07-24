#!/bin/awk -f

BEGIN{row=1; cds=0}

{
    # Look for an entry (chromosome/genome)
    if ($0 ~ /LOCUS/) {
        filename = $2".tsv"
        print "start\tend\tCDS_type\tName" > filename
        next
    }
        
    if ($0 ~ /CDS/) {
        # Save start and end positions
        split(gensub(/<|>/, "", "G" , $2), array, /\.\./)
        for (i=1; i<=2; i++) {
            cds_array[row, i] = array[i]
        }
        cds = 1
        next
    }

    if (cds &&  ($0 ~ /\/gene/ || $0 ~ /\/product/)) {
        # Save type and name of CDS
        split(gensub(/\s+\//, "", 1, $0), array, /=/)
        for (i=3; i<=4; i++) {
            cds_array[row, i] = array[i-2]
        }
        row += 1
        cds = 0
        next
    }

    # Print into file the table
    if ($0 ~ /^\/\//) {
        for (i=1; i<=row; i++){
            print cds_array[i, 1]"\t"cds_array[i, 2]"\t"cds_array[i, 3]"\t"cds_array[i, 4] >> filename
        }
        row = 1
        # Empty the array
        split("", cds_array, ":")
        next
    } 
}