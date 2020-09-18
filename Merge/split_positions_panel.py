# Snakemake to extract positions from a panel
# and write them into files with at most n lines

# max_n_lines = config['merge_sample']['max_n_lines']
# panels = list( config['merge_sample']["panels"].keys() ) 

# just a shortcut 
# panel_dicts = { panel:dictionary for dictionary in config[ 'merge_sample' ][ 'panels' ][ panel ] }


def extract_positions_panel( panel, out_dir = 'positions/' ):
    # A panel can be in a plink format (ped, bed) or vcf/vcf.gz/bcf
    # which would determine the way in which we will extract the positions
    if not os.path.exists( out_dir ):
        os.system('mkdir ' + out_dir)

    panel_extension = panel_dicts[ panel ][ 'type' ]
    panel_path = panel_dicts[ panel ][ 'path' ]

    if panel_extension == 'vcf':
        command = "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' " + panel_path
    elif panel_extension == 'ped':
        command = 'sed "s/ \+/\\t/g" ' + panel_path.replace('.ped', '.map') + '| cut -f 1,4,5,6' 
    elif panel_extension == 'bed':
        command = 'sed "s/ \+/\\t/g" ' + panel_path.replace('.bed', '.bim') + '| cut -f 1,4,5,6' 

    positions_file = out_dir + panel + ".positions"
    command = command + " > " +  positions_file
    print( command )
    os.system( command )

    # Recode reference and alternative alleles A/C/G/T to 0/1/2/3 
    # (useful to get genotype likelihoods in ANGSD)
    ref_alt_file = out_dir + panel + '.refalt'
    command = 'sed "s/\\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g" ' + positions_file + ' > ' + ref_alt_file
    os.system( command )

def split_ref_alt_by_chromosome(panel, chr, out_dir = 'positions/' ):

    ref_alt_file = out_dir + panel + '.refalt'
    sites_file = out_dir + panel + ".sites"

    ref_alt_chr_file = out_dir + panel + '_' + chr + '.refalt'
    sites_chr_file = out_dir + panel + '_' + chr + ".sites"

    command = ' cut -f 1,2 ' +  ref_alt_file +'| grep -P "^' + chr  + '_" '+ '| sed "s/_/\t/"|'+ ' awk \'BEGIN {OFS="\t"} {print $1,$2-1,$2 }\' > ' + sites_chr_file
    os.system( command )
    command = 'grep -P "^' + chr + '\_" ' + ref_alt_file + ' > ' + ref_alt_chr_file
    os.system( command )


def split_chromosome_lines( panel, chr, max_n_lines, out_dir = 'positions/' ):
    ref_alt_chr_file = out_dir + panel + '_' + chr + '.refalt'
    sites_chr_file = out_dir + panel + '_' + chr + ".sites"

    def split_file( file, output_list ): 
        command = 'split -l {n} {file} {file}_'.format( n = max_n_lines, file = file )
        # command = 'split -l ' + str( max_n_lines ) + ' --numeric-suffixes --additional-suffix=".txt" ' + ref_alt_chr_file  + ' ' + ref_alt_chr_file + '_'
        os.system( command )
        splitted_files = sorted( glob.glob( file + "_*" ) )

        with open( output_list, 'w' ) as output:
            [ output.write( line + '\n' ) for line in splitted_files ]
    
    output_list = 'positions/' + panel + '_' + chr + "_ref_alt_list.txt"
    split_file( file = ref_alt_chr_file, output_list = output_list )
    output_list = 'positions/' + panel + '_' + chr + "_sites_list.txt"
    split_file( file = sites_chr_file, output_list = output_list )

# extract_positions_panel(panel = panel, out_dir = 'positions/' )
# [ split_ref_alt_by_chromosome( panel = panel, chr = str( chr ) ) for chr in range(1, 23) ] 
# [ split_chromosome_lines( panel = panel, chr = str( chr ), max_n_lines = max_n_lines ) for chr in range(1, 23) ] 