## Salvelinus alpinus GWAS figure based on S. fontinalis genetic map
B. Sutherland
2016-11-01

### Overview
Steps:
1) Anchor anonymous markers from Salvelinus alpinus onto the genetic map of S. fontinalis
   This will use the recent high-density genetic map of S. fontinalis (Sutherland et al. 2016) in combination with anonymous markers.
2) Combine the positioned markers with Fst values and plot in a GWAS figure.


### 1. Anchor

**Input Data:**
Salp sequence file: 
Sfon MapComp file:

**Pipeline:**
'mapcomp_iterative'

Methods:
# Replace ‘alltags’ with ‘Salp.anon’, and LG ‘0’ to ‘1’ in sequence csv file
sed 's/alltags/Salp.anon/g' list_extracted_tags_recd_all_tags_no_fst_JS_2016-06-15.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv


# Remove ‘>’ and header from Sfon csv file, and add the 0 for the totpos position
sed 's/>//g' Sfon_v4.3_female_map.csv | grep -vE '^species' | awk -F, '{ print $1","$2","$3","0","$4","$5","$6 }' > Sfon_v4.3_female_map_no_symbol.csv


# Information on Sfon and Salp input files
wc -l Sfon_v4.3_female_map_no_symbol.csv salp.anon_markers.csv 
    3826 Sfon_v4.3_female_map_no_symbol.csv
    6230 salp.anon_markers.csv
   10056 total


# Combine Sfon and Salp anonymous markers
cat salp.anon_markers.csv Sfon_v4.3_female_map_no_symbol.csv > salp.anon_sfon_markers.csv


This output file can be used as the input for the markers.fasta file


Copy the salp.anon_sfon*.csv to the mapcomp_iterative folder
# Move to the current iterative mapcomp directory
cd /Users/wayne/Documents/bernatchez/01_Sfon_projects/12_JS_Salp_loci_rel_to_sex/mapcomp_iterative/02_data/


# Copy the current combined .csv file to the mapcomp directory
cp /Users/wayne/Documents/bernatchez/01_Sfon_projects/12_JS_Salp_loci_rel_to_sex/00_resources/salp.anon_sfon_markers.csv .


# Prepare the marker.csv file to a .fasta file
cd ..
01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/salp.anon_sfon_markers.csv 


# Check the markers.fasta 
wc -l 02_data/markers.fasta
   20112 02_data/markers.fasta


# Now MapComp is ready, at:
/Users/wayne/Documents/bernatchez/01_Sfon_projects/12_JS_Salp_loci_rel_to_sex/mapcomp_iterative


# Standard run (normal ./mapcomp single run) - Single iteration, Ssal RGI, 10Mb max dist
wc -l 03_mapped/wanted_loci.info 
    2026 03_mapped/wanted_loci.info


# Make sure to add in python script needed to remove markers
cp /Users/wayne/Desktop/Scripts_Eric/fasta_remove.py ./01_scripts/


Standard Results
Standard run, 10 iterations, 10 Mbp, Ssal genome = 2070 markers total.
