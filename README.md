## Salpinus GWAS figure based on S. fontinalis genetic map
B. Sutherland
2016-11-01

### Overview
Steps:    
A) Anchor anonymous markers from Salvelinus alpinus onto the genetic map of S. fontinalis    
This will use the recent high-density genetic map of S. fontinalis (Sutherland et al. 2016) in combination with anonymous markers.    
B) Combine the positioned markers with Fst values and plot in a GWAS figure.   


### A. Anchor

**Input Data:**
Put the following data into `02_data`    
Salp sequence file: `salp_tags.csv`    
Sfon MapComp file: `Sfon_v4.3_female_map.csv`

Make sure to add in python script from Eric Normandeau to remove markers
`fasta_remove.py`

**Pipeline:**
`mapcomp_iterative`

Data Preparation
```
# Move to the data folder
cd 02_data

# Replace ‘alltags’ with ‘Salp.anon’, and LG ‘0’ to ‘1’ in sequence csv file
sed 's/alltags/Salp.anon/g' salp_tags.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv

# Remove ‘>’ and header from Sfon csv file, and add the 0 for the totpos position
sed 's/>//g' Sfon_v4.3_female_map.csv | grep -vE '^species' | awk -F, '{ print $1","$2","$3","0","$4","$5","$6 }' > sfon_markers.csv

# Confirm information on Sfon and Salp input files
wc -l *markers.csv
`6230 salp.anon_markers.csv`
`3826 sfon_markers.csv`

# Combine Sfon and Salp anonymous markers to make input for `MapComp`
cat salp.anon_markers.csv sfon_markers.csv > salp.anon_sfon_markers.csv

# Copy `salp.anon_markers.csv` to the `mapcomp_iterative` folder and 
cp salp.anon_sfon_markers.csv ./../../mapcomp_iterative/02_data/

# Move to the main directory for mapcomp_iterative
cd ./../../mapcomp_iterative

# Prepare the marker.csv file to a .fasta file
./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/salp.anon_sfon_markers.csv

# Check the markers.fasta 
wc -l 02_data/markers.fasta
`20112 02_data/markers.fasta`

# Prepare MapComp_Iterative 
Set the species name in the iterative mapping script
# e.g.  ANON=”Salp.anon”
vi ./01_scripts/remove_paired_anon_and_pair_again.sh

# Set the max distance in the mapcomp script (1000000)
vi ./01_scripts/mapcomp

# Run MapComp iteratively 
./01_scripts/remove_paired_anon_and_pair_again.sh

# Collect results, this will be used by the GWAS script
awk '{ print $1","$5","$11 }' 03_mapped/pairings_out.txt > 03_mapped/Salp_mname_Sfontotpos.csv

# Copy the result file 03_mapped/Salp_mname_Sfontotpos.csv into the folder salp_anon_to_sfon/02_data
```

Now open the GWAS script and follow instructions there
