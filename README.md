## Salpinus GWAS figure based on S. fontinalis genetic map
B. Sutherland
2016-11-01

### Overview
A) Anchor anonymous markers from Salvelinus alpinus onto the genetic map of S. fontinalis    
This will use the recent high-density genetic map of S. fontinalis (Sutherland et al. 2016) in combination with anonymous markers.    
B) Combine the positioned markers with Fst values and plot in a GWAS figure.   

### Input Data
Put the following data into `02_data`    
Salp sequence file: `salp_tags.csv`    
Sfon MapComp file: `Sfon_v4.3_female_map.csv`
Sfon genetic map information: `LG_plot.RData`

Run MapComp iterative using the MapComp repo:  
`mapcomp` https://github.com/enormandeau/mapcomp   
Following instructions given at the top of the MapComp iterative script:  
`01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`  
To obtain results as in the manuscript, run MapComp iterative with 10 iterations, with the default distance setting, and use the Atlantic Salmon reference genome as the genome intermediate:   
ICSASG_v2 https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1  

### A. Anchor

** Data Preparation **
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
```

** MapComp Iterative Instruction **
```
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

### B. Combine with Fst and plot
Uses the following files from BAYESCAN:
94snps-outliers_2016-11-11.txt # only the outliers
Fst_6147SNPs_2016-11-11.txt # all markers

Open the GWAS script and follow instructions there   
`GWAS_from_MapComp_2016-11-02.R`    
This script will merge the Fst and positional data, and plot using the information on the Brook Charr map    
After running this script, you will have a GWAS figure with Fst by Brook Charr linkage group   
