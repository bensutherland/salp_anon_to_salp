## Salpinus GWAS figure based on S. fontinalis genetic map
B. Sutherland
2017-05-30

### Overview
A) Anchor anonymous markers from Salvelinus alpinus onto the genetic map of S. fontinalis    
This will use the recent high-density genetic map of S. fontinalis (Sutherland et al. 2016) in combination with anonymous markers.    
B) Combine the positioned markers with Fst values and plot in a GWAS figure.   

Clone this repo, run all code from within the main repo   

### Input Data
Put the following data into `02_data`    
From: `https://academic.oup.com/gbe/article-lookup/doi/10.1093/gbe/evw262`   
Download the supplemental files and take the following file   
* Sfon Map file: `additional_fileS3_sfon_female_map.txt`   

From: `https://doi.org/10.6084/m9.figshare.5051821.v1`    
* Salp sequence file: `salp_tags.csv`    
* Sfon genetic map information: `LG_plot.RData`
* Salp outliers file: `94snps-outliers_2016-11-11.txt`   
* Fst values (all markers): `Fst_6147SNPs_2016-11-11.txt`   


### A. Prepare Data  

```
# Move to the data folder
cd 02_data

# Replace ‘alltags’ with ‘Salp.anon’, and LG ‘0’ to ‘1’ in sequence csv file
sed 's/alltags/Salp.anon/g' salp_tags.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv

# Remove header from Sfon map file, and add the 0 for the totpos position
grep -vE '^species' additional_fileS3_sfon_female_map.txt | awk '{ print $1","$2","$3","0","$4","$5 }' > sfon_markers.csv

# Confirm information on Sfon and Salp input files
wc -l *markers.csv
`6230 salp.anon_markers.csv`
`3826 sfon_markers.csv`

# Combine Sfon and Salp anonymous markers to make input for `MapComp`
cat salp.anon_markers.csv sfon_markers.csv > salp.anon_sfon_markers.csv

# Move out of the repo
cd ../../

```
Obtain MapComp iterative through the MapComp repo:  
`mapcomp` https://github.com/enormandeau/mapcomp   

Clone MapComp   
`git clone https://github.com/enormandeau/mapcomp.git`

Follow instructions given at the top of the MapComp iterative script:  
`01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`  

To obtain results as in the manuscript, run MapComp iterative with 10 iterations, with the default distance setting, and use the Atlantic Salmon reference genome as the genome intermediate:   
ICSASG_v2 https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1  


### B. Prepare and Run MapComp Iteratively
```
# Copy `salp.anon_markers.csv` to the `mapcomp/02_data` folder   

# Move to the main mapcomp directory

# Prepare the marker.csv file to a fasta file
./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/salp.anon_sfon_markers.csv

# Check the markers.fasta 
wc -l 02_data/markers.fasta
`20112 02_data/markers.fasta`

# Prepare MapComp
# Set the species name in the iterative mapping script
# e.g.  ANON=”Salp.anon”
vi ./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh

# Set the max distance in the mapcomp script (1000000)
vi ./01_scripts/mapcomp

# Set the path to the genome file in both the following:   
vi ./mapcomp   
vi 01_scripts/01_bwa_align_reads.sh   

# Run MapComp iteratively 
./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh

# Collect results, this will be used by the GWAS script
awk '{ print $1","$5","$11 }' 03_mapped/pairings_out.txt > 03_mapped/Salp_mname_Sfontotpos.csv

# Copy the result file 03_mapped/Salp_mname_Sfontotpos.csv into the folder salp_anon_to_sfon/02_data
```

### C. Combine with Fst and plot
```
# Open the GWAS script in R and follow instructions there   
`GWAS_from_MapComp_2016-11-02.R`    
# This script will merge the Fst and positional data, and plot using the information on the Brook Charr map    
# After running this script, you will have a GWAS figure with Fst by Brook Charr linkage group   
```
