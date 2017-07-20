## Salpinus GWAS figure based on _S. alpinus_ genetic map
B. Sutherland
2017-07-20

### Overview
A) Anchor anonymous markers from _Salvelinus alpinus_ onto the genetic map of _S. alpinus_    
This will use the high-density genetic map of _S. alpinus_ (Nugent et al. 2017) in combination with anonymous markers.    
B) Combine the positioned markers with Fst values and plot in a GWAS figure.   

Clone this repo, run all code from within the main repo   

### Input Data
Put the following data into `02_data`    
From: [Supplemental Data from Nugent et al. 2017, G3](http://www.g3journal.org/content/7/2/543.supplemental)
Download the supplemental files and take the following files   
* Salp marker file: `FileS1.xlsx`    
* Salp map file: `FileS2.xlsx`    

From FileS1.xlsx, save out the sheet labeled 'Map_SNPs' as a .csv file and title it as: `FileS1.csv`    
Collect marker name and tag sequence from this file:
`grep -v 'Polymorphism' FileS1.csv | awk -F, '{ print $1 "," $4 }' > salp_marker_and_seq.csv`   

From FileS2.xlsx, save out the only sheet as a .csv file and title it as `FileS2.csv`    
Collect only the lines with female linkage groups that have markers with positions:    
`awk -F, '{ print $1","$3","$4 }' FileS2.csv | sed 's/,AC-/,AC/g' | sed 's/,-/,empty/g'  | grep -vE 'NA|empty' - | grep -v 'Marker,Female,Map' - > ./salp_female_map.csv`

Finally, to finish preparing the input data, go to R to make some final adjustments to prepare for `MapComp`    
In addition to format adjusting, this will also change linkage groups AC-20 and AC-4 from the current format of split by AC-20a and b to one continuous linkage group with a cumulative cM position.   

##### STILL TO CORRECT #####

#  From: [Figshare data](https://doi.org/10.6084/m9.figshare.5051821.v2)    
#  * Salp sequence file: `salp_tags.csv`    
#  * Sfon genetic map information: `LG_plot.RData`
#  * Salp outliers file: `94snps-outliers_2016-11-11.txt`   
#  * Fst values (all markers): `Fst_6147SNPs_2016-11-11.txt`   

##### END STILL TO CORRECT #####


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
[MapComp](https://github.com/enormandeau/mapcomp)   

Clone MapComp   
`git clone https://github.com/enormandeau/mapcomp.git`

Follow instructions given at the top of the MapComp iterative script:  
`01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`  

To obtain results as in the manuscript, run MapComp iterative with 10 iterations, with the default distance setting, and use the Atlantic Salmon reference genome as the genome intermediate:   
[ICSASG_v2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1)  


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
`01_scripts/GWAS_from_MapComp_2016-11-02.R`    
# This script will merge the Fst and positional data, and plot using the information on the Brook Charr map    
# After running this script, you will have a GWAS figure with Fst by Brook Charr linkage group   
```
