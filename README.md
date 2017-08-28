## Salpinus GWAS figure based on _S. alpinus_ genetic map
B. Sutherland
2017-07-20

### Overview
*Part 1*: Anchor anonymous markers from _Salvelinus alpinus_ onto the genetic map of _S. alpinus_. This will use the high-density genetic map of _S. alpinus_ (Nugent et al. 2017) in combination with anonymous markers from _S. alpinus_.   
    
*Part 2*: Combine the positioned markers with Fst values and plot in a GWAS figure.   

_Clone this repo, run all code from within the main repo_   

## 1. Anchoring anonymous markers   
### A. Obtain input Data
Put the following data into `02_data`    
From: [Figshare data](https://doi.org/10.6084/m9.figshare.5051821.v2) download:     
* Salp anon marker sequence file: `salp_tags.csv`    
From: [Supplemental Data from Nugent et al. 2017, G3](http://www.g3journal.org/content/7/2/543.supplemental), download the supplemental files:    
* Salp map marker file: `FileS1.xlsx`    
* Salp map position file: `FileS2.xlsx`    

Within `FileS1.xlsx`, save the sheet 'Map_SNPs' as a .csv entitled `FileS1.csv`    
Collect marker name and tag sequence from this file:
`grep -v 'Polymorphism' FileS1.csv | awk -F, '{ print $1 "," $4 }' > salp_marker_and_seq.csv`   

Within `FileS2.xlsx`, save the only sheet as a .csv file entitled `FileS2.csv`    

Collect only the lines with female linkage groups that have markers with positions:    
`awk -F, '{ print $1","$3","$4 }' FileS2.csv | sed 's/,AC-/,AC/g' | sed 's/,-/,empty/g'  | grep -vE 'NA|empty' - | grep -v 'Marker,Female,Map' - > ./salp_female_map.csv`

Same as above, but to collect male map:
`awk -F, '{ print $1","$2","$4 }' FileS2.csv | sed 's/,AC-/,AC/g' | sed 's/,-/,empty/g'  | grep -vE 'NA|empty|UNA' - | grep -v 'Marker,Male,Map' | sed 's/m\,/\,/g' | sed 's/\-\,/\,/g' > ./salp_male_map.csv`

Finally, to finish preparing the input data, go to R to make some final adjustments to prepare for `MapComp`. I suggest using Rstudio, setting working directory to this github repo.       
In addition to format adjusting, this will also change linkage groups AC-20 and AC-4 from the current format of split by _a_ and _b_ arms to one continuous linkage group with a cumulative cM position.   

Female
`sed 's/Salp/Salp.fem/g' salp_merged_sorted_clean.csv > salp_fem_sep_merged_sorted_clean.csv`

Male
`sed 's/Salp/Salp.male/g' salp_male_merged_sorted_clean.csv > salp_male_sep_merged_sorted_clean.csv`

Consensus
`cat salp_male_merged_sorted_clean.csv salp_fem_merged_sorted_clean.csv > consensus_merged_sorted_clean.csv`

Bring in Brook Charr map    
`cp /Users/wayne/Documents/bernatchez/01_Sfon_projects/12_JS_Salp_loci_rel_to_sex/salp_anon_to_sfon/02_data/sfon_markers.csv`


### B. Prepare Data For MapComp 

```
# Move to the data folder
cd 02_data

# To the anonymous markers, replace ‘alltags’ with ‘Salp.anon’, and the arbitrarily named variable LG ‘0’ to ‘1’ in sequence csv file to make compatible with MapComp    
sed 's/alltags/Salp.anon/g' salp_tags.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv

# Confirm information on MapComp input files     
wc -l salp.anon_markers.csv salp_merged_sorted_clean.csv
`1656 salp_merged_sorted_clean.csv` (mapped markers)   
`6230 salp.anon_markers.csv`   (anonymous markers)   

# Combine Sfon and Salp anonymous markers to make input for `MapComp`
`cat salp.anon_markers.csv salp_merged_sorted_clean.csv > salp.anon_salp.fem.map.csv` 

# Move out of the repo
cd ../../

```
Obtain MapComp iterative through the MapComp repo below     
`git clone https://github.com/enormandeau/mapcomp.git`

Move into the MapComp repo    
`cd mapcomp`    

Follow instructions given at the top of the MapComp iterative script:  
`01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`  

To obtain results as in the manuscript, run MapComp iterative with 10 iterations, with the default distance setting, and use the Atlantic Salmon reference genome as the genome intermediate:   
[ICSASG_v2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1)  
From: Lien et al., 2016. The Atlantic Salmon genome provides insights into rediploidization. Nature 533: 200–205.     

### B. Prepare and Run MapComp Iteratively
```
# Copy `salp.anon_markers.csv` to the `mapcomp/02_data` folder   
`cp ./../salp_anon_to_salp/02_data/salp.anon_salp.fem.map.csv ./02_data/markers.csv`

# Prepare the marker.csv file to a fasta file
`./01_scripts/00_prepare_input_fasta_file_from_csv.sh`

# Check the markers.fasta 
`wc -l 02_data/markers.fasta`
`15772 02_data/markers.fasta`

# Prepare MapComp variables and parameters
# Set the species name in the iterative mapping script
# e.g.  ANON=”Salp.anon”
vi ./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh

# Set the max distance in the mapcomp script (e.g. 1000000 or 10000000)
vi ./mapcomp

# Set the path to the genome file in both the following:   
vi ./mapcomp   
vi 01_scripts/01_bwa_align_reads.sh   

# Run MapComp iteratively 
./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh


# Collect results, this will be used by the GWAS script
awk '{ print $1","$5","$11 }' 03_mapped/pairings_out.txt > 03_mapped/Salp_mname_Salptotpos.csv

# Copy the result file 03_mapped/Salp_mname_Salptotpos.csv into the folder salp_anon_to_salp/02_data

### C. Combine with Fst and plot
# Open the GWAS script in R and follow instructions there   
`01_scripts/GWAS_from_MapComp_2016-11-02.R`    
# This script will merge the Fst and positional data, and plot using the information on the Brook Charr map    
# After running this script, you will have a GWAS figure with Fst by Brook Charr linkage group   
```
