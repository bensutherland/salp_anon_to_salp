## Arctic Charr anonymous marker placement on genetic map
B. Sutherland    
2017-07-20    
This methods description accompanies the results presented in the following manuscript:    
Moore J.-S., Harris L. N., Le Luyer J., Sutherland B. J. G., Rougemont Q., Tallman R. F., Fisk A. T., Bernatchez L., 2017 Migration harshness drives habitat choice and local adaptation in anadromous Arctic Char: evidence from integrating population genomics and acoustic telemetry. bioRxiv: 1–39.

### Overview
**Part 1** will anchor anonymous markers onto a high-density genetic map of _Salvelinus alpinus_ (Nugent et al. 2017)   
**Part 2** will combine positioned markers with population genetic values (e.g. Fst) and plot in a GWAS figure   

_Clone this repo, run all code from within the main repo_   

## Part 1: Position anonymous markers   
### A) Obtain input data and put in 02_data
From: [Figshare data](https://doi.org/10.6084/m9.figshare.5051821.v2)       
* Salp anon marker sequence file: `salp_tags.csv`    

From: [Supplemental Data from Nugent et al. 2017, G3](http://www.g3journal.org/content/7/2/543.supplemental)    
* Salp map marker file: `FileS1.xlsx`    
* Salp map position file: `FileS2.xlsx`    

i) Within `FileS1.xlsx`, save the sheet 'Map_SNPs' as a .csv entitled `FileS1.csv`    

ii) Collect marker name and sequence from this file:    
`grep -v 'Polymorphism' FileS1.csv | awk -F, '{ print $1 "," $4 }' > salp_marker_and_seq.csv`   

iii) Within `FileS2.xlsx`, save the only sheet as a .csv file entitled `FileS2.csv`    

iv) Collect only the lines with female linkage groups that have markers with positions (and a little formatting):    
`awk -F, '{ print $1","$3","$4 }' FileS2.csv | sed 's/,AC-/,AC/g' | sed 's/,-/,empty/g'  | grep -vE 'NA|empty' - | grep -v 'Marker,Female,Map' - > ./salp_female_map.csv`

v) Same as above, but collect the male map:
`awk -F, '{ print $1","$2","$4 }' FileS2.csv | sed 's/,AC-/,AC/g' | sed 's/,-/,empty/g'  | grep -vE 'NA|empty|UNA' - | grep -v 'Marker,Male,Map' | sed 's/m\,/\,/g' | sed 's/\-\,/\,/g' > ./salp_male_map.csv`

vi) Finish preparing the data using R (i.e. formats, make AC-20 and AC-4 continuous naming and cumulative position instead of LG arm split).   
I suggest opening the following script in RStudio, setting working directory to this repo.       
`01_scripts/salp_collect_information.R`   
This will produce salp_male_merged_sorted_clean.csv and salp_merged_sorted_clean.csv

Female
`sed 's/Salp/Salp.fem/g' salp_merged_sorted_clean.csv > salp_fem_sep_merged_sorted_clean.csv`

Male
`sed 's/Salp/Salp.male/g' salp_male_merged_sorted_clean.csv > salp_male_sep_merged_sorted_clean.csv`

Consensus
`cat salp_male_merged_sorted_clean.csv salp_fem_merged_sorted_clean.csv > consensus_merged_sorted_clean.csv`


### B) Format data for MapComp 

i) Move to the data folder, and replace anonymous title ‘alltags’ with ‘Salp.anon’. Also for anonymous markers make the LG 0 variable all equal LG 1.   
`sed 's/alltags/Salp.anon/g' salp_tags.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv`   

ii) Confirm information on MapComp input files     
wc -l salp.anon_markers.csv consensus_merged_sorted_clean.csv   
`3145 consensus_merged_sorted_clean.csv` (mapped markers)   
`6230 salp.anon_markers.csv`   (anonymous markers)   

iii) Combine all markers to make input for MapComp    
`cat salp.anon_markers.csv consensus_merged_sorted_clean.csv > salp.anon_salp.fem.map.csv`     

iv) Move back up to the main folder, and clone in the MapComp repo   
`cd ../../`    
`git clone https://github.com/enormandeau/mapcomp.git`    

v) Move into the MapComp repo   
`cd mapcomp`    

vi) Follow instructions given at the top of the MapComp iterative script:  
`01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`  

To obtain results as in the manuscript, run MapComp iterative with 10 iterations, with the default distance setting, and use the Atlantic Salmon reference genome as the genome intermediate:   
[ICSASG_v2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1)  
From: Lien et al., 2016. The Atlantic Salmon genome provides insights into rediploidization. Nature 533: 200–205.     

### C. Prepare and run MapComp iteratively
i) Copy the combined output from above into the `mapcomp/02_data` folder   
`cp ./../salp_anon_to_salp/02_data/<combined_data_here> ./02_data/markers.csv`

ii) Prepare the marker.csv file to a fasta file
`./01_scripts/00_prepare_input_fasta_file_from_csv.sh`

iii) Check the markers.fasta 
`grep -c '>' 02_data/markers.fasta`
`9375 02_data/markers.fasta`

iv) Prepare MapComp variables and parameters
Using vi, or similar, set the species name in the iterative mapping script   
e.g.  `ANON=”Salp.anon”`    
`vi ./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`

Set the max distance in the mapcomp script (e.g. 1000000 or 10000000)
`vi ./mapcomp`

Set the path to the genome file in both the following:   
`vi ./mapcomp`   
`vi 01_scripts/01_bwa_align_reads.sh`   

v) Run MapComp iteratively 
`./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh`

vi) Collect results, this will be used in Part 2.    
`awk '{ print $1","$5","$11 }' 03_mapped/pairings_out.txt > 03_mapped/Salp_mname_Salptotpos.csv`

vii) Copy the result file from the previous step into the folder `salp_anon_to_salp/02_data`   


## Part 2: Combine with population genetic values and plot     
i) Open the script entitled `01_scripts/GWAS_from_MapComp_2016-11-02.R` in R and follow instructions there. This will produce Figure 7 from the associated manuscript.       
