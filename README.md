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
Salp sequence file: `salp_tags.csv`
Sfon MapComp file: `Sfon_v4.3_female_map.csv`

Make sure to add in python script needed to remove markers
cp /Users/wayne/Desktop/Scripts_Eric/fasta_remove.py ./01_scripts/

**Pipeline:**
'mapcomp_iterative'

Methods:
1. Move the input data into the folder 02_data and then move into this folder

2. Replace ‘alltags’ with ‘Salp.anon’, and LG ‘0’ to ‘1’ in sequence csv file
```
sed 's/alltags/Salp.anon/g' salp_tags.csv | sed 's/anon,0/anon,1/g' > salp.anon_markers.csv
```

3. Remove ‘>’ and header from Sfon csv file, and add the 0 for the totpos position
```sed 's/>//g' Sfon_v4.3_female_map.csv | grep -vE '^species' | awk -F, '{ print $1","$2","$3","0","$4","$5","$6 }' > sfon_markers.csv```

4. Confirm information on Sfon and Salp input files
```wc -l *markers.csv``` 
    `6230 salp.anon_markers.csv`
    `3826 sfon_markers.csv`

5. Combine Sfon and Salp anonymous markers to make input for `MapComp`
```cat salp.anon_markers.csv sfon_markers.csv > salp.anon_sfon_markers.csv```

6. Copy `salp.anon_markers.csv` to the `mapcomp_iterative` folder and move to the main directory for mapcomp_iterative
```cp salp.anon_sfon_markers.csv ./../../mapcomp_iterative/02_data/```
```cd ./../../mapcomp_iterative```

7. Prepare the marker.csv file to a .fasta file
```./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/salp.anon_sfon_markers.csv``` 

8. Check the markers.fasta 
```wc -l 02_data/markers.fasta```
`20112 02_data/markers.fasta`

9. Prepare MapComp_Iterative 
Set the species name in the iterative mapping script
e.g.  ANON=”Salp.anon”
```vi ./01_scripts/remove_paired_anon_and_pair_again.sh```
Set the max distance in the mapcomp script (1000000)
```vi ./01_scripts/mapcomp```

10. Run MapComp iteratively 
./01_scripts/remove_paired_anon_and_pair_again.sh

11. Collect results, this will be used by the GWAS script
```awk '{ print $1","$5","$11 }' 03_mapped/pairings_out.txt > 03_mapped/Salp_mname_Sfontotpos.csv```

