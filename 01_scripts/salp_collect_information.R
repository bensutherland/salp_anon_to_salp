# R script to clean up the input files to be used in MapComp (iterative)
setwd("~/Documents/bernatchez/01_Sfon_projects/14_JS_migration/salp_anon_to_salp")

# used to combine the following files:
marker.file <- "02_data/salp_female_map.csv" # markers
seq.file <- "02_data/salp_marker_and_seq.csv" # seq

# import seq file
seq <- read.csv(file = seq.file, header = F, col.names = c("mname","seq"))
head(seq)
str(seq)

# import marker file
markers <- read.csv(file = marker.file, header = F, col.names = c("mname","LG","pos"))
head(markers)
str(markers)

# Match marker names between the two files:
length(intersect(markers$mname, seq$mname)) #1656 matching

# fix the LG issues (AC20 and AC4) in the markers file
# find all rows with AC20b
markers[which(markers$LG == "AC20a"),] 
max(markers[which(markers$LG == "AC20a"),3])
# total length of AC20a is 108.2 cM, so add this to AC20b

which(markers$LG == "AC20b")
markers[which(markers$LG == "AC20b"),3] + 108.2

# adds the first part of LG to the second (cumulative)
# AC20
add.cumul <- max(markers[which(markers$LG == "AC20a"),3]) # number to add
markers[which(markers$LG == "AC20b"),3]
markers[which(markers$LG == "AC20b"),3] <- markers[which(markers$LG == "AC20b"),3] + add.cumul
markers[which(markers$LG == "AC20b"),3]

# Change names from 20b and 20a to 20 now that cumulative positions have been given
markers$LG <- gsub(pattern = "20b", replacement = "20", x = markers$LG)
markers$LG <- gsub(pattern = "20a", replacement = "20", x = markers$LG)
unique(markers$LG)

# AC4
add.cumul <- max(markers[which(markers$LG == "AC4p"),3])
markers[which(markers$LG == "AC4q"),3]
markers[which(markers$LG == "AC4q"),3] <- markers[which(markers$LG == "AC4q"),3] + add.cumul
markers[which(markers$LG == "AC4q"),3]

markers$LG <- gsub(pattern = "4p", replacement = "4", x = markers$LG)
markers$LG <- gsub(pattern = "4q", replacement = "4", x = markers$LG)
unique(markers$LG)

## Change LG to numeric
markers$LG <- gsub("[[:punct:]]", "", markers$LG) #remove plus sign
markers$LG <- gsub(pattern = "121", replacement = "1", x = markers$LG) # remove extra 21
markers$LG <- as.numeric(gsub(pattern = "AC", replacement = "", x = markers$LG))
unique(markers$LG)
# Note: this makes so that AC1+21 is now only labeled as AC1

## Sort marker object
markers.sorted <- markers[order(markers$LG, markers$pos),] # double check
unique(markers.sorted$LG)

## Merge marker objects
salp <- merge(x = markers.sorted, y = seq, by = "mname")
str(salp)
# double check

# Add extra columns
salp$species <- rep(x = "Salp", times = length(salp[,1]))
salp$cM.null <- rep(x = 0, times = length(salp[,1]))

head(salp)

# Reorder columns
salp <- salp[,c(5,2,3,6,1,4)]
head(salp)

# sort again
salp <- salp[order(salp$LG, salp$pos),]
head(salp)

# Write out to use in MapComp iterative
write.table(salp, file = "02_data/salp_merged_sorted_clean.csv"
            , row.names = F, col.names = F, sep = ","
            , quote = F)

