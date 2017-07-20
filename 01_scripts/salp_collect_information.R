
setwd(dir = "/Users/wayne/Documents/bernatchez/05_continuing_salmonid_MapComp/Salpinus_supp_files_from_nugent_etal2016")



# used to combine the following files:
marker.file <- "salp_female_map.csv" # markers
seq.file <- "salp_marker_and_seq.csv" # seq

# import seq file
seq <- read.csv(file = seq.file, header = F, col.names = c("mname","seq"))
head(seq)
str(seq)


## MARKERS ##
# import marker file
markers <- read.csv(file = marker.file, header = F, col.names = c("mname","LG","pos"))
head(markers)
str(markers)

# number matching between the two:
length(intersect(markers$mname, seq$mname)) #1656

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

markers$LG <- gsub(pattern = "20b", replacement = "20", x = markers$LG)
markers$LG <- gsub(pattern = "20a", replacement = "20", x = markers$LG)
unique(markers$LG)

#AC4
add.cumul <- max(markers[which(markers$LG == "AC4p"),3])
markers[which(markers$LG == "AC4q"),3]
markers[which(markers$LG == "AC4q"),3] <- markers[which(markers$LG == "AC4q"),3] + add.cumul
markers[which(markers$LG == "AC4q"),3]

markers$LG <- gsub(pattern = "4p", replacement = "4", x = markers$LG)
markers$LG <- gsub(pattern = "4q", replacement = "4", x = markers$LG)
unique(markers$LG)

# change LG to numeric
markers$LG <- gsub("[[:punct:]]", "", markers$LG) #remove plus sign
markers$LG <- gsub(pattern = "121", replacement = "1", x = markers$LG) # remove extra 21
markers$LG <- as.numeric(gsub(pattern = "AC", replacement = "", x = markers$LG))
unique(markers$LG)

# sort marker object
markers.sorted <- markers[order(markers$LG, markers$pos),] # double check
unique(markers.sorted$LG)

#merge
salp <- merge(x = markers.sorted, y = seq, by = "mname")
str(salp)

# add extra columns
salp$species <- rep(x = "Salp", times = length(salp[,1]))
salp$cM.null <- rep(x = 0, times = length(salp[,1]))

head(salp)

# reorder columns
salp <- salp[,c(5,2,3,6,1,4)]
head(salp)

# sort again
salp <- salp[order(salp$LG, salp$pos),]
head(salp)

write.table(salp, file = "salp_merged_sorted_clean.csv", row.names = F, col.names = F, sep = ",")



