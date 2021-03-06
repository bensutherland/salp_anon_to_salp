# Connect total position and Fst to plot in a GWAS-style plot
# B. Sutherland
# 2017-08-20

setwd("~/Documents/bernatchez/01_Sfon_projects/14_JS_migration/salp_anon_to_salp/")

#rm(list=ls())

#### Data Import ####
# Set your input files file names
index.file = "02_data/tagID_sequences_wIndex_6136SNPsSexRem_30Jun2016.txt" # index
input.pos.csv = "02_data/Salp_mname_Salptotpos_consensus.csv" # position data
res.bayescan = "02_data/batch_1_miss25MigrantsRepAcousHistSexRem_EKASURcomb_bayescan_9jun2016_fst.txt" 
res.PCAdapt = "02_data/PCAdapt_qvaluesALLsnps_K2_imputed_14sep2016.txt"
res.LFMM = "02_data/lfmm_adjPvalues_imputed_14Sep2016.txt"


# Load files
index <- read.table(file = index.file, header = T)

pos <- read.csv(input.pos.csv, header = F, col.names = c("sp","mname","totpos"))
head(pos)
str(pos)
max(pos$totpos) # finds totpos of last marker in map

bayescan.fst <- read.table(file = res.bayescan
                           , header = T
                           , col.names = c("Index", "BS.prob", "BS.log10.PO", "BS.qval", "BS.alpha", "BS.fst"))

PCAdapt.fst <- read.table(file = res.PCAdapt
                           , header = T
                           , col.names = c("PA.chr","PA.pos","PA.pval","PA.qval","Index"))

LFMM.fst <- read.table(file = res.LFMM
                          , header = T
                          , col.names = c("LF.adj.pval","LF.Z","Index"))

# Transform bayescan 0 values to 0.00001 (sig fig basement of bayescan)
bayescan.fst$BS.qval[bayescan.fst$BS.qval < 0.00001] <- 0.00001


# Merge data
# Merge test statistic dataframes together by 'Index'
bayescan.PCAdapt <- merge(x = bayescan.fst, y = PCAdapt.fst, by = "Index")
all.test.stats <- merge(x = bayescan.PCAdapt, y = LFMM.fst, by = "Index")
colnames(all.test.stats)
str(all.test.stats)

# Rename main dataframe
fst <- all.test.stats

# Merge test statistic data with index file that contains tag 'ID' using index column
colnames(index)
colnames(fst)
length(intersect(fst$Index, index$Index)) # check concordance first

fst.indexed <- merge(fst, index, by = "Index") # merge
colnames(fst.indexed)

# Check that the two have been fully added together
dim(fst)
dim(index)
dim(fst.indexed)

# Merge indexed test statistics file with position file
# First check correspondence
x <- fst.indexed$ID
y <- pos$mname
str(x)
str(y)
length(intersect(x,y)) 
# Note: not all positioned data have corresponding values in test statistic,
# because the sex-linked SNPs were removed from the analysis (n = 94 tags)

# Second, sort both files to be merged
# Sort data file by ID
fst.sorted <- fst.indexed[order(fst.indexed$ID),]
head(fst.sorted)
# rename ID to mname for consistency with position file
colnames(fst.sorted)[which(colnames(fst.sorted)=="ID")] <- "mname"

# Sort pos file by mname (aka ID)
pos.sorted <- pos[order(pos$mname),]
head(pos.sorted)

# Finally, merge test statistic file with position file
gwas.fst <- merge(fst.sorted, pos.sorted, by = "mname")
head(gwas.fst)
dim(gwas.fst)


#### Map Information Import ####
# To plot stats using map coordinates
map.file = "02_data/consensus_merged_sorted_clean.csv"

# Import
map <- read.csv(file = map.file, header = F
                , col.names = c("species","LG","pos","blank","marker","seq"))
head(map) # note is not ordered (for LG10)

# Sort map file
map.sorted <- map[order(map$LG, map$pos),]


# piece that works:
#tail(map$pos[which(map$LG == 3)], n = 1)

# Collect the length of each chr (Note: requires that map is sorted)
len.chr = NULL
for(i in 1:37){
  len.chr = c(len.chr, 
              tail(map.sorted$pos[which(map.sorted$LG == i)], n = 1))  
}

len.chr

# Find cumulative length of each chromosome adding the one before
cumul.leng <- cumsum(len.chr)
chr.end.df <- as.data.frame(cumul.leng)

# Create positions of grey boxes for plotting
left.grey <- c(0, cumul.leng[c(FALSE, TRUE)]) # find x-axis positions, left side
right.grey <- c(cumul.leng[c(TRUE, FALSE)])   # find x-axis positions, right side

# Find position for labels
position <- ((right.grey[1:18] - left.grey[1:18])/2) + left.grey[1:18]


##### Plot the three test statistics ####
colnames(gwas.fst)
plot.statistics <- gwas.fst[,c(#7 # currently BS.fst
                               5 # BS.qval
                               ,11 # PA.qval
                               ,12 # LF.pval
                               )]
head(plot.statistics)
#neglog10.plot.stats <- cbind(plot.statistics[,1], -log10(plot.statistics)[,c(2,3)])
neglog10.plot.stats <- -log10(plot.statistics)
head(neglog10.plot.stats)


# Find max position for each test statistic for the y-axis of the plot
max.pos <- NULL
for(i in 1:3){
  max.pos = c(max.pos, max(neglog10.plot.stats[,i], na.rm = T)
              )  
}
max.pos

max.pos.buffer <- max.pos+max.pos*0.1 # add 10 % to the max position

# Information for plot
yaxis.titles <- c("bayescan -log10(qval)","PCAdapt -log10(qval)", "LFMM -log10(adj pval)")
low.sig.thresh <- -log10(0.05)
high.sig.thresh <- -log10(0.01)


## Plot
par(mfrow=c(3,1), mar= c(4,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

for(i in 1:3){
  plot(x = c(0, max(pos$totpos))
       #, y = c(0, max.pos.buffer[i]) # dynamic position
       , y = c(0, 6.5) # fixed position
       , type = "n" # plot empty chart, will fill in with points later
       , xaxt = 'n'
       , xlab = "Arctic Char Linkage Group"
       , ylab = yaxis.titles[i],
       las = 1, cex.axis = 0.8)
  # Add grey boxes for LGs
  rect(xleft = left.grey[0:18], xright= right.grey[0:19],
       ybottom = 0, ytop = 6.5,
       col = c("lightgrey"),
       border=NA)
  # Add LG labels
  labels.no21 <- c(seq(from = 1, to = 20), seq(from = 22, to = 37))
  labels.to.add <- labels.no21[c(TRUE, FALSE)]
  axis(side = 1, at = position, labels = labels.to.add
       , cex.axis = 0.8 )
  # Plot points
  points(gwas.fst$totpos, neglog10.plot.stats[,i], type = "p", cex = 0.8)
  
  # Plot signif thresholds
  abline(h = low.sig.thresh, lty = 2)
  abline(h = high.sig.thresh, lty = 3)
}

# save as 10 x 6



##### Data exploration (not required) #####
# The rest is for data exploration only
# More detail on markers of interest

# gwas.fst$totpos[which(gwas.fst$fst > 0.03)] # positions
# gwas.fst$mname[which(gwas.fst$fst > 0.03)] # marker names

# Find markers by LG
chr.end.df

# Create dataframe
LG.names <- c(1:20, 22:37)
chr.info <- cbind(chr.end.df,
               as.numeric(LG.names))
colnames(chr.info) <- c("cumul.length", "LG.names")


# Select numeric LG (e.g. 30)
COI <- 34

# Run following to identify the position range
chr.start <- NULL ; chr.end <- NULL
chr.start <- chr.info$cumul.length[chr.info$LG.names == COI-1] + 0.1
chr.end <- chr.info$cumul.length[chr.info$LG.names == COI]

# Identify the marker
gwas.fst$mname[which(
                       -log10(gwas.fst$BS.qval) > 1.3 & 
                       gwas.fst$totpos > chr.start & 
                       gwas.fst$totpos < chr.end)]

# gwas.fst.sig$mname[which(gwas.fst.sig$totpos > chr.start & gwas.fst.sig$totpos < chr.end)]
# 
# gwas.fst$totpos[gwas.fst$mname==298178] #find out where a marker is

