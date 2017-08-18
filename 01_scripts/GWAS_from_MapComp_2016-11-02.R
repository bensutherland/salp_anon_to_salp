# Connect total position and Fst to plot in a GWAS-style plot
# B. Sutherland
# 2017-08-10

setwd("~/Documents/bernatchez/01_Sfon_projects/14_JS_migration/salp_anon_to_salp/")

#rm(list=ls())

#### Data Import ####
# Set your input files file names
input.fst.txt = "02_data/batch_1_miss25MigrantsRepAcousHistSexRem_EKASURcomb_bayescan_9jun2016_fst.txt"
input.pos.csv = "02_data/Salp_mname_Salptotpos.csv"
index.file = "02_data/tagID_sequences_wIndex_6136SNPsSexRem_30Jun2016.txt"

# Load data
pos = read.csv(input.pos.csv, header = F, col.names = c("sp","mname","totpos"))
head(pos)
str(pos)

# Import various test statistics
bayescan.fst <- read.table(file = "02_data/batch_1_miss25MigrantsRepAcousHistSexRem_EKASURcomb_bayescan_9jun2016_fst.txt"
                           , header = T
                           , col.names = c("Index", "BS.prob", "BS.log10.PO", "BS.qval", "BS.alpha", "BS.fst"))
# fix bayescan zero values, set them to the next lowest decimal place possible (five sig figs)
bayescan.fst$BS.qval[bayescan.fst$BS.qval < 0.00001] <- 0.00001

PCAdapt.fst <- read.table(file = "02_data/PCAdapt_qvaluesALLsnps_K2_imputed_14sep2016.txt"
                           , header = T
                           , col.names = c("PA.chr","PA.pos","PA.pval","PA.qval","Index"))

LFMM.fst <- read.table(file = "02_data/lfmm_adjPvalues_imputed_14Sep2016.txt"
                          , header = T
                          , col.names = c("LF.adj.pval","LF.Z","Index"))

colnames(bayescan.fst)
head(bayescan.fst, n = 2)
colnames(PCAdapt.fst)
head(PCAdapt.fst, n = 2)
colnames(LFMM.fst)
head(LFMM.fst, n = 2)

# Join the files together by column 'Index'
bayescan.PCAdapt <- merge(x = bayescan.fst, y = PCAdapt.fst, by = "Index")
all.test.stats <- merge(x = bayescan.PCAdapt, y = LFMM.fst, by = "Index")
colnames(all.test.stats)
str(all.test.stats)
fst <- all.test.stats


# fst = read.table(file = input.fst.txt, header = T
#                  #, col.names = c("Index","prob","log10.PO", "qval", "alpha", "fst")
#                  )
# colnames(fst)
# head(fst)
# str(fst)

index <- read.table(file = index.file, header = T
                    #col.names = c("ID")
                    )
head(index)
str(index)

# Attach ID column to fst file from index file
colnames(index)
colnames(fst)

# merge index and fst files
# check concordance first:
length(intersect(fst$Index, index$Index))
fst.indexed <- merge(fst, index, by = "Index")

dim(fst)
dim(index)
dim(fst.indexed)

colnames(fst.indexed)

# merge fst and position file
# check the concordance between the vectors
x <- fst.indexed$ID
y <- pos$mname
str(x)
str(y)
length(intersect(x,y)) # are all from Fst given positions (note: almost)

# Sort fst and pos by mname
fst.sorted <- fst.indexed[order(fst.indexed$ID),]
head(fst.sorted)

# rename ID to mname for consistency
colnames(fst.sorted)[which(colnames(fst.sorted)=="ID")] <- "mname"

pos.sorted <- pos[order(pos$mname),]
head(pos.sorted)

# Merge
gwas.fst <- merge(fst.sorted, pos.sorted, by = "mname")
head(gwas.fst)
dim(gwas.fst)


###### MAP #######
# Plot the image using the map coordinates
# Find map coordinates for image (Arctic Charr)
# bring in map file:
map <- read.csv(file = "02_data/consensus_merged_sorted_clean.csv", header = F
                , col.names = c("species","LG","pos","blank","marker","seq"))
head(map)

# piece that works:
#tail(map$pos[which(map$LG == 3)], n = 1)

# Collect the length of each chr
len.chr = NULL
for(i in 1:37){
  len.chr = c(len.chr, 
              tail(map$pos[which(map$LG == i)], n = 1))  
}

len.chr

# Find cumulative length
cumul.leng <- cumsum(len.chr)
chr.end.df <- as.data.frame(cumul.leng)

# Create Positions
left.grey <- c(0, cumul.leng[c(FALSE, TRUE)]) # find x-axis positions, left side
right.grey <- c(cumul.leng[c(TRUE, FALSE)])   # find x-axis positions, right side

# Find position for labels
position <- ((right.grey[1:18] - left.grey[1:18])/2) + left.grey[1:18]

#### Plot SINGLE FST ####
# Plot an empty plot
#par(mfrow=c(1,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
plot(x = c(0, max(cumul.leng)), y = c(0, max(gwas.fst$fst+0.025)), type = "n"
     , xaxt = 'n'
     , xlab = "Arctic Char Linkage Group", ylab = "Fst",
     las = 1, cex.axis = 0.8)
# Add grey boxes for LGs
rect(xleft = left.grey[0:18], xright= right.grey[0:19],
     ybottom = 0, ytop = 11,
     col = c("lightgrey"),
     border=NA)

# Find position for labels
position <- ((right.grey[1:18] - left.grey[1:18])/2) + left.grey[1:18]

# Add axis and points
labels.no21 <- c(seq(from = 1, to = 20), seq(from = 22, to = 37))
labels.to.add <- labels.no21[c(TRUE, FALSE)]
axis(side = 1, at = position, labels = labels.to.add
     , cex.axis = 0.8 )
points(gwas.fst$totpos, gwas.fst$fst, type = "p", cex = 0.8)

# save as 10 x 6



##### PLOT ALL THREE TEST STATS ####
# Plot an empty plot
#par(mfrow=c(1,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

# these are the characteristics we want to plot:
# plot.statistics <- as.data.frame(c(
#   #  -log10(gwas.fst$BS.qval)
#     gwas.fst$BS.fst
#   , -log10(gwas.fst$PA.qval)
#   , -log10(gwas.fst$LF.adj.pval))
# )


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

# find max position for each test statistic for the y-axis of the plot
max.pos <- NULL
for(i in 1:3){
  max.pos = c(max.pos, max(neglog10.plot.stats[,i], na.rm = T)
              )  
}
max.pos

max.pos.buffer <- max.pos+max.pos*0.1 # add 10 % to the max position

#yaxis.titles <- c("bayescan.fst","PCAdapt -log10(qval)", "LFMM -log10(adj pval)")
yaxis.titles <- c("bayescan -log10(qval)","PCAdapt -log10(qval)", "LFMM -log10(adj pval)")
low.sig.thresh <- -log10(0.05)
high.sig.thresh <- -log10(0.01)

## Plot all three
par(mfrow=c(3,1), mar= c(4,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

for(i in 1:3){
  plot(x = c(0, max(cumul.leng))
       #, y = c(0, max.pos.buffer[i]) # dynamic position
       , y = c(0, 6.5) # fixed position
       , type = "n"
       , xaxt = 'n'
       , xlab = "Arctic Char Linkage Group"
       , ylab = yaxis.titles[i],
       las = 1, cex.axis = 0.8)
  # Add grey boxes for LGs
  rect(xleft = left.grey[0:18], xright= right.grey[0:19],
       ybottom = 0, ytop = 11,
       col = c("lightgrey"),
       border=NA)
  # Add LG labels
  labels.no21 <- c(seq(from = 1, to = 20), seq(from = 22, to = 37))
  labels.to.add <- labels.no21[c(TRUE, FALSE)]
  axis(side = 1, at = position, labels = labels.to.add
       , cex.axis = 0.8 )
  # Plot points
  points(gwas.fst$totpos, neglog10.plot.stats[,i], type = "p", cex = 0.8)
  
  # Plot APPROX. signif thresholds
  abline(h = low.sig.thresh, lty = 2)
  abline(h = high.sig.thresh, lty = 3)
}

# save as 10 x 6





###### OLD MATERIAL #####
# ##### OBTAIN ONLY SIGNIFICANT FST MNAMES ####
# outliers94 <- read.table(file = "../00_resources/94snps-outliers_2016-11-11.txt", header = T)
# head(outliers94)
# names(outliers94) #"SNP" is the mname
# 
# outliers94$SNP
# significantSNPs <- outliers94$SNP
# 
# gwas.fst.sig <- gwas.fst[gwas.fst$mname %in% significantSNPs, ] # select only the significant snps
# length(gwas.fst.sig$mname) # how many outliers were present in the positioned markers?
# 
# points(gwas.fst.sig$totpos, gwas.fst.sig$fst, pch = 3, cex = 1)
# 
# 
# # export markers with positions
# write.csv(gwas.fst$mname, file = "salp_markers_with_Sfon_positions_2016-11-13.csv")

gwas.fst$mname[which(gwas.fst$BS.fst)]

##### DATA EXPLORATION #####
# More detail on markers of interest
gwas.fst$totpos[which(gwas.fst$fst > 0.03)] # positions
gwas.fst$mname[which(gwas.fst$fst > 0.03)] # marker names

# find markers in BC38 with Fst > 0.1
chr.end.df
chr.start <- 1721.3238
chr.end <- 1977.5399
gwas.fst$mname[which(gwas.fst$fst > 0.005 & 
                       gwas.fst$totpos > chr.start & 
                       gwas.fst$totpos < chr.end)]

gwas.fst.sig$mname[which(gwas.fst.sig$totpos > chr.start & gwas.fst.sig$totpos < chr.end)]

gwas.fst$totpos[gwas.fst$mname==298178] #find out where a marker is
