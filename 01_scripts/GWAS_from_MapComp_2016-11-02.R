# Connect total position and Fst to plot in a GWAS-style plot
# B. Sutherland
# 2016-11-02

setwd("~/Documents/bernatchez/01_Sfon_projects/12_JS_Salp_loci_rel_to_sex/salp_anon_to_sfon/")

#rm(list=ls())

#### Data Import ####
# Set your input files file names
input.fst.txt = "02_data/Fst_6147SNPs_2016-11-11.txt"
input.pos.csv = "02_data/Salp_mname_Sfontotpos.csv"

# Load data
pos = read.csv(input.pos.csv, header = F, col.names = c("sp","mname","totpos"))
head(pos)
str(pos)

fst = read.table(file = input.fst.txt, header = T
                 , col.names = c("stacks.chr","pos","mname", "prob", "log_PO", "q_val", "alpha", "fst")
                 )
head(fst)
str(fst)

# check the concordance between the vectors
x <- fst$mname
y <- pos$mname
str(x)
str(y)
length(intersect(x,y)) # are all from Fst given positions (almost)

# Sort fst and pos by mname
fst.sorted <- fst[order(fst$mname),]
head(fst.sorted)

pos.sorted <- pos[order(pos$mname),]
head(pos.sorted)

# Merge
gwas.fst <- merge(fst.sorted, pos.sorted, by = "mname")
head(gwas.fst)
dim(gwas.fst)

#### Brook Charr LG info ####
# Obtain information of positions of LGs
library(qtl)
load(file = "02_data/LG_plot.RData")

# Collect the length of each chr
len.chr = NULL
for(i in names(sfon$geno)){
  len.chr = c(len.chr, max(sfon$geno[[i]]$map))  
}

len.chr

# Find cumulative length
cumul.leng <- cumsum(len.chr)
chr.end.df <- as.data.frame(cumul.leng)

# Create Positions
left.grey <- c(0, cumul.leng[c(FALSE, TRUE)]) # find x-axis positions, left side
right.grey <- c(cumul.leng[c(TRUE, FALSE)])   # find x-axis positions, right side


#### Plot ####
# Plot an empty plot
#par(mfrow=c(1,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
plot(x = c(0, max(cumul.leng)), y = c(0, max(gwas.fst$fst+0.025)), type = "n"
     , xaxt = 'n'
     , xlab = "Brook Char LG", ylab = "Fst",
     las = 1, cex.axis = 0.8)
# Add grey boxes for LGs
rect(xleft = left.grey[0:21], xright= right.grey[0:21],
     ybottom = 0, ytop = 11,
     col = c("lightgrey"),
     border=NA)

# Find position for labels
position <- ((right.grey[1:21] - left.grey[1:21])/2) + left.grey[1:21]

# Add axis and points
axis(side = 1, at = position, labels = seq(from = 1, to = 41, by = 2)
     , cex.axis = 0.8 )
points(gwas.fst$totpos, gwas.fst$fst, type = "p", cex = 0.7)

# save as 10 x 6


##### OBTAIN ONLY SIGNIFICANT FST MNAMES ####
outliers94 <- read.table(file = "../00_resources/94snps-outliers_2016-11-11.txt", header = T)
head(outliers94)
names(outliers94) #"SNP" is the mname

outliers94$SNP
significantSNPs <- outliers94$SNP

gwas.fst.sig <- gwas.fst[gwas.fst$mname %in% significantSNPs, ] # select only the significant snps

points(gwas.fst.sig$totpos, gwas.fst.sig$fst, pch = 3, cex = 1)


# export markers with positions
write.csv(gwas.fst$mname, file = "salp_markers_with_Sfon_positions_2016-11-13.csv")



##### DATA EXPLORATION #####
# More detail on markers of interest
gwas.fst$totpos[which(gwas.fst$fst > 0.03)] # positions
gwas.fst$mname[which(gwas.fst$fst > 0.03)] # marker names

# find markers in BC38 with Fst > 0.1
chr.end.df
chr.start <- 5125.2076
chr.end <- 5284.0927
gwas.fst$mname[which(gwas.fst$fst > 0.005 & 
                       gwas.fst$totpos > chr.start & 
                       gwas.fst$totpos < chr.end)]

gwas.fst$totpos[gwas.fst$mname==298178] #find out where a marker is
