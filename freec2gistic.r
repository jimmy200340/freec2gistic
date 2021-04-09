library(tidyverse)
library(GenomicRanges)
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
	make_option("--work-dir", type = "character", help = "R dir")
)

opt <- parse_args(OptionParser(option_list=option_list))

setwd(opt$'work-dir')
files <- list.files(path = "./", pattern="freec_segments.bed", full.names=TRUE)

segData <- NULL
for( i in files ) {
  print(i)
  tab <- read.table(i)
  id <- gsub(pattern = ".//", replacement = "", i)
  id <- gsub(pattern = ".freec_segments.bed", replacement = "", id)
  tab$Sample <- id
  tab$Dummy <- 'NA'
  tab <- tab[, c("Sample", "V1", "V2", "V3", "Dummy", "V4")]
  # convert ratio to log10(ratio)-1
  tab$V4 <- log2(tab$V4) - 1
  segData <- rbind(segData, tab)
}

segData <- segData %>% 
  dplyr::filter( V1 != 'chrX' & V1 != 'chrY') %>% 
  dplyr::rename( Chromosome = V1, Start = V2, End = V3, log2ratio = V4 )

segData$Chromosome <- gsub(pattern = "chr", replacement = "", segData$Chromosome)

######
common.cnvs <- GRanges( read.delim( file = "/home/jimmy200340/opt/FREEC-11.6/CNV.hg19.bypos.111213.txt" ) )
# get index of all CNVs without a hit in common.cnv
#subsetByOverlaps( GRanges(cnvdata[4, 1:4]), common.cnvs ) 
idx <- which( countOverlaps( GRanges(segData[, -1]), common.cnvs ) == 0 )
segData <- segData[idx, ]
rownames(segData) <- NULL

######

files <- list.files(path = "./", pattern="bam_sample.cpn", full.names=TRUE)
cpnData <- NULL
for( i in files ) {
  print(i)
  tab <- read.table(i)
  id <- gsub(pattern = ".//", replacement = "", i)
  id <- gsub(pattern = ".mark.sort.recal.bam_sample.cpn", replacement = "", id)
  tab$Sample <- id
  tab <- tab[, c("Sample", "V1", "V2", "V3", "V4", "V5")]
  cpnData <- rbind(cpnData, tab)
}

colnames(cpnData) <- c("Sample", "Chromosome", "Start", "End", "Number", "Range")

segData_num <- NULL
for ( i in unique( cpnData$Sample ) ) {
  print(i)
  cpn.s <- cpnData %>% dplyr::filter( Sample == i ) %>% dplyr::select( -Sample ) %>% GRanges
  seg.s <- segData %>% dplyr::filter( Sample == i ) %>% dplyr::select( -Sample ) %>% dplyr::rename(Num_Probes=Dummy) %>% GRanges
  
  # sum features per interval
  for ( j in 1:nrow( data.frame(seg.s) ) ) {
    Num_Probes <- sum( data.frame( cpn.s[ c(data.frame( findOverlaps(seg.s[j,], cpn.s) )$subjectHits), ] )$Number )
    mcols( seg.s[j, ] )$Num_Probes <- Num_Probes
  }
  
  seg.s <- data.frame(seg.s)
  seg.s$Sample <- i
  segData_num <- rbind(segData_num, seg.s)
  
}

segData_num <- segData_num %>% 
dplyr::select( Sample, seqnames, start, end, Num_Probes, log2ratio)

write.table(x = segData, file = "gistic.segments.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(x = segData_num, file = "gistic.segments.numProbes.txt", sep = "\t", row.names = F, col.names = F, quote = F)