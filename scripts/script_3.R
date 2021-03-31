#!/usr/bin/env Rscript

library("icesTAF")
library("DescTools") 
library("GenomicRanges")
library("tidyr")
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#args[1] path to the input and the output files (i.e., "/home/")

##Remove redundancy and assign an ID for each TE and each method from the the bed files (REPET, TIDAL, TEMP)
#Read bed file from REPET info
parent.o <- args[1]
x <- read.table(file=paste(parent.o, "11_strains_repet_concatenated_indexed_sorted_bedtools.txt", sep=""), header=F, sep="\t")

#Count the number of strains in which each TE is present
a <- str_count(x[,4],",")
b <- x
b[,4] <- as.character(b[,4])
x[,4] <- as.character(x[,4])
f <- NULL
for (i in 1:length(a)) {
	if (a[i] > 11) {
		b[i,4] <- as.character(paste(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T))), collapse=","))
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} else {
	 	b[i,4] <- x[i,4]
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} 
}
 
b.f <- cbind(b[ ,1:3], f, b[ ,4:5])
colnames(b.f) <-  c("Chr_TE_ID", "Start", "End", "Freq", "Strain", "ID")
write.table(b.f, file=paste(parent.o, "11_strains_repet_TE_frequency_no_redundant.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

#Read bed file from TIDAL info
x <- read.table(file=paste(parent.o, "11_strains_tidal_concatenated_indexed_sorted_bedtools.txt", sep=""), header=F, sep="\t")

#Count the number of strains in which each TE is present
a <- str_count(x[,4],",")
b <- x
b[,4] <- as.character(b[,4])
x[,4] <- as.character(x[,4])
f <- b.f <- NULL
f <- NULL
for (i in 1:length(a)) {
	if (a[i] > 11) {
		b[i,4] <- as.character(paste(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T))), collapse=","))
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} else {
	 	b[i,4] <- x[i,4]
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} 
}

b.f <- cbind(b[ ,1:3], f, b[ ,4:5])
colnames(b.f) <-  c("Chr_TE_ID", "Start", "End", "Freq", "Strain", "ID")
write.table(b.f, file=paste(parent.o, "11_strains_tidal_TE_frequency_no_redundant.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

#Read bed file from TEMP info
x <- read.table(file=paste(parent.o, "11_strains_temp_concatenated_indexed_sorted_bedtools.txt", sep=""), header=F, sep="\t")

#Count the number of strains in which each TE is present
a <- str_count(x[,4],",")
b <- x
b[,4] <- as.character(b[,4])
x[,4] <- as.character(x[,4])
f <- b.f <- NULL
f <- NULL
for (i in 1:length(a)) {
	if (a[i] > 11) {
		b[i,4] <- as.character(paste(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T))), collapse=","))
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} else {
	 	b[i,4] <- x[i,4]
		f <- c(f, length(unique(unlist(strsplit(as.character(b[i,4]), split=",", fixed=T)))))
	} 
}

b.f <- cbind(b[ ,1:3], f, b[ ,4:5])
colnames(b.f) <-  c("Chr_TE_ID", "Start", "End", "Freq", "Strain", "ID")
write.table(b.f, file=paste(parent.o, "11_strains_temp_TE_frequency_no_redundant.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")	

##Find TEs in each category
repet <- read.table(file=paste(parent.o, "11_strains_repet_TE_frequency_no_redundant.txt", sep=""), header=T, sep="\t")
tidal <- read.table(file=paste(parent.o, "11_strains_tidal_TE_frequency_no_redundant.txt", sep=""), header=T, sep="\t")
temp <- read.table(file=paste(parent.o, "11_strains_temp_TE_frequency_no_redundant.txt", sep=""), header=T, sep="\t")

#Create an ID for each TE from each method (REPET, TIDAL, TEMP)
r1 <- rep("REP", length(repet[,1]))
r2 <- seq(1, length(repet[,1]), 1)
r3 <- paste(r1, r2, sep="")

ti1 <- rep("TI", length(tidal[,1]))
ti2 <- seq(1, length(tidal[,1]), 1)
ti3 <- paste(ti1, ti2, sep="")

te1 <- rep("TE", length(temp[,1]))
te2 <- seq(1, length(temp[,1]), 1)
te3 <- paste(te1, te2, sep="")

#Assing the ID
repet.1 <- tidal.1 <- temp.1 <- all <- NULL
repet.1 <- cbind(repet, r3)
tidal.1 <- cbind(tidal, ti3)
temp.1 <- cbind(temp, te3)
			
colnames(repet.1)[7] <- "File"
colnames(tidal.1)[7] <- "File"
colnames(temp.1)[7] <- "File"

#Merge thes TEs intervals from REPET, TIDAL and TEMP
all <- NULL
all <- rbind(repet.1, tidal.1, temp.1)
write.table(all, file=paste(parent.o, "all_tes_non_redundant_indexed.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

