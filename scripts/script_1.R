#!/usr/bin/env Rscript

library("icesTAF")
library("DescTools") 
library("GenomicRanges")
library("tidyr")
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the REPET input files (i.e., "/home/")
#args[2]: path to the TIDAL input files (i.e., "/home/")
#args[3]: path to the TEMP input files (i.e., "/home/")
#args[4]: path to the output files (i.e., "/home/")
 
##Concatenate TEs info from strains and filter out those in heterochromatin region
#TEs detected using REPET
parent.r <- args[1]
list.r <- list.files(path=parent.r, pattern="*.bed")
parent.o <- args[4]

x <- y <- z <- NULL
for (i in 1:length(list.r)) {
	x <- read.table(file=paste(parent.r, list.r[i], sep=""), header=F, sep="\t")
	x[,1] <- gsub("2L", "chr2L", x[,1])
	x[,1] <- gsub("3L", "chr3L", x[,1])
	x[,1] <- gsub("2R", "chr2R", x[,1])
	x[,1] <- gsub("3R", "chr3R", x[,1])
	x[,1] <- gsub("X", "chrX", x[,1])
	name <- strsplit(list.r[i], split="_", fixed=T)[[1]][1]
	r1 <- rep(name, length(x[,1]))
	r2 <- seq(1, length(x[,1]), 1)
	r3 <- paste(r1, r2, sep="_ID")
	x <- cbind(x, rep(strsplit(list.r[i], split="_", fixed=T)[[1]][1], length(x[,1])), r3)
	y <- rbind(y, x)  
}
colnames(y) <- c("ChrREF","StartREF","EndREF","TEName","REFTE","Chr","Start","End","TEName.1","TESize","Strand","Family",
						"SuperFamily","Order","Class","SourceAnnot","TELenRatio","ClosestGenes","DistanceClosestGene","NumClosestGene","GeneDown",
						"GeneDownDist","GeneUp","GeneUpDist","Strain", "ID")

z <- y %>% unite(Chr_TE_ID, c(ChrREF, Family), sep = "&&", remove = TRUE)	
write.table(z, file=paste(parent.o, "11_strains_repet_concatenated.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")	

#TEs detected using TIDAL
parent.ti <- args[2]
list.ti <- list.files(path=parent.ti, pattern="*.txt")

x <- y <- z <- NULL
for (i in 1:length(list.ti)) {
	x <- read.csv(file=paste(parent.ti, list.ti[i], sep=""), header=T, sep="\t")
	name <- strsplit(list.ti[i], split="_", fixed=T)[[1]][1]
	r1 <- rep(name, length(x[,1]))
	r2 <- seq(1, length(x[,1]), 1)
	r3 <- paste(r1, r2, sep="_ID")
	x <- cbind(x, rep(strsplit(list.ti[i], split="_", fixed=T)[[1]][1], length(x[,1])), r3)
	y <- rbind(y, x)  
}

colnames(y)[20] <- "Strain"
colnames(y)[21] <-	"ID"

y[,2] <- as.character(y[,2])
eu <- NULL
for (i in 1:length(y[,1])) {
    if (y[i,2] == "chr2L" & y[i,3] > 530000 & y[i,4] < 18870000) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,2] == "chr2R" & y[i,3] > 5982495 & y[i,4] < 24972477) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,2] == "chr3L" & y[i,3] > 750000 & y[i,4] < 19026900) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,2] == "chr3R" & y[i,3] > 6754278 & y[i,4] < 31614278) {
	  	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,2] == "chrX" & y[i,3] > 1325967 & y[i,4] < 21338973) {
	   	eu  <- rbind(eu , y[i, ])
	}
}	

z <- eu %>% unite(Chr_TE_ID, c(Chr, TE), sep = "&&", remove = TRUE)	
write.table(z[-1], file=paste(parent.o, "11_strains_tidal_concatenated.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")	

#TEs detected using TEMP
parent.te <- args[3]
list.te <- list.files(path=parent.te, pattern="*.txt")
x <- y <- z <- NULL
for (i in 1:length(list.te)) {
	pr <- read.csv(file=paste(parent.te, list.te[i], sep=""), header=T, sep="\t")
	x <- pr[pr[,6] == "1p1", ]
	name <- strsplit(list.te[i], split="_", fixed=T)[[1]][1]
	r1 <- rep(name, length(x[,1]))
	r2 <- seq(1, length(x[,1]), 1)
	r3 <- paste(r1, r2, sep="_ID")
	x <- cbind(x, rep(strsplit(list.te[i], split="_", fixed=T)[[1]][1], length(x[,1])), r3)
	y <- rbind(y, x)  
}

colnames(y)[15] <- "Strain"
colnames(y)[16] <-	"ID"

y[,1] <- as.character(y[,1])
eu <- NULL
for (i in 1:length(y[,1])) {
    if (y[i,1] == "chr2L" & y[i,2] > 530000 & y[i,3] < 18870000) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,1] == "chr2R" & y[i,2] > 5982495 & y[i,3] < 24972477) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,1] == "chr3L" & y[i,2] > 750000 & y[i,3] < 19026900) {
	   	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,1] == "chr3R" & y[i,2] > 6754278 & y[i,3] < 31614278) {
	  	eu  <- rbind(eu , y[i, ])
	}
	if (y[i,1] == "chrX" & y[i,2] > 1325967 & y[i,3] < 21338973) {
	   	eu  <- rbind(eu , y[i, ])
	}
}	

z <- eu %>% unite(Chr_TE_ID, c(Chr, TransposonName), sep = "&&", remove = TRUE)	
write.table(z, file=paste(parent.o, "11_strains_temp_concatenated.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")	
