#!/usr/bin/env Rscript

library("icesTAF")
library("DescTools") 
library("GenomicRanges")
library("tidyr")
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#args[1] path to the input and the output files (i.e., "/home/")
#The file Reference_tes.txt should be placed also in the same folder as the input/output files

##Get those TEs found with one, two and three methods (REPET, TIDAL and TEMP)
#Read bed file with all TEs detected using the three methods (REPET, TIDAL and TEMP)
parent.o <- args[1]
y <- read.table(file=paste(parent.o, "all_tes_non_redundant_indexed_sorted_bedtools.txt", sep=""), header=F, sep="\t")
all <- read.table(file=paste(parent.o, "all_tes_non_redundant_indexed.txt", sep=""), header=T, sep="\t")

#Separte those TE intervals with only one ID from those with multiple IDs
y1 <- y[grepl(",", y[ ,4]), ] #TE intervals with more than one ID
y2 <- y[!(grepl(",", y[ ,4])), ] #TE intervals with one ID

#Get the original info (Chr, TE, Start, End and Frequency) for each TE interval and separate those TE intervals with only one ID in REPET, TIDAL and TEMP
y3 <- all[as.character(all[,7]) %in% as.character(y2[,4]), ]
re <- y3[grepl("REP", y3[,7]), ]
ti <- y3[grepl("TI", y3[,7]), ]
te <- y3[grepl("TE", y3[,7]), ]

#Get the maximun number of IDs for each TE interval from the vector with multiple IDs
a1 <- NULL
for (i in 1:length(y1[,1])) {
	a1 <- c(a1, length(unlist(strsplit(as.character(y1[i,4]), split=","))))
}

#Get the original info (chr, TE, Start, End and Frequency) for each TE interval from IDs 
all.1 <- as.matrix(all)
b <- matrix( ,ncol=max(a1)*7, nrow=length(y1[,1])) 
for (i in 1:length(y1[,1])) {
	a <- length(unlist(strsplit(as.character(y1[i,4]), split=",")))
	pi <- sprintf("i is = %i", i)
	print(pi)
	for (j in 1:a) {
		if (j == 1) {
			b[i, 1:7] <- all.1[all.1[,7] %in% strsplit(as.character(y1[i,4]), split=",")[[1]][j], 1:7]
		} else {
			j1 <- (((j-1)*7)+1)
			j2 <- (((j-1)*7)+7)
			b[i, j1:j2] <- all.1[all.1[,7] %in% strsplit(as.character(y1[i,4]), split=",")[[1]][j], 1:7]
		}
	}	
}

#Assing to each vector their corresponding method of detection (for those TE intervals with only one ID)
re.1 <- re
re.1$File <- as.character(re.1$File)
re.1$File[grepl("REP", re.1$File)] <- "REP"

ti.1 <- ti
ti.1$File <- as.character(ti.1$File)
ti.1$File[grepl("TI", ti.1$File)] <- "TI"

te.1 <- te
te.1$File <-as.character(te.1$File)
te.1$File[grepl("TE", te.1$File)] <- "TE"

#Generate a vector with the info where the name of the origin of the file (REP, TI, TE) is	
b.1 <- b
co <- NULL
for (i in 1:(max(a1))) { co <- c(co, 7*i) }

#Generate a vector with the info where the chr and TE info is 
co.1 <- NULL
for (i in 1:(max(a1)-1)) { co.1 <- c(co.1, (7*i)+1) }

co.2 <- setdiff(seq(1, 7*max(a1), 1), co.1)
for (i in 1:length(co)) {
	b.1[,co[i]]	<- as.character(b.1[,co[i]])
}

for (i in 1:length(co)) {
	b.1[,co[i]][ grepl("REP",b.1[,co[i]]) ] <- "REP"
	b.1[,co[i]][ grepl("TI",b.1[,co[i]]) ] <- "TI"
	b.1[,co[i]][ grepl("TE",b.1[,co[i]]) ] <- "TE"
}
b.2 <- b.1[ ,co.2]

#Get the info for each comparison (TEs found by two and three methods)
rep.ti <- b.2[!(grepl("TE", b.2[ ,7])) & !(grepl("TE", b.2[ ,13])) & !(grepl("TE", b.2[ ,19])) & !(grepl("TE", b.2[ ,25])) & !(grepl("TE", b.2[ ,31])) & !(grepl("TE", b.2[ ,37])), ] 
rep.te <- b.2[!(grepl("TI", b.2[ ,7])) & !(grepl("TI", b.2[ ,13])) & !(grepl("TI", b.2[ ,19])) & !(grepl("TI", b.2[ ,25])) & !(grepl("TI", b.2[ ,31])) & !(grepl("TI", b.2[ ,37])), ] 
ti.te <- b.2[!(grepl("REP", b.2[ ,7])) & !(grepl("REP", b.2[ ,13])) & !(grepl("REP", b.2[ ,19])) & !(grepl("REP", b.2[ ,25])) & !(grepl("REP", b.2[ ,31])) & !(grepl("REP", b.2[ ,37])), ] 

rep.ti.te <- NULL
for (i in 1:length(b.2[,1])) {
	if ( ("REP" %in% b.2[i, ]) & ("TI" %in% b.2[i, ]) & ("TE" %in% b.2[i, ]) ) {
		rep.ti.te <- rbind(rep.ti.te, b.2[i, ])
	}
}

#Separate the info (in the first column) about the Chr and the TE family
#REPET
all.re.1 <- NULL
all.re.1 <- cbind(do.call(rbind, str_split(re.1$Chr_TE_ID, '&&')), re.1[, 2:7])

#TIDAL
all.ti.1 <- NULL
all.ti.1 <- cbind(do.call(rbind, str_split(ti.1$Chr_TE_ID, '&&')), ti.1[, 2:7])

#TEMP
all.te.1 <- NULL
all.te.1 <- cbind(do.call(rbind, str_split(te.1$Chr_TE_ID, '&&')), te.1[, 2:7])

#REPET_TIDAL
all.rep.ti.1 <- NULL
all.rep.ti.1 <- cbind(do.call(rbind, str_split(rep.ti[,1], '&&')), rep.ti[, 2:37])

#REPET_TEMP
all.rep.te.1 <- NULL
all.rep.te.1 <- cbind(do.call(rbind, str_split(rep.te[,1], '&&')), rep.te[, 2:37])

#TIDAL-TEMP
all.ti.te.1 <- NULL
all.ti.te.1 <- cbind(do.call(rbind, str_split(ti.te[,1], '&&')), ti.te[, 2:37])

#REPET-TIDAL-TEMP
all.rep.ti.te.1 <- NULL
all.rep.ti.te.1 <- cbind(do.call(rbind, str_split(rep.ti.te[,1], '&&')), rep.ti.te[, 2:37])

colnames(all.re.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File")
colnames(all.ti.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File")  
colnames(all.te.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File")  
colnames(all.rep.ti.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File", rep(c("Start","End","Freq","Strain","ID","File"), 5))  
colnames(all.rep.te.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File", rep(c("Start","End","Freq","Strain","ID","File"), 5))
colnames(all.ti.te.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File", rep(c("Start","End","Freq","Strain","ID","File"), 5)) 
colnames(all.rep.ti.te.1) <- c("Chr","TE_ID","Start","End","Freq","Strain","ID","File", rep(c("Start","End","Freq","Strain","ID","File"), 5))  

write.table(all.re.1, file=paste(parent.o, "repet_only_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.ti.1, file=paste(parent.o, "tidal_only_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.te.1, file=paste(parent.o, "temp_only_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.ti.1, file=paste(parent.o, "repet_tidal_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.te.1, file=paste(parent.o, "repet_temp_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.ti.te.1, file=paste(parent.o, "tidal_temp_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.ti.te.1, file=paste(parent.o, "repet_tidal_temp_euch_no_redundant.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

##Count the number of TEs in each category (TEs found by one, two and three methods)
dim(repet)[1]
dim(tidal)[1]
dim(temp)[1]
dim(re)[1]
dim(ti)[1]
dim(te)[1]

#REPET
rep.ti.te.s <- rowSums(rep.ti.te == 'REP', na.rm=T)
ct <- NULL
for (i in min(rep.ti.te.s):max(rep.ti.te.s)) {
	ct <- c(ct, rep.ti.te.s[rep.ti.te.s == i])
}
tot.rep.ti.te.r <- sum(ct)

rep.ti.s <- rowSums(rep.ti == 'REP', na.rm=T)
ct <- NULL
for (i in min(rep.ti.s):max(rep.ti.s)) {
	ct <- c(ct, rep.ti.s[rep.ti.s == i])
}
tot.rep.ti.r <- sum(ct)

rep.te.s <- rowSums(rep.te == 'REP', na.rm=T)
ct <- NULL
for (i in min(rep.te.s):max(rep.te.s)) {
	ct <- c(ct, rep.te.s[rep.te.s == i])
}
tot.rep.te.r <- sum(ct)
all.rep <- tot.rep.ti.te.r + tot.rep.ti.r + tot.rep.te.r + dim(re)[1]

#TIDAL
rep.ti.te.s <- rowSums(rep.ti.te == 'TI', na.rm=T)
ct <- NULL
for (i in min(rep.ti.te.s):max(rep.ti.te.s)) {
	ct <- c(ct, rep.ti.te.s[rep.ti.te.s == i])
}
tot.rep.ti.te.ti <- sum(ct)

rep.ti.s <- rowSums(rep.ti == 'TI', na.rm=T)
ct <- NULL
for (i in min(rep.ti.s):max(rep.ti.s)) {
	ct <- c(ct, rep.ti.s[rep.ti.s == i])
}
tot.rep.ti.ti <- sum(ct)

ti.te.s <- rowSums(ti.te == 'TI', na.rm=T)
ct <- NULL
for (i in min(ti.te.s):max(ti.te.s)) {
	ct <- c(ct, ti.te.s[ti.te.s == i])
}
tot.ti.te.ti <- sum(ct)
all.ti <- tot.rep.ti.te.ti + tot.rep.ti.ti + tot.ti.te.ti + dim(ti)[1]

#TEMP
rep.ti.te.s <- rowSums(rep.ti.te == 'TE', na.rm=T)
ct <- NULL
for (i in min(rep.ti.te.s):max(rep.ti.te.s)) {
	ct <- c(ct, rep.ti.te.s[rep.ti.te.s == i])
}
tot.rep.ti.te.te <- sum(ct)

rep.te.s <- rowSums(rep.te == 'TE', na.rm=T)
ct <- NULL
for (i in min(rep.te.s):max(rep.te.s)) {
	ct <- c(ct, rep.te.s[rep.te.s == i])
}
tot.rep.te.te <- sum(ct)

ti.te.s <- rowSums(ti.te == 'TE', na.rm=T)
ct <- NULL
for (i in min(ti.te.s):max(ti.te.s)) {
	ct <- c(ct, ti.te.s[ti.te.s == i])
}
tot.ti.te.te <- sum(ct)
all.te <- tot.rep.ti.te.te + tot.rep.te.te + tot.ti.te.te + dim(te)[1]

table.final <- cbind(
c("repet", "only repet", "repet tidal -based on repet-", "repet temp -based on repet-", "repet tidal temp -based on repet-",
"tidal", "only tidal", "repet tidal -based on tidal-", "tidal temp -based on tidal-", "repet tidal temp -based on tidal-",
"temp", "only temp", "repet temp -based on temp-", "tidal temp -based on temp-", "repet tidal temp -based on temp-"), 
c(all.rep, dim(re)[1], tot.rep.ti.r,  tot.rep.te.r, tot.rep.ti.te.r,  
  all.ti,  dim(ti)[1], tot.rep.ti.ti, tot.ti.te.ti, tot.rep.ti.te.ti, 
  all.te,  dim(te)[1], tot.rep.te.te, tot.ti.te.te, tot.rep.ti.te.te)
)
write.table(table.final, file=paste(parent.o, "repet_tidal_temp_TE_counts_for_each_category.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

##Count the number of polymorphic TEs in each category (TEs found by one, two and three methods)

all.re.1.p <- all.re.1[as.numeric(all.re.1[,5]) != 1, ]
all.ti.1.p <- all.ti.1[as.numeric(all.ti.1[,5]) != 1, ]
all.te.1.p <- all.te.1[as.numeric(all.te.1[,5]) != 1, ]

all.rep.ti.2 <- as.data.frame(all.rep.ti.1)
all.rep.te.2 <- as.data.frame(all.rep.te.1)
all.ti.te.2 <- as.data.frame(all.ti.te.1)
all.rep.ti.te.2 <- as.data.frame(all.rep.ti.te.1)

all.rep.ti.1.p <- NULL
for (i in 1:length(all.rep.ti.2[,1])) {
	if ( (as.numeric(all.rep.ti.2[i,5]) != 1 | is.na(all.rep.ti.2[i,5])) & 
					(as.numeric(all.rep.ti.2[i,11]) != 1 | is.na(all.rep.ti.2[i,11])) &  
					(as.numeric(all.rep.ti.2[i,17]) != 1 | is.na(all.rep.ti.2[i,17])) &  
					(as.numeric(all.rep.ti.2[i,23]) != 1 | is.na(all.rep.ti.2[i,23])) & 
					(as.numeric(all.rep.ti.2[i,29]) != 1 | is.na(all.rep.ti.2[i,29])) &  
					(as.numeric(all.rep.ti.2[i,35]) != 1 | is.na(all.rep.ti.2[i,35]))  ) {
		all.rep.ti.1.p <- rbind(all.rep.ti.1.p, all.rep.ti.2[i,  ])			
	}	
}				
					
all.rep.te.1.p <- NULL
for (i in 1:length(all.rep.te.2[,1])) {
	if ( (as.numeric(all.rep.te.2[i,5]) != 1 | is.na(all.rep.te.2[i,5])) & 
					(as.numeric(all.rep.te.2[i,11]) != 1 | is.na(all.rep.te.2[i,11])) &  
					(as.numeric(all.rep.te.2[i,17]) != 1 | is.na(all.rep.te.2[i,17])) &  
					(as.numeric(all.rep.te.2[i,23]) != 1 | is.na(all.rep.te.2[i,23])) & 
					(as.numeric(all.rep.te.2[i,29]) != 1 | is.na(all.rep.te.2[i,29])) &  
					(as.numeric(all.rep.te.2[i,35]) != 1 | is.na(all.rep.te.2[i,35]))  ) {
		all.rep.te.1.p <- rbind(all.rep.te.1.p, all.rep.te.2[i,  ])			
	}	
}

all.ti.te.1.p <- NULL
for (i in 1:length(all.ti.te.2[,1])) {
	if ( (as.numeric(all.ti.te.2[i,5]) != 1 | is.na(all.ti.te.2[i,5])) & 
					(as.numeric(all.ti.te.2[i,11]) != 1 | is.na(all.ti.te.2[i,11])) &  
					(as.numeric(all.ti.te.2[i,17]) != 1 | is.na(all.ti.te.2[i,17])) &  
					(as.numeric(all.ti.te.2[i,23]) != 1 | is.na(all.ti.te.2[i,23])) & 
					(as.numeric(all.ti.te.2[i,29]) != 1 | is.na(all.ti.te.2[i,29])) &  
					(as.numeric(all.ti.te.2[i,35]) != 1 | is.na(all.ti.te.2[i,35]))  ) {
		all.ti.te.1.p <- rbind(all.ti.te.1.p, all.ti.te.2[i,  ])			
	}	
}	

all.rep.ti.te.1.p <- NULL
for (i in 1:length(all.rep.ti.te.2[,1])) {
	if ( (as.numeric(all.rep.ti.te.2[i,5]) != 1 | is.na(all.rep.ti.te.2[i,5])) & 
					(as.numeric(all.rep.ti.te.2[i,11]) != 1 | is.na(all.rep.ti.te.2[i,11])) &  
					(as.numeric(all.rep.ti.te.2[i,17]) != 1 | is.na(all.rep.ti.te.2[i,17])) &  
					(as.numeric(all.rep.ti.te.2[i,23]) != 1 | is.na(all.rep.ti.te.2[i,23])) & 
					(as.numeric(all.rep.ti.te.2[i,29]) != 1 | is.na(all.rep.ti.te.2[i,29])) &  
					(as.numeric(all.rep.ti.te.2[i,35]) != 1 | is.na(all.rep.ti.te.2[i,35]))  ) {
		all.rep.ti.te.1.p <- rbind(all.rep.ti.te.1.p, all.rep.ti.te.2[i,  ])			
	}	
}

write.table(all.re.1.p, file=paste(parent.o, "repet_only_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.ti.1.p, file=paste(parent.o,"tidal_only_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.te.1.p, file=paste(parent.o,"temp_only_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.ti.1.p, file=paste(parent.o,"repet_tidal_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.te.1.p, file=paste(parent.o,"repet_temp_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.ti.te.1.p, file=paste(parent.o,"tidal_temp_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
write.table(all.rep.ti.te.1.p, file=paste(parent.o,"repet_tidal_temp_euch_no_redundant_pol.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

##Count the number of polymorphic TEs in each category (TEs in more than one strain)
dim(repet[repet[,4] != 1, ])[1]
dim(tidal[tidal[,4] != 1, ])[1]
dim(temp[temp[,4] != 1, ])[1]
dim(all.re.1.p)[1]
dim(all.ti.1.p)[1]
dim(all.te.1.p)[1]

#REPET
all.rep.ti.te.1.p.s <- rowSums(all.rep.ti.te.1.p == 'REP', na.rm=T)
ct <- NULL
for (i in min(all.rep.ti.te.1.p.s):max(all.rep.ti.te.1.p.s)) {
	ct <- c(ct, all.rep.ti.te.1.p.s[all.rep.ti.te.1.p.s == i])
}
tot.all.rep.ti.te.1.p.r <- sum(ct)

all.rep.ti.1.p.s <- rowSums(all.rep.ti.1.p == 'REP', na.rm=T)
ct <- NULL
for (i in min(all.rep.ti.1.p.s):max(all.rep.ti.1.p.s)) {
	ct <- c(ct, all.rep.ti.1.p.s[all.rep.ti.1.p.s == i])
}
tot.all.rep.ti.1.p.r <- sum(ct)

all.rep.te.1.p.s <- rowSums(all.rep.te.1.p == 'REP', na.rm=T)
ct <- NULL
for (i in min(all.rep.te.1.p.s):max(all.rep.te.1.p.s)) {
	ct <- c(ct, all.rep.te.1.p.s[all.rep.te.1.p.s == i])
}
tot.all.rep.te.1.p.r <- sum(ct)
all.rep.p <- tot.all.rep.ti.te.1.p.r + tot.all.rep.ti.1.p.r + tot.all.rep.te.1.p.r + dim(re[re[,4] != 1, ])[1]

#TIDAL
all.rep.ti.te.1.p.s <- rowSums(all.rep.ti.te.1.p == 'TI', na.rm=T)
ct <- NULL
for (i in min(all.rep.ti.te.1.p.s):max(all.rep.ti.te.1.p.s)) {
	ct <- c(ct, all.rep.ti.te.1.p.s[all.rep.ti.te.1.p.s == i])
}
tot.all.rep.ti.te.1.p.ti <- sum(ct)

all.rep.ti.1.p.s <- rowSums(all.rep.ti.1.p == 'TI', na.rm=T)
ct <- NULL
for (i in min(all.rep.ti.1.p.s):max(all.rep.ti.1.p.s)) {
	ct <- c(ct, all.rep.ti.1.p.s[all.rep.ti.1.p.s == i])
}
tot.all.rep.ti.1.p.ti <- sum(ct)

all.ti.te.1.p.s <- rowSums(all.ti.te.1.p == 'TI', na.rm=T)
ct <- NULL
for (i in min(all.ti.te.1.p.s):max(all.ti.te.1.p.s)) {
	ct <- c(ct, all.ti.te.1.p.s[all.ti.te.1.p.s == i])
}
tot.all.ti.te.1.p.ti <- sum(ct)
all.ti.p <- tot.all.rep.ti.te.1.p.ti + tot.all.rep.ti.1.p.ti + tot.all.ti.te.1.p.ti + dim(ti[ti[,4] != 1, ])[1]

#TEMP
all.rep.ti.te.1.p.s <- rowSums(all.rep.ti.te.1.p == 'TE', na.rm=T)
ct <- NULL
for (i in min(all.rep.ti.te.1.p.s):max(all.rep.ti.te.1.p.s)) {
	ct <- c(ct, all.rep.ti.te.1.p.s[all.rep.ti.te.1.p.s == i])
}
tot.all.rep.ti.te.1.p.te <- sum(ct)

all.rep.te.1.p.s <- rowSums(all.rep.te.1.p == 'TE', na.rm=T)
ct <- NULL
for (i in min(all.rep.te.1.p.s):max(all.rep.te.1.p.s)) {
	ct <- c(ct, all.rep.te.1.p.s[all.rep.te.1.p.s == i])
}
tot.all.rep.te.1.p.te <- sum(ct)

all.ti.te.1.p.s <- rowSums(all.ti.te.1.p == 'TE', na.rm=T)
ct <- NULL
for (i in min(all.ti.te.1.p.s):max(all.ti.te.1.p.s)) {
	ct <- c(ct, all.ti.te.1.p.s[all.ti.te.1.p.s == i])
}
tot.all.ti.te.1.p.te <- sum(ct)
all.te.p <- tot.all.rep.ti.te.1.p.te + tot.all.rep.te.1.p.te + tot.all.ti.te.1.p.te + dim(te[te[,4] != 1, ])[1]

table.final.p <- NULL
table.final.p <- cbind(
c("repet", "only repet", "repet tidal temp -based on repet-", "repet tidal -based on repet-", "repet temp -based on repet-",
"tidal", "only tidal", "repet tidal temp -based on tidal-", "repet tidal -based on tidal-", "tidal temp -based on tidal-",
"temp", "only temp", "repet tidal temp -based on temp-", "repet temp -based on temp-", "tidal temp -based on temp-"), 
c(all.rep.p, dim(re[re[,4] != 1, ])[1], tot.all.rep.ti.te.1.p.r, tot.all.rep.ti.1.p.r, tot.all.rep.te.1.p.r, 
  all.ti.p, dim(ti[ti[,4] != 1, ])[1], tot.all.rep.ti.te.1.p.ti, tot.all.rep.ti.1.p.ti, tot.all.ti.te.1.p.ti, 
  all.te.p, dim(te[te[,4] != 1, ])[1], tot.all.rep.ti.te.1.p.te, tot.all.rep.te.1.p.te, tot.all.ti.te.1.p.te)
)

write.table(table.final.p, file=paste(parent.o, "repet_tidal_temp_TE_counts_for_each_category_pol.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

##Include TEs TYPE (those found in TIDAL and TEMP are named REF and the rest are named REPET
tes <- read.table(file=paste(parent.o, "Reference_tes.txt", sep=""), header=F, sep="\t")

all <- NULL
list.f <- list.files(path=parent.o, pattern="*_no_redundant.txt")

for (i in 1:length(list.f)) {
		all <- NULL
		x <- read.table(file=paste(parent.o ,list.f[i], sep=""), header=T, sep="\t")
		name <- strsplit(list.f[i], split=".txt")
		ref <- x[x[,2] %in% tes[,1], ]
		rep <- x[!(x[,2] %in% tes[,1]), ]
		ref.1 <- cbind(ref, rep("REF", length(ref[,1])))
		rep.1 <- cbind(rep, rep("REPET", length(rep[,1])))
		colnames(ref.1)[length(x[1,])+1] <- "Type"
		colnames(rep.1)[length(x[1,])+1] <- "Type"
		all <- rbind(ref.1, rep.1)
		write.table(all, file=paste(parent.o, name, "_info.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}

##Include TEs TYPE for polymorphic TEs (those found in TIDAL and TEMP are named REF and the rest are named REPET
list.f <- list.files(path=parent.o, pattern="*_no_redundant_pol.txt")
tes <- read.table(file=paste(parent.o, "Reference_tes.txt", sep=""), header=F, sep="\t")
all <- NULL

for (i in 1:length(list.f)) {
		all <- NULL
		x <- read.table(file=paste(parent.o ,list.f[i], sep=""), header=T, sep="\t")
		name <- strsplit(list.f[i], split=".txt")
		ref <- x[x[,2] %in% tes[,1], ]
		rep <- x[!(x[,2] %in% tes[,1]), ]
		ref.1 <- cbind(ref, rep("REF", length(ref[,1])))
		rep.1 <- cbind(rep, rep("REPET", length(rep[,1])))
		colnames(ref.1)[length(x[1,])+1] <- "Type"
		colnames(rep.1)[length(x[1,])+1] <- "Type"
		all <- rbind(ref.1, rep.1)
		write.table(all, file=paste(parent.o, name, "_info.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}


##Include information about TE ID (in REPET) and the length of TEs
#Only REPET
x <- read.table(file=paste(parent.o, "11_strains_repet_concatenated.txt", sep=""), header=T, sep="\t")
y <- read.table(file=paste(parent.o, "repet_only_euch_no_redundant_info.txt", sep=""), header=T, sep="\t")
y1 <- y[grepl(",", y[ ,7]), ] #TE intervals with more than one ID
y2 <- y[!(grepl(",", y[ ,7])), ] #TE intervals with one ID

z <- NULL
z <- x[x[ ,25] %in% y2[, 7], c(4,10)]
z1 <- x[x[ ,25] %in% y2[, 7], c(4,10, 25)]
z2 <- y[y[ ,7] %in% z1[ ,3], ]

z1o <- z1[order(z1$ID), ]
z2o <- z2[order(z2$ID), ]
y3 <- NULL
y3 <- cbind(z2o[,1], z1o[,1:2],z2o[,2:9])

colnames(y3)[2] <- "TE_name"
colnames(y3)[3] <- "TE_length"
colnames(y3)[1] <- "Chr" 

all.id <- all.ln <- NULL
for (j in 1:dim(y1)[1]) {
	pj <- sprintf("j is = %i", j)
	print(pj)
	a <- length(unlist(strsplit(as.character(y1[j,7]), split=",")))
	id <- ln <- NULL
	for (k in 1:a) {
		id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y1[j, 7]), split=",")[[1]][k]), 4]))
		ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y1[j, 7]), split=",")[[1]][k]), 10]))
	}
	all.id <- rbind(all.id, paste(id, collapse=","))	
	all.ln <- rbind(all.ln, paste(ln, collapse=";"))
}
all.id.ln <- cbind(all.id, all.ln)
y4 <- NULL
y4 <- cbind(y1[,1], all.id.ln, y1[,2:9])
colnames(y4)[2] <- "TE_name"
colnames(y4)[3] <- "TE_length"
colnames(y4)[1] <- "Chr"
y5 <- NULL
y5 <- rbind(y3,y4)
write.table(y5, file=paste(parent.o, "repet_only_euch_no_redundant_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TIDAL		 
y <- read.table(file=paste(parent.o, "repet_tidal_euch_no_redundant_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
			}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_tidal_euch_no_redundant_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TEMP		 
y <- read.table(file=paste(parent.o, "repet_temp_euch_no_redundant_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
			}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}	
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_temp_euch_no_redundant_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TIDAL-TEMP
y <- read.table(file=paste(parent.o, "repet_tidal_temp_euch_no_redundant_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
			}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_tidal_temp_euch_no_redundant_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

##Include information about TE ID (in Repet) and the length of polymorphic TEs 
#Only REPET

x <- read.table(file=paste(parent.o, "11_strains_repet_concatenated.txt", sep=""), header=T, sep="\t")
y <- read.table(file=paste(parent.o, "repet_only_euch_no_redundant_pol_info.txt", sep=""), header=T, sep="\t")
y1 <- y[grepl(",", y[ ,7]), ] #TE intervals with more than one ID
y2 <- y[!(grepl(",", y[ ,7])), ] #TE intervals with one ID

z <- NULL
z <- x[x[, 25] %in% y2[, 7], c(4,10)]
y3 <- NULL
y3 <- cbind(y2[,1],z, y2[ ,2:9])
colnames(y3)[2] <- "TE_name"
colnames(y3)[3] <- "TE_length"
colnames(y3)[1] <- "Chr" 

all.id <- all.ln <- NULL
for (j in 1:dim(y1)[1]) {
	pj <- sprintf("j is = %i", j)
	print(pj)
	a <- length(unlist(strsplit(as.character(y1[j,7]), split=",")))
	id <- ln <- NULL
	for (k in 1:a) {
		id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y1[j, 7]), split=",")[[1]][k]), 4]))
		ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y1[j, 7]), split=",")[[1]][k]), 10]))
	}
	all.id <- rbind(all.id, paste(id, collapse=","))	
	all.ln <- rbind(all.ln, paste(ln, collapse=";"))
}
all.id.ln <- cbind(all.id, all.ln)
y4 <- NULL
y4 <- cbind(y1[,1], all.id.ln, y1[,2:9])
colnames(y4)[2] <- "TE_name"
colnames(y4)[3] <- "TE_length"
colnames(y4)[1] <- "Chr"
y5 <- NULL
y5 <- rbind(y3,y4)
write.table(y5, file=paste(parent.o, "repet_only_euch_no_redundant_pol_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TIDAL		 
y <- read.table(file=paste(parent.o, "repet_tidal_euch_no_redundant_pol_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
			}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_tidal_euch_no_redundant_pol_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TEMP	 
y <- read.table(file=paste(parent.o, "repet_temp_euch_no_redundant_pol_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
			}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}	
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_temp_euch_no_redundant_pol_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

#REPET-TIDAL-TEMP
y <- read.table(file=paste(parent.o, "repet_tidal_temp_euch_no_redundant_pol_info.txt", sep=""), header=T, sep="\t")
c <- rowSums(y == 'REP', na.rm=T) #how many REP are in each row
if (max(c) <= 1) {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
		id <- ln <- NULL
		for (j in 1:a) {
			id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 4]))
			ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][j]), 10]))
		}
		all.id <- rbind(all.id, paste(id, collapse=","))	
		all.ln <- rbind(all.ln, paste(ln, collapse=";"))
	}
	y6 <- cbind(y[,1], all.id, all.ln, y[, 2:39]) 	
	colnames(y6)[2] <- "TE_name"
	colnames(y6)[3] <- "TE_length"
	colnames(y6)[1] <- "Chr"	
} else {
	all.id <- all.ln <- NULL
	for (i in 1:length(c)) {
		if (c[i] == 1) {
			a <- length(unlist(strsplit(as.character(y[i,abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")))
			id <- ln <- NULL
			for (m in 1:a) {
				id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 4]))
				ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, abs(1 - as.numeric(which(y[i, ] == "REP")))]), split=",")[[1]][m]), 10]))
}
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		} else {
			b <- abs(1 - as.numeric(which(y[i, ] == "REP")))
			id <- ln <- NULL
			for (n in 1:length(b)) { 
				a <- length(unlist(strsplit(as.character(y[i, b[n]]), split=",")))
				for (o in 1:a) {
					id <- c(id, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 4]))
					ln <- c(ln, as.character(x[x[, 25] %in% as.character(strsplit(as.character(y[i, b[n]]), split=",")[[1]][o]), 10]))
				}
			}	
			all.id <- rbind(all.id, paste(id, collapse=","))	
			all.ln <- rbind(all.ln, paste(ln, collapse=";"))
		}		
	}
}
y7 <- NULL
y7 <- cbind(y[,1], all.id, all.ln, y[,2:39])
colnames(y7)[2] <- "TE_name"
colnames(y7)[3] <- "TE_length"
colnames(y7)[1] <- "Chr"
write.table(y7, file=paste(parent.o, "repet_tidal_temp_euch_no_redundant_pol_info_length_id.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

	