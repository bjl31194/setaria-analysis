library(ggplot2)

snps <- read.table(file = '~/setaria/results_files/SI_raw_snps_table.tsv', sep = '\t', header = TRUE)

snps <- subset(snps, select = -c(SI_1716_A.GT,SI_1716_G.GT,SI_1716_Y.GT,SI_575.GT))

###GENOTYPING###

#break read support into separate columns
snps$albino.ref <- gsub(",.*", "", snps$SI_1716_A.AD)
snps$albino.alt <- gsub(".*,", "", snps$SI_1716_A.AD)

snps$green.ref <- gsub(",.*", "", snps$SI_1716_G.AD)
snps$green.alt <- gsub(".*,", "", snps$SI_1716_G.AD)

snps$yellow.ref <- gsub(",.*", "", snps$SI_1716_Y.AD)
snps$yellow.alt <- gsub(".*,", "", snps$SI_1716_Y.AD)

snps$SI_575.ref <- gsub(",.*", "", snps$SI_575.AD)
snps$SI_575.alt <- gsub(".*,", "", snps$SI_575.AD)

#genotype snps
genotype_snps <- function(x,y) {
  if(as.integer(x) == 0 && as.integer(y) == 0) {
    return('-')
  }
  if(as.integer(x) == 0 && as.integer(y) > 0) {
    return('B')
  }
  else if(as.integer(x) / as.integer(y) >= 10) {
    return('A')
  }
  else if(as.integer(y) / as.integer(x) >= 10) {
    return('B')
  }
  else if(as.integer(x) / as.integer(y) < 10 && as.integer(x) / as.integer(y) > 4) {
    return('D')
  }
  else if(as.integer(y) / as.integer(x) < 10 && as.integer(y) / as.integer(x) > 4) {
    return('C')
  }
  else {
    return('H')
  }
}
for(j in 1:nrow(snps)) {
  snps[j,"albino.GT"] <- genotype_snps(snps[j,"albino.ref"],snps[j,"albino.alt"])
}

for(j in 1:nrow(snps)) {
  snps[j,"green.GT"] <- genotype_snps(snps[j,"green.ref"],snps[j,"green.alt"])
}

for(j in 1:nrow(snps)) {
  snps[j,"yellow.GT"] <- genotype_snps(snps[j,"yellow.ref"],snps[j,"yellow.alt"])
}

for(j in 1:nrow(snps)) {
  snps[j,"SI_575.GT"] <- genotype_snps(snps[j,"SI_575.ref"],snps[j,"SI_575.alt"])
}

#write.csv(snps, "~/setaria/snps_genotyped_readfilter.csv", row.names = FALSE)

###PLOTTING & READ DEPTH FILTERING### 

#splits AD fields by comma and puts sum in new column
if(is.character(snps$SI_1716_A.AD)) {
  snps$albino.total <- sapply(strsplit(snps$SI_1716_A.AD,","), function(x) sum(as.integer(x)))
}

if(is.character(snps$SI_1716_G.AD)) {
  snps$green.total <- sapply(strsplit(snps$SI_1716_G.AD,","), function(x) sum(as.integer(x)))
}

if(is.character(snps$SI_1716_Y.AD)) {
  snps$yellow.total <- sapply(strsplit(snps$SI_1716_Y.AD,","), function(x) sum(as.integer(x)))
}

if(is.character(snps$SI_575.AD)) {
  snps$SI_575.total <- sapply(strsplit(snps$SI_575.AD,","), function(x) sum(as.integer(x)))
}

#filters out snps with low read support in green/575 and extremely high support in all
snps <- subset(snps, snps$green.total > 2 & snps$SI_575.total > 2)
snps <- subset(snps, snps$albino.total < 30 & snps$green.total < 30 & snps$yellow.total < 30 & snps$SI_575.total < 30)

#divides alt reads by total per site
snps$albino.alt <- as.numeric(snps$albino.alt)
snps$albino.total <- as.numeric(snps$albino.total)
snps$green.alt <- as.numeric(snps$green.alt)
snps$green.total <- as.numeric(snps$green.total)
snps$yellow.alt <- as.numeric(snps$yellow.alt)
snps$yellow.total <- as.numeric(snps$yellow.total)
snps$SI_575.alt <- as.numeric(snps$SI_575.alt)
snps$SI_575.total <- as.numeric(snps$SI_575.total)

snps$green.AF <- snps$green.alt / snps$green.total
snps$albino.AF <- snps$albino.alt / snps$albino.total
snps$yellow.AF <- snps$yellow.alt / snps$yellow.total
snps$SI_575.AF <- snps$SI_575.alt / snps$SI_575.total


#make new column of %alt(A)-%alt(G)
snps$Aminus575.AF <- snps$albino.AF - snps$SI_575.AF

#plots A-SI_575 with mean line
ggplot(snps, aes(x = POS, y = Aminus575.AF)) +
  geom_point(size = 0.001) +
  labs(title = "Difference in ALT Allele Frequencies Between Albino and SI_575 Pools", x = "Genomic Position", y = "Difference in Allele Frequency (Albino - SI_575)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(vars(CHROM)) +
  stat_summary_bin(fun.y = 'mean', bins=500, color='orange', size=0.5, geom='line')

###FILTERING###
snps_filtered <- snps

#filter out non bialleleic snps
snps_filtered$bialleleic <- sapply(strsplit(snps_filtered$ALT,","), function(x) length(unique(x)) < 2)
snps_filtered <- subset(snps_filtered, snps_filtered$bialleleic == 'TRUE')

#filter snps that have same genotype in A as G or 575
snps_filtered <- subset(snps_filtered, (as.character(snps_filtered$albino.GT) != as.character(snps_filtered$green.GT)) & 
                          (as.character(snps_filtered$albino.GT) != as.character(snps_filtered$SI_575.GT)))

#find snps that are B/C/- in albino and A/D/H in green and 575
snps_filtered <- subset(snps_filtered, as.character(snps_filtered$albino.GT) == 'B' | as.character(snps_filtered$albino.GT) == 'C' | as.character(snps_filtered$albino.GT) == '-')

snps_filtered <- subset(snps_filtered, (as.character(snps_filtered$green.GT) == 'A') | (as.character(snps_filtered$green.GT) == 'D') | (as.character(snps_filtered$green.GT) == 'H'))
                                        
snps_filtered <- subset(snps_filtered, (as.character(snps_filtered$SI_575.GT) == 'A') | (as.character(snps_filtered$SI_575.GT) == 'D') | (as.character(snps_filtered$SI_575.GT) == 'H'))

#filters out snps that are the same genotype in albinos as yellows
snps_filtered <- subset(snps_filtered, (as.character(snps_filtered$albino.GT) != as.character(snps_filtered$yellow.GT)))

###ALTERNATE FILTERING STRATEGY - MORE STRINGENT

#only keeps snps that are B/- in albinos and A in 575
#snps_confident_filtered <- subset(snps_filtered, as.character(snps_filtered$albino.GT) == 'B' | as.character(snps_filtered$albino.GT) == '-')
#snps_confident_filtered <- subset(snps_filtered, (as.character(snps_filtered$SI_575.GT) == 'A')) 


#write filtered variants to new file
write.csv(snps_filtered, "~/setaria/results_files/snps_filtered_genotyped_full.csv", row.names = FALSE)
write.csv(snps_confident_filtered, "~/setaria/results_files/snps_BdelAlbino_A575.csv", row.names = FALSE)


###FINDING PREDICTED REFERENCE SEQUENCE ERRORS

#returns snps found at any frequency in the A/G/Y pool and 575
snps_genotyped_readfilter <- read.csv("~/setaria/results_files/snps_genotyped_readfilter.csv")

snps_refseq_errors <- subset(snps_genotyped_readfilter,
    (as.character(snps_genotyped_readfilter$albino.GT)) != "A" |
    (as.character(snps_genotyped_readfilter$green.GT)) != "A" |
    (as.character(snps_genotyped_readfilter$yellow.GT)) != "A")

snps_refseq_errors <- subset(snps_genotyped_readfilter,
    (as.character(snps_genotyped_readfilter$SI_575.GT)) != "A")

#filter out non bialleleic snps
snps_refseq_errors$bialleleic <- sapply(strsplit(snps_refseq_errors$ALT,","), function(x) length(unique(x)) < 2)
snps_refseq_errors <- subset(snps_refseq_errors, snps_refseq_errors$bialleleic == 'TRUE')

#write filtered variants to new file
write.csv(snps_refseq_errors, "~/setaria/results_files/snps_refseq_errors.csv", row.names = FALSE)

#write non biallelic snps to file
snps_non_biallelic <- subset(snps_refseq_errors, snps_refseq_errors$bialleleic == FALSE)
write.csv(snps_non_biallelic, "~/setaria/results_files/snps_non_biallelic.csv")
