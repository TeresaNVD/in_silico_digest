
#Install library
#source("https://bioconductor.org/biocLite.R")
#biocLite("cleaver")

#Load the library
library("cleaver")

#Set-up missed cleavages
cleave("LAAGKVEDSD", enzym="asp-n", missedCleavages=0)
cleave("LAAGKVEDSD", enzym="asp-n", missedCleavages=1)
cleave("LAAGKVEDSD", enzym="asp-n", missedCleavages=c(0:2))

#Other information
cleavageRanges("LAAGKVEDSD", enzym="asp-n")
cleavageSites("LAAGKVEDSD", enzym="asp-n")

#Digest data from a collection of strings
#Load the library
library(Biostrings)
p <- AAStringSet(c(gaju="LAAGKVEDSD", pnm="AGEPKLDAGV"))
cleave(p, enzym="asp-n")


#Digest data from UniProtKB (This approach is recommended for a few hundred proteins)
#Load the library
library(UniProt.ws)
#fetch organism data
UniProt.ws <- UniProt.ws(taxId=9606)
#fetch all uniprotids for 9606
all_uniprotID <- UniProt.ws@taxIdUniprots
#fetch sequences for the uniprotids and name the sequence by uniprotid
s <- select(UniProt.ws, keys=all_uniprotID, columns=c("SEQUENCE"))
sequences <- setNames(s$SEQUENCE, s$UNIPROTKB)
sequences <- gsub(pattern="[[:space:]]", replacement="", x=sequences)
#perform cleavage
cleavage_zero <- cleave(sequences, enzym="asp-n", missedCleavages=0)
cleavage_zero_one <- cleave(sequences, enzym="asp-n", missedCleavages=0:1)

#Digest data from downloaded proteome
library(Biostrings)
p <- readAAStringSet("/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/uniprot-proteome%3AUP000005640+reviewed%3Ayes_cannonical_20-11-2017.fasta", "fasta")
pp <- data.frame(p)
sequences <- setNames(pp$p, row.names(pp))

#******ASP-N ANALYSIS FOR ZERO MISSED CLEAVAGE****
digest <- cleave(sequences, enzym="asp-n", missedCleavages=0)

#Statistics for CK2 sites
#all digestion peptides
a <- lapply(digest, function(x){lapply(x, function(y){nchar(y)})})
length(unlist(a))
a <- unlist(a)
hist(a)
hist(a[a<=50])

#specific length digested peptides
int_pept <- lapply(digest, function(x){lapply(x, function(y){y[which(nchar(y) >= 8 & nchar(y) <= 37)]})})
length(unlist(int_pept))
int_pept_ms <- unlist(int_pept)
int_pept_ms_length <- sapply(int_pept_ms, function(x){nchar(x)})
hist(int_pept_ms_length)

int_pept_ms_ck2 <- int_pept_ms[grep("[ST][DE][A-Z][DE]", int_pept_ms)]
length(int_pept_ms_ck2)

int_pept_ms_ck2_df <- data.frame(name = names(int_pept_ms_ck2), peptide = int_pept_ms_ck2)
dim(int_pept_ms_ck2_df)

int_pept_ms_ck2_df_not_kr <- int_pept_ms_ck2[grep("[KR]", int_pept_ms_ck2, invert = TRUE)]
int_pept_ms_ck2_df_not_kr <- data.frame(name = names(int_pept_ms_ck2_df_not_kr), peptide = int_pept_ms_ck2_df_not_kr)
dim(int_pept_ms_ck2_df_not_kr)

#******ASP-N ANALYSIS FOR ONE MISSED CLEAVAGES****
digest <- cleave(sequences, enzym="asp-n", missedCleavages=1)

#Statistics for CK2 sites
#all digestion peptides
a <- lapply(digest, function(x){lapply(x, function(y){nchar(y)})})
length(unlist(a))
a <- unlist(a)
hist(a)
hist(a[a<=50])

#specific length digested peptides
int_pept <- lapply(digest, function(x){lapply(x, function(y){y[which(nchar(y) >= 8 & nchar(y) <= 37)]})})
length(unlist(int_pept))
int_pept_df <- data.frame(name = names(unlist(int_pept)), peptide = unlist(int_pept))
rownames(int_pept_df) <- NULL
dim(int_pept_df)
write.csv(int_pept_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/asp-n_1missed_8-37pep.csv")

int_pept_ms <- unlist(int_pept)
int_pept_ms_length <- sapply(int_pept_ms, function(x){nchar(x)})
hist(int_pept_ms_length)

int_pept_ms_ck2 <- int_pept_ms[grep("[ST][DE][A-Z][DE]", int_pept_ms)]
length(int_pept_ms_ck2)
int_pept_ms_ck2_df <- data.frame(name = names(int_pept_ms_ck2), peptide = int_pept_ms_ck2)
rownames(int_pept_ms_ck2_df) <- NULL
dim(int_pept_ms_ck2_df)
write.csv(int_pept_ms_ck2_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/asp-n_1missed_8-37pep_ck2.csv")


int_pept_ms_ck2_df_not_kr <- int_pept_ms_ck2[grep("[KR]", int_pept_ms_ck2, invert = TRUE)]
int_pept_ms_ck2_df_not_kr <- data.frame(name = names(int_pept_ms_ck2_df_not_kr), peptide = int_pept_ms_ck2_df_not_kr)
rownames(int_pept_ms_ck2_df_not_kr) <- NULL
dim(int_pept_ms_ck2_df_not_kr)
write.csv(int_pept_ms_ck2_df_not_kr, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/asp-n_1missed_8-37pep_ck2_nokr.csv")


#******Glu-C ANALYSIS FOR ONE MISSED CLEAVAGES****
digest <- cleave(sequences, enzym="glutamyl endopeptidase", missedCleavages=1)

#Statistics for CK2 sites
#all digestion peptides
a <- lapply(digest, function(x){lapply(x, function(y){nchar(y)})})
length(unlist(a))
a <- unlist(a)
hist(a)
hist(a[a<=50])

#specific length digested peptides
int_pept <- lapply(digest, function(x){lapply(x, function(y){y[which(nchar(y) >= 8 & nchar(y) <= 37)]})})
length(unlist(int_pept))
int_pept_df <- data.frame(name = names(unlist(int_pept)), peptide = unlist(int_pept))
rownames(int_pept_df) <- NULL
dim(int_pept_df)
write.csv(int_pept_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_1missed_8-37pep.csv")

int_pept_ms <- unlist(int_pept)
int_pept_ms_length <- sapply(int_pept_ms, function(x){nchar(x)})
hist(int_pept_ms_length)

int_pept_ms_ck2 <- int_pept_ms[grep("[ST][DE][A-Z][DE]", int_pept_ms)]
length(int_pept_ms_ck2)
int_pept_ms_ck2_df <- data.frame(name = names(int_pept_ms_ck2), peptide = int_pept_ms_ck2)
rownames(int_pept_ms_ck2_df) <- NULL
dim(int_pept_ms_ck2_df)
write.csv(int_pept_ms_ck2_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_1missed_8-37pep_ck2.csv")


int_pept_ms_ck2_df_not_kr <- int_pept_ms_ck2[grep("[KR]", int_pept_ms_ck2, invert = TRUE)]
int_pept_ms_ck2_df_not_kr <- data.frame(name = names(int_pept_ms_ck2_df_not_kr), peptide = int_pept_ms_ck2_df_not_kr)
rownames(int_pept_ms_ck2_df_not_kr) <- NULL
dim(int_pept_ms_ck2_df_not_kr)
write.csv(int_pept_ms_ck2_df_not_kr, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_1missed_8-37pep_ck2_nokr.csv")

#******Glu-C ANALYSIS FOR ZERO MISSED CLEAVAGES****
digest <- cleave(sequences, enzym="glutamyl endopeptidase", missedCleavages=0)

#Statistics for CK2 sites
#all digestion peptides
a <- lapply(digest, function(x){lapply(x, function(y){nchar(y)})})
length(unlist(a))
a <- unlist(a)
hist(a)
hist(a[a<=50])

#specific length digested peptides
int_pept <- lapply(digest, function(x){lapply(x, function(y){y[which(nchar(y) >= 8 & nchar(y) <= 37)]})})
length(unlist(int_pept))
int_pept_df <- data.frame(name = names(unlist(int_pept)), peptide = unlist(int_pept))
rownames(int_pept_df) <- NULL
dim(int_pept_df)
write.csv(int_pept_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_0missed_8-37pep.csv")

int_pept_ms <- unlist(int_pept)
int_pept_ms_length <- sapply(int_pept_ms, function(x){nchar(x)})
hist(int_pept_ms_length)

int_pept_ms_ck2 <- int_pept_ms[grep("[ST][DE][A-Z][DE]", int_pept_ms)]
length(int_pept_ms_ck2)
int_pept_ms_ck2_df <- data.frame(name = names(int_pept_ms_ck2), peptide = int_pept_ms_ck2)
rownames(int_pept_ms_ck2_df) <- NULL
dim(int_pept_ms_ck2_df)
write.csv(int_pept_ms_ck2_df, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_0missed_8-37pep_ck2.csv")


int_pept_ms_ck2_df_not_kr <- int_pept_ms_ck2[grep("[KR]", int_pept_ms_ck2, invert = TRUE)]
int_pept_ms_ck2_df_not_kr <- data.frame(name = names(int_pept_ms_ck2_df_not_kr), peptide = int_pept_ms_ck2_df_not_kr)
rownames(int_pept_ms_ck2_df_not_kr) <- NULL
dim(int_pept_ms_ck2_df_not_kr)
write.csv(int_pept_ms_ck2_df_not_kr, file = "/home/teresanvd/MyFolder/MyPHD/Y2017/June2017/in_silico_digestion/glu-c_0missed_8-37pep_ck2_nokr.csv")

