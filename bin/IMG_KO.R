################################################################################
#
# IMG Genome KEGG KO:
#       Find and Organize KEGG KOs for Genomes in IMG
#
# Written by: Mario Muscarella
#
# Last Update: 2017 04 17
#
# Goals:
#       1.
#       2.
#       3.
#
################################################################################

# System Specific Information
sys <- Sys.info()
if(sys[1] == "Darwin"){
  setwd("~/GitHub/microbial-innovations/bin")
  output_location <- "~/Desktop/output/"
} else {
  setwd("~/microbial-innovations/bin")
  output_location <- "~/JGI/output/"
}

# Genome Location
output_location <- "~/JGI/output/"


# Import Genome IDs
genomes.raw <- list.files(path = paste(output_location, "img/", sep = ""), 
                          pattern = "tar.gz")
genomes <- as.data.frame(matrix(NA, ncol = 2, nrow = length(genomes.raw)))
colnames(genomes) <- c("Genome", "OID")
for (i in 1:length(genomes.raw)){
  temp <- unlist(strsplit(genomes.raw[i], split = "[.]"))
  genomes$Genome[i] <- gsub("_FD", "", temp[1])
  genomes$OID[i] <- gsub("^[0-9]*\\-", "", temp[2])
}

# Import KEGG Database
gene.database.temp <- read.table("../data/ko00001.keg", header = T,
                                 fill = T, comment.char = "#",
                                 quote = "", sep = "\t")

# Oranize KEGG Genes and KO
D.level.index <- grep("D      K", gene.database.temp$X.D)
D.levels <- gsub("D      K", "K",
                 as.character(gene.database.temp$X.D[D.level.index]))
D.ko <- sapply(strsplit(D.levels , "  "), function(x) x[1])
D.gene <- sapply(strsplit(D.levels , "  "), function(x) x[2])
KEGG.KO <- sort(unique(D.ko), decreasing = F)


# Create KEGG Database for Genomes
files <- list.files(path = output_location, pattern = "ko.tab.txt", recursive= T)
img.kegg <-data.frame(matrix(NA, ncol = length(KEGG.KO) + 2, nrow = length(files)))
colnames(img.kegg) <- c("Genome", "OID", KEGG.KO)

for (i in 1:length(files)){
  temp <- read.delim(paste(output_location, files[i], sep = ""), quote = "")
  temp.ko <- gsub("KO:", "", temp[,10])
  abund.ko <- table(temp.ko)
  genome <- genomes[which(genomes$OID == unlist(strsplit(files[i], "/"))[1]), ]
  print(dim(genome))
  if(dim(genome)[1] == 0){
    genome <- c("Genome" = NA, "OID" = NA)
  }
  print(genome)
  img.kegg[i, 1:2] <- as.vector(genome)
  img.ko <- as.vector(abund.ko[match(KEGG.KO, names(abund.ko))])
  img.ko[is.na(img.ko)] <- 0
  print(dim(img.ko))
  img.kegg[i , 3:dim(img.kegg)[2]] <- img.ko
}

# Save KEGG Database
write.table(img.kegg, file = "../data/JGI.KEGG.txt", quote = F, sep = "\t",
            row.names = F, col.names = T)

# Create 16S Database for Genomes
files <- list.files(path = output_location, pattern = ".gff", recursive= T)
img.16S <-data.frame(matrix(NA, ncol = 3, nrow = length(files)))
colnames(img.kegg) <- c("Genome", "OID", "16S")

for (i in 1:length(files)){
  temp <- read.delim(paste(output_location, files[i], sep = ""), quote = "",
                     head = F, skip = 1, sep = "\n")
  grep("rRNA.*product=16S", temp, perl = T)
  grep("rRNA.*product=16S", readLines(paste(output_location, files[i], sep = "")))
  temp[817 - 1, ]
  temp[898 - 1, ]
