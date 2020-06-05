#Assignment 3 Part 1
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv", destfile = "gene.expression.tsv")

#Question 1:Read in the file, making the gene accession numbers the row names. Show a table of values for the first six genes. 
x<- read.table("gene.expression.tsv")
x <-read.table("gene.expression.tsv",header= TRUE, row.names = 1)
x[1:7, ]
head(x)
str(x)

#Question 2:Make a new column which is the mean of the other columns. Show a table of values for the first six genes. 
x$Mean1 <- rowMeans(x)
head(x)

#Question 3:  List the 10 genes with the highest mean expression 
order(-x$Mean1)
x[order(-x$Mean1), ]
x_sorted <- x[order(-x$Mean1), ]
head(x_sorted,10)

#Question 4:  Determine the number of genes with a mean <10 
a <- subset(x, Mean1<10) 
head(a)
str(a)

#Question 5: Make a histogram plot of the mean values in png format and paste it into your report.
hist(x$Mean1)
hist(x$Mean1, breaks = 20, border="red", col = "green")

#Question 6: Import this csv file into an R object. What are the column names?
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv", destfile = "growth_data.csv")
y<- read.csv("growth_data.csv", header = TRUE)
head(y)
str(y)
#The column names are Site, TreeID, Circumf_2004_cm, Circumf_2009_cm, Circumf_2014_cm and Circumf_2019_cm.

#Question 7: Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites.
mean(y$Circumf_2004_cm)
subset(y, Site=="northeast")
head(y)
str(y)
subset(y, Site=="southwest")
head(y)
str(y)
ne <- subset(y, Site=="northeast")
head(ne)
str(ne)
sw <- subset(y, Site=="southwest")
head(sw)
str(sw)
mean(ne$Circumf_2004_cm)
mean(sw$Circumf_2004_cm)
mean(ne$Circumf_2019_cm)
mean(sw$Circumf_2019_cm)
sd(sw$Circumf_2019_cm)
sd(ne$Circumf_2019_cm)
sd(sw$Circumf_2004_cm)
sd(ne$Circumf_2004_cm)
head(ne)
head(sw)

# Question 8: Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm)
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm,sw$Circumf_2004_cm,sw$Circumf_2019_cm,names= c("ne2004","ne2019","sw2004","sw2019"), ylab="Circumference (cm)", xlab="site and years", main="Growth at plantation site")






#Assignment 3 Part 2

#downloading library
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1: Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. Use themakeblast() function to create a blast database. How many sequences are present in the E.coli set?
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")


R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl","-parse_seqids")
# 4140 sequences are present in the E.coli set

#Question 2: Download the sample fasta sequences and read them in as above. For your allocated sequence,determine the length (in bp) and the proportion of GC bases.
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile= "sample.fa")
s <- read.fasta("sample.fa") # "s" is the object for the sample fasta file.
head(s)
str(s)
#Allocated sequence is "25"
c <- s[[25]] # "c" is the object for our group sequence 25
c
seqinr::getLength(c)

seqinr::GC(c)
#length for sequence"25" is 1047 and proportion of GC bases is 0.5635148

#Question 3: You will be provided with R functions to create BLAST databases and perform blast searches. Use blast to identify what E. coli  gene your sequence matches best. Show a table of the top 3 hits including percent identity, E-value and bit scores. 
myblastn <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(myseq,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 1 -evalue 0.05 -ungapped >"
                       , "tmp"))
  
  res <- NULL
  if (file.info("tmp")$size > 0 ) {
    res <- readLines("tmp")
  }
  
  res
}

myblastn_tab <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(myseq,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 6 -evalue 0.05 -ungapped >"
                       , mytmpfile2))
  
  res <- NULL
  if (file.info(mytmpfile2)$size > 0 ) {
    res <- read.csv(mytmpfile2,sep="\t",header=FALSE)
    colnames(res) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                       "qstart","qend","sstart","send","evalue","bitscore")
  }
  unlink(c(mytmpfile1,mytmpfile2))
  if (!is.null(res)  ) {
    res <- res[order(-res$bitscore),]
  }
  res
}
mutator <- function(myseq,nmut) {
  myseq_mod <- myseq
  mypos<-sample(seq_along(myseq),nmut)
  myseq_mod[mypos] <- sample(c("a","c","g","t"),length(mypos),replace = TRUE)
  return(myseq_mod)
}

res <- myblastn_tab(myseq = c, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(res)
str(res)

mysequence <- as.character(res$sseqid[1:1047])


View(res)
res
#The E.coli gene which matches our sequence the best is "CP4-6 prophage; putative ferric transporter subunit"

#Question 4:You will be provided with a function that enables you to make a set number of point mutations to your sequence of interest. Run the function and write an R code to check the number of mismatches between the original and mutated sequence. 
C25_mut <-mutator(myseq = c, 100) #C25_mut is the object for the mutated group sequence for 100 mutations.

C25_mut

myblastn_tab(C25_mut,db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
res <- myblastn_tab(C25_mut,db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(res)

#Pairwise alignment created
C25_mut_1 <- DNAString(c2s(C25_mut)) #"C25_mut_1" is the DNA string object of "C25_mut"

c_1 <- DNAString(c2s(c)) #"c_1" is the DNA string object of "c"

aln <- pairwiseAlignment(c_1,C25_mut_1) #aln is the object of the pairwise alignment between "c_1 & C25_mut_1"

#Calculation of percent sequence identity for "aln"
pid(aln)

nmismatch(aln) #The number of mismatches found were 57.