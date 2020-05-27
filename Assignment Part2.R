
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")



download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")


R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl","-parse_seqids")
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile= "sample.fa")
s <- read.fasta("sample.fa")
head(s)
str(s)
c <- s$`25`
head(c)
str(c)

seqinr::getLength(c)

seqinr::GC(c)

myblastn <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(myseq,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 1 -evalue 0.05 -ungapped >" 
                       , "tmp"))
  # probably need to add headers here
  res <- NULL
  if (file.info("tmp")$size > 0 ) {
    res <- readLines("tmp")
  }
  #unlink(c(mytmpfile1,mytmpfile2))
  #res <- res[order(-res$bitscore),]
  res
}

myblastn_tab <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(myseq,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 6 -evalue 0.05 -ungapped >" 
                       , mytmpfile2))
  # probably need to add headers here
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
V


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
#search ensembl Bacteria -> ecoli AAC73365 in google and include name of gene in report

