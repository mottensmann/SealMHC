#' convert fasta file to matrix
#' @description reads and aligns sequences within a single fasta file; then converts to
#' a matrix
#' @param fasta file contaning otus in fasta format
#' @return matrix with aligned sequences
fasta2mat <- function(fasta = "miseq_reads/DQB-Pool/clustered_reads/dqb_pct_1.0_a_1.0.fixed.otus.fa") { 
  
  library(magrittr)
  library(ShortRead)
  library(stringr)
  
  # read fasta sequences
  data <- readFasta(fasta)
  
  # shorten names
  short_name <- 
    as.character(id(data)) %>%
    lapply(., function(x) str_replace(x, pattern = "otu", replacement = "")) %>%
    unlist()
  
  
  # Extract sequence and split
  seqs <- 
    sread(data) %>%
    clustal(.) %>%
    set_rownames(.,short_name) %>%
    as.matrix()
  
  return(seqs)
}

#' export sequences as a fasta file
fastq2disk = function(data = NULL, file = NULL) {
  if (file.exists(file)) {
    sink(file, append = F)
    sink()
  }
  ShortRead::writeFastq(object = data,
                        file = file,
                        full = T,
                        compress = F,
                        mode = "a")
}

fasta2disk = function(data = NULL, file = NULL) {
  if (file.exists(file)) {
    sink(file, append = F)
    sink()
  }
  ShortRead::writeFasta(object = data,
                        file = file)
}


#' subset fastq file 
#' 
#' @param read
#' path to read
#' 
#' @param n
#' number of reads to keep
#' 
#' @param out
#' suffix of output
#' 
#' @export
#'
fastq_subset <- function(read = "miseq_reads/DQB-Pool/parsed_barcodes/reads1.fastq", n = 1000000, out = "sub") {
  
  loadR <- function(x, n = n, out = out) {
    source("R/fastq2disk.R")
    fastq <- ShortRead::readFastq(x)
    cat("\nWrite to file...\n")
    out_name <- paste0(strsplit(x, ".fastq")[[1]][1],"_",out,".fastq")
    fastq2disk(data = fastq[1:n], file = out_name)
  }
  
  cat("Loading...\n\t", read)
  loadR(x = read, n = n, out = out)
  cat("done")  
  
} 

