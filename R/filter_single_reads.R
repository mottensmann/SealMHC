#' Filters unpaired FASTQ files prior to clustering sequencing into OTUs
#' 
#' @param reads
#' path to a file containing raw reads in fastq format with file ending .fastq
#' 
#' @param barcodes
#' path to a file containing barcodes for \code{reads} with file ending .fastq
#' 
#' @param mapping_file
#' mapping file that contains barcode sequences for all samples
#' 
#' @param forward_primer
#' character giving the forward primer sequence in 5'-3' orientation
#' 
#' @param reverse_primer
#' character giving the reverse primer sequence in 3'-5' orientation 
#' 
#' @param max.mismatch
#' Interger giving the maximum number of allowed mismatches between sequences and primer. By default 1.
#' 
#' @param with.indels
#' Boolean. By default no indels in primer sequences are allowed.
#' 
#' @param suffix 
#' string added to file names for generating the output. Note, existing files may be overwritten.
#' 
#' @param illumina_adapter
#' illumina adapter sequence written at the end of the identifier line. Same for all samples.
#' 
#' @param pcr_size
#' expected fragment length. If specified all reads of differing size will be filtered out.
#' 
#' @param min_size
#' minimum fragment length
#' 
#' @export
#' 
filter_single_reads <- function(reads = NULL, barcodes = NULL, mapping_file = NULL, forward_primer = NULL, reverse_primer = NULL, max.mismatch = 1, with.indels = F, suffix = "_filtered", illumina_adapter = "GCCAAT", pcr_size = NULL, min_size = NULL, splits = 100) {
  
  suppressPackageStartupMessages(library(ShortRead))
  #source("R/fasta_fastq_functions.R")
  
  reads_path <- reads
  barcodes_path <- barcodes
  
  ## check file size > 2GB
  cat("Load data ...")
  seq_num <- as.numeric(countLines(reads_path)/4)
  seq_num2 <- as.numeric(countLines(barcodes_path)/4)
  cat(seq_num, "reads\n")
  if (seq_num2 > seq_num) {
    cat("Subset barcodes based on reads ...")
    reads <- readFastq(reads_path)
    reads <- id(reads)
    barcodes <- readFastq(barcodes_path)
    barcodes <- barcodes[id(barcodes) %in% reads]
    rm("reads")
    barcodes_path <- strsplit(barcodes_path, "/")[[1]]
    barcodes_path <- 
      file.path(paste0(barcodes_path[1:length(barcodes_path) - 1], collapse = "/"), "barcodes_filt.fastq")
    if (file.exists(barcodes_path)) file.remove(barcodes_path)

    writeFastq(object = barcodes,
               file = barcodes_path,
               compress = F)
    rm("barcodes")
    rm("seq_num2")
  }
  
  steps <- floor(seq(0, seq_num, length.out = splits + 1))
  
  # read mapping file
  mapping <- read.table(file = mapping_file, header = T)[,2]
  
  ## define output names
  make_name <- function(x, suffix = suffix) {
    path <- strsplit(x, "/")[[1]]
    name <- strsplit(path[length(path)], ".fastq")[[1]]
    prefix <- paste0(path[-length(path)], collapse = "/")
    return(paste0(prefix, "/", name, suffix, ".fastq"))
  }
  
  out_names <- unlist(lapply(X = c(reads_path, barcodes_path), make_name, suffix))
  out <- lapply(out_names, function(x) {
    if (file.exists(x)) file.remove(x)
  })
  
  rm(list = c("seq_num", "mapping_file", "suffix"))
  
  for (i in 1:(length(steps) - 1)) {
  cat("\n", "Step", i, "of", length(steps) - 1)
    
  cat("\nRead files ... ")
  reads <- readFastq(reads_path)[(steps[i] + 1):steps[i + 1]]
  barcodes <- readFastq(barcodes_path)[(steps[i] + 1):steps[i + 1]]
  cat("done\n")
  
  # extract sequences from fastq files
  cat("Extract sequences ... ")
  reads_seq = sread(reads)
  barcodes_seq <- sread(barcodes)
  cat("done\n")
  
  # search for primer
  if (!is.null(forward_primer)) {
    cat("Extract forward primer region ... ")
    primer_f = DNAStringSet(substr(reads_seq, 1, nchar(forward_primer)))  
    cat("done\nDetect matches with the primer sequence ... ")
    hits_f = vcountPattern(forward_primer, 
                           primer_f, 
                           max.mismatch = max.mismatch,
                           with.indels = with.indels)
    hits <- as.logical(hits_f)
    rm(list = c("primer_f", "hits_f"))
    cat("\n kept", round(length(hits[hits == TRUE])/length(hits)*100, 2), "%\n")
  } else {
    cat("Extract reverse primer region ... ")
    # search for reverse primer in the reverse complement
    cat("Extract reverse primer region ... ")
    primer_r = DNAStringSet(substr(reverseComplement(reads_seq), 1, nchar(reverse_primer))) 
    cat("done\n")
    cat("Detect matches with the primer sequence ... ")
    hits_r = vcountPattern(reverse_primer,
                           primer_r, 
                           max.mismatch = max.mismatch,
                           with.indels = with.indels)
    cat("done\n")
    hits <- as.logical(hits_r)
  }
  
  # filter the reads
  reads <- reads[hits]
  
  # subset barcodes file based on ids matching the reads
  barcodes <- barcodes[id(barcodes) %in% id(reads)]
  
  cat("\nAdd barcodes to sequence headers ... ")
  
  # clean workspace
  
  # for (i in 1:(length(steps) - 1)) {
  #   cat("\n", "Step", i, "of", length(steps) - 1)
    # prepare lablels
    x <- as.character(reads@id)

    # remove whitespaces
    x <- sub(pattern = " ", replacement = "", x = x)
    
    # subsitute illumina label by 'barcodelabel='
    x <- sub(pattern = illumina_adapter, replacement = "barcodelabel=", x = x)
    x <- paste0(x, as.character(barcodes@sread))
    
    ## polish header
    reads <- 
      ShortReadQ(sread(reads), quality(reads), BStringSet(x))
    
    barcodes <- 
      ShortReadQ(sread(barcodes), quality(barcodes), BStringSet(x))

    #  extract sequences
    barcodes_seq <- sread(barcodes)
    
    ## match barcodes to samples
    indices <- which(barcodes_seq %in% mapping)
    
    reads <- reads[indices]
    barcodes <- barcodes[indices]
    
    writeFastq(reads, out_names[1] , "a", compress = F)  
    writeFastq(barcodes, out_names[2] , "a", compress = F)    
    
    rm(list = c("barcodes", "barcodes_seq",
                "reads", "reads_seq",
                "x", "indices", "hits"))

  }
  
}
 


