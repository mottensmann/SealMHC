#' Filters merged reads prior to (i) clustering Zotus and (ii) assigning raw reads to them
#' 
#' @description 
#' This script does the following tasks:
#' 1.) Ensure that reads contain expected primer pair
#' 2.) Discard reads with unexpected length
#' 3.) Adds barcodes to headers of the fastq file
#' 4.) Removes reads that contain unexpected barcodes
#' 5.) Export subsetted to new FASTQ files.
#' 
#' @param reads
#' path to a file containing merged reads in fastq format with file ending '.fastq'
#' 
#' @param barcodes
#' path to a file containing barcodes for \code{reads} with file ending '.fastq'
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
#' Interger giving the maximum number of allowed mismatches between sequences and primer. By default 1
#' 
#' @param with.indels
#' Boolean. By default no indels in primer sequences are allowed.
#' 
#' @param suffix 
#' string added to file names for generating the output. Note, existing files may be overwritten.
#' 
#' @param illumina_adapter
#' illumina adapter sequence written at the end of the identifier line. Same for all samples. This string will be replaced by 'barcodelabel=' to allow demultiplexing by usearch
#' 
#' @param pcr_size
#' expected fragment length. If specified all reads of differing size will be filtered out.
#' 
#' @param min_size
#' minimum fragment length. If specifed all reads that are shorter will be discarded
#' 
#' @export
#' 
filter_merged_reads <- function(reads = NULL, barcodes = NULL, mapping_file = NULL, forward_primer = NULL, reverse_primer = NULL, max.mismatch = 1, with.indels = F, suffix = "_filtered", illumina_adapter = "GCCAAT", pcr_size = NULL, min_size = NULL) {

  suppressPackageStartupMessages(library(ShortRead))
  source("R/fasta_fastq_functions.R")

# checks  
  if (!file.exists(reads)) stop("Specify a valid path to the reads file")
  if (!file.exists(barcodes)) stop("Specify a valid path to the barcodes file")
  if (!file.exists(mapping_file)) stop("Specify a valid path to the mapping file")
  if (!is.character(forward_primer)) stop("Give a forward primer sequence")
  if (!is.character(reverse_primer)) stop("Give a reverse primer sequece")
  
  ## keep path information for later use
  reads_path <- reads
  barcodes_path <- barcodes
  
  ## get reads and associated barcodes
  cat("Read files ... ")
  reads <- readFastq(reads)
  initial_size <- length(reads)
  barcodes <- readFastq(barcodes)
  cat("done\n")
  
  # extract sequences from fastq files
  cat("Extract sequences ... ")
  reads_seq = sread(reads)
  barcodes_seq <- sread(barcodes)
  cat("done\n")
  
  # search for forward primer
  cat("Extract forward primer region ... ")
  primer_f = DNAStringSet(substr(reads_seq, 1, nchar(forward_primer)))  
  cat("done\nDetect matches with the primer sequence ... ")
  hits_f = vcountPattern(forward_primer, 
                         primer_f, 
                         max.mismatch = max.mismatch,
                         with.indels = with.indels)
  cat("done\n")
  perc_forward <- round(sum(hits_f/length(hits_f))*100,2)
  size_forward <- sum(hits_f)
  cat(perc_forward, "% of all reads contain the expected forward primer\n")
  
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
  perc_reverse <- round(sum(hits_r/length(hits_r))*100,2) 
  size_reverse <- sum(hits_r)
  cat(perc_reverse, "% of all raw reads contain the expected reverser primer\n")
  
  # merge hits and turn into a logical
  hits <- hits_f + hits_r 
  hits <- ifelse(hits == 2, 1, 0) 
  hits <- as.logical(hits)
  
  # keep reads with bother primers present
  reads <- reads[hits]
  primers_truncated <- length(reads)
  size_primers <- length(reads)
  perc_primers <- round(size_primers/initial_size*100,2)
  cat(perc_primers, "% of all raw reads contain both expected primers\n")
  
  # size selection. Retain only reads of the expected length
  if (!is.null(pcr_size)) {
    cat("Size selection with expected n =", pcr_size," ... ")
    width <- reads@sread@ranges@width
    width.matched <- which(width %in% seq(from = floor(pcr_size*0.99), to = floor(pcr_size*1.01)))
    # apply selection
    reads <- reads[width.matched]
    cat("done\n")
    cat("Retained", round(100*(length(width.matched)/initial_size),2),"% of all reads")
  }
  
  if (!is.null(min_size)) {
    cat("Size selection with expected length >= ", min_size," ... ")
    width <- reads@sread@ranges@width
    width.matched <- which(width >= min_size)
    # apply selection
    reads <- reads[width.matched]
    #barcodes <- barcodes[width.matched]
    cat("done\n")
    cat("Retained", round(100*(length(width.matched)/initial_size),2),"% of all reads")
  }


  # get sequence header
  reads_id <- id(reads)
  barcodes_id <- id(barcodes)
  # subset barcodes based on id of reads
  barcodes <- barcodes[barcodes_id %in% reads_id]
  
  ## define output names
  make_name <- function(x, suffix = suffix) {
    path <- strsplit(x, "/")[[1]]
    name <- strsplit(path[length(path)], ".fastq")[[1]]
    prefix <- paste0(path[-length(path)], collapse = "/")
    return(paste0(prefix, "/", name, suffix, ".fastq"))
  }
  
  out_names <- unlist(lapply(X = c(reads_path, barcodes_path), make_name, suffix))
  
  cat("\nAdd barcodes to sequence headers ... ")
  
  # preapre labels
  x <- as.character(reads@id)
  # remove whitespaces
  x <- sub(pattern = " ", replacement = "", x = x)
  # subsitute illumina label by 'barcodelabel='
  x <- sub(pattern = illumina_adapter, replacement = "barcodelabel=", x = x)
  x <- paste0(x, as.character(barcodes@sread))
  
  ## polish header
  reads <- ShortReadQ(sread(reads), quality(reads), BStringSet(x))
  barcodes <- ShortReadQ(sread(barcodes), quality(barcodes), BStringSet(x))
  cat("done\n")
  
  #  extract sequences
  cat("Filter for expected barcodes ... ")
  barcodes_seq <- sread(barcodes)
  
  # read mapping file
  mapping <- read.table(file = mapping_file, header = T)[,2]
  
  ## match barcodes to samples
  indices <- which(barcodes_seq %in% mapping)

  reads <- reads[indices]
  barcodes_selected <- length(reads)
  
  ## take a look at unmatched barcodes
  unexp_barcode <- barcodes[-indices]
  
  barcodes <- barcodes[indices]
  cat("done\n")
  
 
  cat("Write to disk ... ")
  fastq2disk(data = reads, file = out_names[1])
  fastq2disk(data = barcodes, file = out_names[2])
  cat("done\n")
  cat(round(100 * barcodes_selected/initial_size,2), " % of all reads are retained\n")
  
  ## write log file
  log_out <- strsplit(reads_path, "/")[[1]]
  log_out <- paste0(log_out[-length(log_out)], collapse = "/")
  log_out <- paste0(log_out, "/log", suffix, ".txt")
  if (file.exists(log_out)) {
    sink(log_out, append = F)
    sink()
  }

  sink(log_out)
  cat("Input file:\t", reads_path, "\n\t", initial_size,"raw reads\n")
  cat("Filter for primers:","\n\t",
      size_forward, paste0("(", perc_forward,"%) matches with forward primer sequence\n\t"),
      size_reverse, paste0("(", perc_reverse,"%) matches with reverse primer sequence\n\t"),
      size_primers, paste0("(", perc_primers,"%) matches with both primer sequences\n"))
  if (!is.null(pcr_size)) {
    cat("Filter for expected size of", pcr_size,"bp:\n\t",
        length(width.matched), paste0("(", round(length(width.matched)/initial_size*100,2),"%) sequences are retained\n"))
  }
  if (!is.null(min_size)) {
    cat("Filter for expected size >", min_size,"bp:\n\t",
        length(width.matched), paste0("(", round(length(width.matched)/initial_size*100,2),"%) sequences are retained\n"))
    
  }
  cat("Filter for expected barcodes:\n\t",
      barcodes_selected, paste0("(", round(100 * barcodes_selected/initial_size,2),"%) of all raw reads are retained"))
      
  sink()
  
}