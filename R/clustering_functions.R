#######################################################
## Functions used for clustering/processing of reads ##
#######################################################

#' Put barcode counts into data frame
#' @param file barcode counts returned by mhc_cluster2.py
#' 
#' 
barcode.counts2df <- function(file) {
  header <- strsplit(readLines(file, n = 1), "\t")[[1]]
  lines <- R.utils::countLines(file)
  steps <- seq(2,lines,2)
  out <- lapply(steps, function(x) {
    mat <- matrix(data = NA, nrow = 1, ncol = 3)
    mat[1, 1] <- as.character(read.table(file, skip = x - 1, nrows = 1, header = F)[[1]])  
    mat[1, 2] <- as.character(read.table(file, skip = x, nrows = 1, header = F)[[1]]) 
    mat[1, 3] <- as.character(read.table(file, skip = x, nrows = 1, header = F)[[2]]) 
    return(as.data.frame(mat))
  })
  df <- do.call("rbind", out)
  names(df) <- header
  return(df)
}


#' Computes critical abundance skew (Edgar 2016)
#' @param d nucleotide differences
#' @param a alpha parameter
#' @return critical skew(M,C)
beta <- function(d, a = 2) {
 1/(2^(a*d + 1))
}

#' demultiplex pooled fastq files uses baroces given in the header
#' 
#' @description splits a fastq file containing data of pooled samples into individual files
#' using the barocelabels provided in sequence headers. Files called 'seqs' will be deposited
#' in subfoldes out the output directory names with the barcode strings
#' @param reads path to fastq file
#' @param out path to outputfolder
#' 
#' 
demultiplex_fastq <- function(reads = "miseq_reads/DQB-Pool/merged_reads/merged_filtered.fastq", out = "miseq_reads/DQB-Pool/merged_reads/demultiplexed", outname = "seq.fastq") {
  
  suppressPackageStartupMessages(library(ShortRead))
  source("R/fasta_fastq_functions.R")
  reads_path <- reads
  
  ## check file size > 2GB
  seq_num <- as.numeric(countLines(reads_path)/4)
  steps <- floor(seq(0, seq_num, length.out = 10))
  
  for (i in 1:(length(steps) - 1)) {
  cat("\n", "Step", i, "of", length(steps) - 1)
  reads <- readFastq(reads_path)[(steps[i] + 1):steps[i + 1]]
    
  # extract header from fastq files
  header <- as.character(id(reads)) 
  
  # get barcodes
  num_cols <- stringr::str_count(header[1], pattern = ":")
  labels <- lapply(header, function(x) strsplit(x, ":")[[1]][num_cols + 1]) 
  labels <- unlist(stringr::str_replace(labels,"barcodelabel=","")) 

  # get uniques
  barcodes <- unique(labels)
 
  if (i == 1) {
  ## erase files if they exist
  silent <- lapply(barcodes, function(x) {
    if (file.exists(file.path(out, x, outname))) file.remove(file.path(out, x, outname))
    
  })
  }
   
  # write seqs to file 
  write_fa <- function(barcode) {
    ind <- which(labels %in% barcode)
    sub <- reads[ind]
    ifelse(!dir.exists(file.path(out, barcode)), 
           dir.create(file.path(out, barcode), showWarnings = F),
           FALSE)
    
    writeFastq(sub, file.path(out, barcode, outname), "a", compress = F)  
    # fastq2disk(data = sub, file = file.path(out, barcode, outname))
  }
  x <- lapply(barcodes, write_fa)
  
  }
}


#' DensityAmplicon (Sommer et al. 2013, S4)
#' @param reads observed number of reads within amplicon
#' @param proba efficiency values for alleles with non-null read numbers
#' @param logOutput: if TRUE return the logdensity, if FALSE return the density.
DensityAmplicon <- function(reads, proba, logOutput) {
  allele.indexes <- which(reads != 0L)
  return(dmultinom(x = reads[allele.indexes], prob = proba[allele.indexes],
                   log = logOutput))
}

#' Compute the loglikelihood of the dataset (Sommer et al. 2014, S4)
#' @param proba: the efficiency values for all alleles.
#' @param data: the name of the dataset.
#' @return the loglikelihood of the dataset
LoglikData <- function(proba, data) {
  densities <- apply(data, 1,
                     function(x) DensityAmplicon(x, proba, logOutput = T))
  return(sum(densities))
}

#' pool otus generated for individual samples into a single file
#' 
#' @description loads sequences of all files matching filen within subfolders of the 
#' parentfolder. The sequences are de-replicated and written to 'pooled_zotus.fa' 
#' within the parentfolder
#' 
#' @param parentfolder path
#' @param filen filename 
pool_zotus <- function(parentfolder = NULL, filen = NULL) {
  suppressPackageStartupMessages(library(ShortRead))
  library(magrittr)

  ## list directories  
  dirs <- list.dirs(parentfolder, full.names = F, recursive = F)
  
  ## read zotus
  zotus <- lapply(dirs, function(x) {
    readFasta(file.path(parentfolder, x, filen)) %>%
      sread() %>%
      as.character()})
  variants <- unlist(zotus) %>%
    unique()
  
  
  sink(file = paste0(parentfolder,filen,"_pooled_zotus.fa"))
  for (i in 1:length(variants)) {
    cat(paste0(">", "Zotu", i, "\n"))
    cat(variants[i], "\n")
  }
  sink()
}

# This function reads and formats clustering results and exports them as RData 
process_otus <- function(locus = c('dqb', 'drb'), fname = NULL) {
  library(magrittr)
  
  locus <- match.arg(locus)
  
  ## load mapping file
  if (locus == 'dqb') {
    load("miseq_reads/DQB-Pool/RData/dqb_mapping_file.RData")
    mapping_file <- dqb_mapping_file
  } else {
    load("miseq_reads/DRB-Pool/RData/drb_mapping_file.RData")  
    mapping_file <- drb_mapping_file
  }
  
  mapping_file$pos <- as.character(mapping_file$pos)
  
  ## load otus
  otu_tab <-
    read.table(paste0(fname,".otu_table.txt"), header = T)
  
  ## extract barcodes
  otu_barcode <-
    as.character(names(otu_tab)[-1])
  
  ## sort barcodes with respect to the OTU table
  mapping_file <-
    mapping_file[match(otu_barcode, mapping_file$BarcodeSequence),]
  
  ## count reads per barcode
  x <- 
    apply(otu_tab[,-1], 2, sum) %>%
    as.data.frame() %>%
    cbind(rownames(.),.) %>%
    set_colnames(., c("BarcodeSequence", "mapped"))
  
  ## format otu table replace OTU barcodes by rack location
  names(otu_tab)[-1] <- 
    as.character(mapping_file$pos)
  
  # sort rows by allele number
  otu_number <- lapply(otu_tab$OTUId, function(x) strsplit(as.character(x), split = "Zotu")[[1]][2])
  otu_number <- as.numeric(unlist(otu_number))
   
  otu_tab <- 
    otu_tab[order(otu_number),] # match(paste0("Zotu", as.character(1:nrow(otu_tab))), otu_tab[,1]),
  rownames(otu_tab) <- as.character(otu_tab$OTUId)
  otu_tab <- otu_tab[,-1]
  
  ## summarise barcode frequency
  bc_counts <-
    barcode.counts2df(file = paste0(fname,".barcode.counts.txt")) %>%
    dplyr::left_join(x = ., y = mapping_file, by = "BarcodeSequence") %>%
    dplyr::left_join(x = ., y = x, by = "BarcodeSequence")
  
  bc_counts[,c(2:3,9)] <-  
    apply(bc_counts[,c(2:3,9)], MARGIN =  2, FUN =  function(x) as.numeric(as.character(x)))
  
  otu_freq <- 
    data.frame(allele = rownames(otu_tab),
               freq = apply(otu_tab, 1, sum))
  
  otu_freq <- otu_freq[order(otu_freq$freq, decreasing = T),]
  otu_freq$allele <- factor(otu_freq$allele, levels = as.character(otu_freq$allele))
  
  
  # normalise otu
  otu_tab2 <- apply(otu_tab, 2, function(x) x/sum(x))
  
  ## create list
  out <- list(otu_tab = otu_tab,
              otu_freq = otu_freq,
              bc_counts = bc_counts)
}

# count alleles per sample 
count_alleles <- function(x = NULL,
                          type = c("freq", "size")) {
  type = match.arg(type)
  if (type == "freq") {
    freq_x <- x/max(x)
    keep <- which(freq_x > 0.25)
  } else {
    keep <- which(x > 1250)
  }
  length(keep)
}

total_cov <- function(x,
                      type = c("freq", "size")) {
  type = match.arg(type)
  if (type == "freq") {
    freq_x <- x/max(x)
    keep <- which(freq_x > 0.25)
  } else {
    keep <- which(x > 1250)
  }
  sum(x[keep])
}

assign_alleles <- function(mat, transpose = FALSE,
                           type = c("freq", "size")) {
  type = match.arg(type)
  
  allele_list <- lapply(1:ncol(mat), function(fx) {
    if (type == "freq") {
      freq_x <- mat[,fx]/max(mat[,fx])
      keep <- which(freq_x > 0.25)
    } else {
      keep <- which(mat[,fx] > 1250)
    }
    rownames(mat)[keep]
  })
  names(allele_list) <- colnames(mat)
  allele_list
}

erase_noise <- function(matrix,
                        type = c("freq", "size")) {
  type = match.arg(type)
  allele_list <- lapply(1:ncol(matrix), function(fx) {
    
    if (type == "freq") {
      freq_x <- matrix[,fx]/max(matrix[,fx])
      keep <- which(freq_x > 0.25)
    } else {
      keep <- which(matrix[,fx] > 1250)
    }
    keep
  })
  
  # set noise to zero
  out <- matrix
  
  for (i in 1:ncol(out)) {
    out[-allele_list[[i]],i] <- 0
  }
  return(out)
}

#' hamming distance among pairs of sequences
#' @param  char character vector of sequences
hamming.distance <- function(char) {
  # helper function
  pair.diff <- function(v1,v2) {
    diff <- 0
    for (i in 1:length(v1))  diff <- diff + ifelse(v1[i] == v2[i], 0, 1)
    diff
  }
# split strings and convert to matrix  
x <- lapply(char, function(x) strsplit(x, "")[[1]]) 
x <- t(as.data.frame(x))
mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
rownames(mat) <- names(char)
colnames(mat) <- names(char)

for (k in 1:(nrow(mat) - 1)) {
  for (l in (k + 1):nrow(mat)) {
         mat[k,l] <- pair.diff(x[k,], x[l,])
         mat[l, k] <- mat[k, l]
       }
}
# set above diagonal to NA

for (i in 1:nrow(mat)) mat[i, i:ncol(mat)] <- NA

return(mat)
}


#' calculate polymorphism per site (Reche & Reinherz 2003)
#' @param vec_x vector of amino acids across sequences per site
shannon_entropy <- function(vec_x) {
  # count amino acid occurrences
  aminos <- summary(as.factor(vec_x[!is.na(vec_x)]))
  
  # estimate proportions
  prop <- aminos/sum(aminos)
  
  # estimate variation
  var <- sum(prop*log2(prop))*-1
  return(var)
} 



#' remove files created by mhc_cluster2.py for demultiplexed samples and only keep
#' raw sequences and list of fixed otus
remove_mhc_cluster_files <- function(parentfolder = "miseq_reads/DQB-Pool/merged_reads/demultiplexed/",
                                     fastq = "seq.fastq") {
  dirs <- list.dirs(parentfolder, full.names = T, recursive = F)  
  deleteR <- function(dir, keep_pattern = "\\.fixed.otus.fa", fastq = fastq) {
    files_keep <- c(list.files(dir, pattern = keep_pattern), fastq)
    files_all <- list.files(dir)
    files_remove <- subset(files_all, !(files_all %in% files_keep))
    file.remove(file.path(dir, files_remove))
  }
  out <- lapply(dirs, deleteR, fastq = fastq)
}  



# subsample and call alleles
rarefaction <- function(m = NULL, n = NULL, gain = 0.05, doc_min = 40, depth_min = 0.7) {
  out <- 
    rrarefy(m, n) %>%
    apply(., 1, function(x) {
      out <- get_genotypes(x, gain = gain, doc_min = doc_min, depth_min = depth_min)
      length(out[["alleles"]])
    })
  out
}


purify_otus <- function(x, y) {
  sample <- rownames(x)
  for (i in 1:length(sample)) {
    purified <- rep(0, ncol(x))
    purified[which(colnames(x) %in% y[[sample[i]]])] <- 
      x[which(rownames(x) == sample[i]),y[[sample[i]]]]
    x[i,] <- purified
  }
  x
}


run_genotyping <- function(fname = NULL, locus = "dqb", gain = 0.05, doc_min = 40, depth_min = 0.7) {
  
  ## get data
  mhc_data <- process_otus(locus = locus, fname = fname)
  
  ## call genotypes
  genotypes_list <- apply(mhc_data$otu_tab, 2,
                          get_genotypes,
                          names = rownames(mhc_data$otu_tab),
                          gain = gain,
                          doc_min = doc_min, 
                          depth_min = depth_min) 
  
  ## count alleles and discard bad amplicons
  allele_num_df <- lapply(genotypes_list, function(x) x[["df"]]) %>%
    do.call("rbind",.) %>%
    subset(., quality == "High")
  allele_num_df$row <- rownames(allele_num_df)
  
  ## subset otu table using amplicon quality classification
  otu_table <- 
    mhc_data$otu_tab[rownames(allele_num_df[allele_num_df$quality == "High",])] %>%
    t()
  
  ## call alleles
  genotypes <- apply(otu_table, 1, function(x) {
    out <- get_genotypes(x)
    out[["alleles"]]
  })
  
  ## unlist and summarise by observed abundance
  df <- unlist(genotypes) %>%
    as.factor() %>%
    summary()
  
  ## order alleles 
  allele_order <- lapply(names(df), function(x) strsplit(x, split = "Zotu")[[1]][2]) %>%
    unlist()  %>%
    as.numeric() %>%
    order(., decreasing = T)
  
  df <- data.frame(x = names(df), y = df)
  df$x <- factor(df$x, levels = df$x[allele_order])
  df$z <- ifelse(df$y == 1,  "Putative Artefact", "Putative Allele")
  
  zotus <- df
  temp <- summary(as.factor(zotus$z))
  
  ## check if there is no artefact
  if (length(temp) == 1 && names(temp) == "Putative Allele") {
    temp <- c(temp[1], 0)
  } 

  zotu_summary <- data.frame(putatitive_artefact = temp[2],
                             putatitive_allele = temp[1],
                             Artefact  = nrow(mhc_data$otu_tab) - sum(temp))
  
  output <- list(called_alleles = genotypes,
                 zotus = zotus,
                 zotu_summary = zotu_summary)
  
  return(output)
  
}
