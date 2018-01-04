#'Checks primer pairs for specifity
#'
#' @description 
#' primer_check inspects blast hits derived from mapping primers to a database in order to determine
#' primer pairs that may be inspecific. Investigates all instances were both primers of a pair bind to
#' the same Contig. Reported are all alignments that might yield to successful amplification
#' 
#' @param Blats
#' output of blastn in format -outfmt 6
#' 
#' @param threshold
#' Fragment size of primers paris. Only fragments below the threshold are reported
#' 
#' @export
#' 
primer_check <- function(Blast = NULL , threshold = 500) {
  
  ## read resutls
  data <- read.table(file = Blast)[,c(1,2,4,9,10)]
  names(data) <- c("primer","contig","length","start","end")
  data[["primer"]] <- as.character(data[["primer"]])
  data[["id"]] <- as.character(remove_F_or_R(data[["primer"]]))
  data[["type"]] <- as.character(primer_type(data[["primer"]]))
  
  pairs <- unique(data[["id"]])
  out <- matrix(data = NA, nrow = 1, ncol = ncol(data))
  colnames(out) <- names(data)
  for (n in 1:length(pairs)) {
    temp <- subset(data, data[["id"]] == pairs[n]) 
    candidates <- unique(temp[["contig"]][which(duplicated(temp[["contig"]]))])
    for (i in 1:length(candidates)) {
      temp2 <- temp[temp[["contig"]] == candidates[i],] 
      if (length(unique(temp2[["type"]])) == 2) {       
        if (length(dist(temp2[["start"]])[abs(dist(temp2[["start"]]) < threshold)] != 0)) { 
          out <- rbind(out, temp2)
        }
      }
    }
  }
  out <- out[-1,]
  rownames(out) <- NULL
  return(out)
}

#' @description 
#' Extracts sequences from  a bed file and allows to add flanking sequences. If segments are too
#' short for giving the requested flanking site, the maximum possible number of bases is returned.
#' 
#' @param file 
#' fasta file generated with bedtools
#' @param cols 
#' columns of the input file. No need to change defaults
#' @param flanking 
#' number of nucleotides on either site of the target
#' @param dir
#' folder containing the specified fasta file
#' 
target_extract <- function(file, cols = c(2,9,10), flanking = 150, dir = "blast") {
  if (length(strsplit(file, "transcriptome")[[1]]) == 2) flanking <- 0
  file <- paste0(dir, file)
  data <- read.table(file, header = F)
  data <- data[,cols] 
  data <- data[!duplicated(data), ]
  for (i in 1:nrow(data)) {
    if (data[i,3] < data[i,2]) { # Sort start and end positions
      x <- data[i,2]
      y <- data[i,3]
      data[i,2] <- y
      data[i,3] <- x
    }
    data[i,2] <- data[i,2] - flanking
    data[i,3] <- data[i,3] + flanking
    data[data[,2] < 0,2] <- 0
  }
  file <- strsplit(file,split = ".fasta")[[1]]
  write.table(x = data, file = paste0(file, ".bed"),
              row.names = F,col.names = F,quote = FALSE,sep = "\t")
}

#' remove '_F' or '_R' of a string
remove_F_or_R <- function(primer_list){
  for (i in 1:length(primer_list)) {
    temp <- strsplit(primer_list[i],split = "_")[[1]]
    primer_list[i] <- stringr::str_c(temp[1:(length(temp) - 1)],collapse = "_")
  }
  primer_list
}

#' get type of a primer by extracting last letter of name
primer_type <- function(primer_list) {
  for (i in 1:length(primer_list)) { 
    temp <- strsplit(primer_list[i],split = "_")[[1]]
    primer_list[i] <- temp[length(temp)]
  }
  primer_list
} 





