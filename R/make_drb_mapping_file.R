## create RData mapping file for demultiplexing clustering results
## 
## Prior to running this script, all suffixes 'mum' and 'pup' were replaced by 'M' and 'P' respectively

## read mapping file containing sampleID and BarcodeSequence
mapping_file <- as.data.frame(readxl::read_xlsx("documents/drb_sample_identifiers.xlsx", sheet = 1)[,1:2])

## Obtain original identifiers by removing suffix (e.g. '.A10') 
mapping_id <- unlist(lapply(mapping_file$SampleID, function(x) {
  t <- strsplit(as.character(x), "\\.")[[1]]
  if (length(t) == 3) t <- paste0(t[1:2], collapse = "/")  
  if (length(t) == 2) t <- t[1]
  t
})) 

## read sample identifiers
names <- as.data.frame(readxl::read_xlsx("documents/drb_sample_identifiers.xlsx", sheet = 2))

## match the order of rows
names <- names[match(mapping_id, names$ID),]

# combine files
drb_mapping_file <- cbind(names,
                          mapping_file,
                          pos = unlist(lapply(mapping_file$SampleID, function(x) {
                            temp <- strsplit(x,"\\.")[[1]]
                            ifelse(length(temp) == 2, temp[2], temp[3])
                          } )))

save(drb_mapping_file, file = "miseq_reads/DRB-Pool/RData/drb_mapping_file.RData")
