source("R/filter_single_reads.R")
filter_single_reads(reads = "miseq_reads/DRB-Pool/parsed_barcodes/reads1_qual_filt.fastq",
                    barcodes = "miseq_reads/DRB-Pool/parsed_barcodes/barcodes.fastq",
                    mapping_file = "miseq_reads/DRB-Pool/drb_barcodes.txt",
                    forward_primer = "GGTGACCGGATCCTCTCTG",
                    max.mismatch = 0,
                    with.indels = F,  
                    suffix = "_cont_filt",
                    illumina_adapter = "GCCAAT",
                    splits = 100)