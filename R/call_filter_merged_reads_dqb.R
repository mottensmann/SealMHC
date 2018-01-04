source("R/filter_merged_reads.R")
filter_merged_reads(reads = "miseq_reads/DQB-Pool/merged_reads/merged.fastq",
                    barcodes = "miseq_reads/DQB-Pool/parsed_barcodes/barcodes.fastq",
                    mapping_file = "miseq_reads/DQB-Pool/dqb_barcodes.txt",
                    forward_primer = "GCTGTTGGTTGGGCTGAG",
                    reverse_primer = "CCACCTCAGCAGGAACAGTG",
                    max.mismatch = 0,
                    with.indels = F,  
                    suffix = "_filtered",
                    illumina_adapter = "CTTGTA",
                    pcr_size = 459)

