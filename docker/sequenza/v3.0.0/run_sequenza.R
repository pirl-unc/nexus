# The purpose of this Rscript is to run Sequenza for a pair of tumor and normal samples


library(sequenza)
library(copynumber) # necessary for hg38
library(optparse)


# Step 1. Parse input parameters
options.list <- list(
    make_option(c("-s", "--sample-id"),
                type="character",
                dest="sample.id",
                help="Sample ID"),
    make_option(c("-i", "--small-seqz-file"),
                type="character",
                dest="small.seqz.file",
                help="Input seqz file"),
    make_option(c("-a", "--assembly"),
                type="character",
                dest="assembly",
                help="Genome assembly ('hg19' or 'hg38')"),
    make_option(c("-c", "--chromosomes"),
                type="character",
                dest="chromosomes",
                help="Chromosomes (e.g. 'chr1,chr2,chr3'"),
    make_option(c("-o", "--output-path"),
                type="character",
                dest="output.path",
                help="Output path")
)
parser <- OptionParser(usage="%prog [options]",
                       option_list = options.list)
args <- parse_args(parser, args = commandArgs(trailingOnly=TRUE))

# Step 2. Run Sequenza
chromosome.list <- strsplit(args$chromosomes, split = "\\,")[[1]]
test <- sequenza.extract(args$small.seqz.file, verbose = TRUE,
                         assembly = args$assembly,
                         chromosome.list = chromosome.list)
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test,
                 cp.table = CP,
                 sample.id = args$sample.id,
                 out.dir = args$output.path)
