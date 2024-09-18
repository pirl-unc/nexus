# The purpose of this Rscript is to run Bambu for one long-read RNA-seq sample


library(bambu)
library(optparse)


# Step 1. Parse input parameters
options.list <- list(
  make_option(c("-i", "--bam-file"),
              type="character",
              dest="bam.file",
              help="Input BAM file"),
  make_option(c("-f", "--fasta-file"),
              type="character",
              dest="fasta.file",
              help="Input FASTA file"),
  make_option(c("-g", "--gtf-file"),
              type="character",
              dest="gtf.file",
              help="Input GTF file"),
  make_option(c("-o", "--output-path"),
              type="character",
              dest="output.path",
              help="Output path")
)
parser <- OptionParser(usage="%prog [options]",
                       option_list = options.list)
args <- parse_args(parser, args = commandArgs(trailingOnly=TRUE))

# Step 2. Input validation
if (!file.exists(args$bam.file)) {
  stop("BAM file not found!")
}
if (!file.exists(args$fasta.file)) {
  stop("FASTA file not found!")
}
if (!file.exists(args$gtf.file)) {
  stop("GTF file not found!")
}
if (!dir.exists(args$output.path)) {
  dir.create(args$output.path, recursive = TRUE)
}

# Step 3. Run Bambu
bambuAnnotations <- prepareAnnotations(args$gtf.file)
se <- bambu(reads = args$bam.file,
            annotations = bambuAnnotations,
            genome = args$fasta.file,
            discovery = TRUE,
            quant = TRUE)
writeBambuOutput(se, path = args$output.path)
