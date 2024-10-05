library("derfinder")
library("derfinderData")
library("GenomicRanges")
library("knitr")


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
   fileDir <- paste(pfix, parameters["fileDir", 2], sep="/")
files <- rawFiles(fileDir,
    samplepatt = "bw", fileterm = NULL
)
names(files) <- gsub(".bw", "", names(files))

   fullCov <- readRDS(paste(pfix, parameters["fullCov", 2], sep="/"))
   regionMat <- readRDS(paste(pfix, parameters["regionMat", 2], sep="/"))

#####################################################################################
# RAIL
meanCov <- Reduce("+", fullCov$chr21) / ncol(fullCov$chr21)
createBw(list("chr21" = DataFrame("meanChr21" = meanCov)),
    keepGR =
        FALSE
)
   summaryFile <- "meanChr21.bw"
sampleFiles <- files
    regionMat.rail <- railMatrix(
        chrs = "chr21", summaryFiles = summaryFile,
        sampleFiles = sampleFiles, L = as.integer(parameters["L", 2]), cutoff = as.integer(parameters["cutoff", 2]), maxClusterGap = 3000L
    )
write.csv(identical(regionMat, regionMat.rail), outputfile)
}
