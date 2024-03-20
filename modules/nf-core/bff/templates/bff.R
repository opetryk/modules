#!/usr/bin/env Rscript

library(DropletUtils)
library(Seurat)
library(ggplot2)
library(cowplot)
if(!require("cellhashR")){
    devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')
    library("cellhashR")
}

library(cellhashR)
library(here)
library(dplyr)
library(argparse)

if(!require("tidyverse")){
      install.packages("tidyverse")
      library("tidyverse")
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
    fileHto = '$fileHto',
    methods = ifelse('$methods' == 'NULL', "combined_bff", '$methods'),
    methodsForConsensus = ifelse('$methodsForConsensus' == 'NULL', NULL, '$methodsForConsensus'),
    cellbarcodeWhitelist = ifelse('$cellbarcodeWhitelist' == 'NULL', NULL, '$cellbarcodeWhitelist'),
    metricsFile =ifelse('$metricsFile' == 'NULL', "metrics_cell_hash_r.csv", '$metricsFile'),
    doTSNE = ifelse('$doTSNE' == 'NULL', TRUE, '$doTSNE'),
    preprocess_bff = ifelse('$preprocess_bff' == "FALSE", FALSE, TRUE),
    barcodeWhitelist = ifelse('$barcodeWhitelist' == 'NULL', NULL, '$barcodeWhitelist'),
    doHeatmap = ifelse('$doHeatmap' == 'NULL', TRUE, '$doHeatmap'),
    perCellSaturation = ifelse('$perCellSaturation' == 'NULL', NULL, '$perCellSaturation'),
    majorityConsensusThreshold = ifelse('$majorityConsensusThreshold' == 'NULL', NULL, '$majorityConsensusThreshold'),
    chemistry = ifelse('$chemistry' == 'NULL', "10xV3", '$chemistry'),
    callerDisagreementThreshold = ifelse('$callerDisagreementThreshold' == 'NULL', NULL, '$callerDisagreementThreshold'),
    assignmentOutBff = ifelse('$assignmentOutBff' == 'NULL', 'bff', '$assignmentOutBff'), 
    outputdir = 'bff_${meta.id}'
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)) {
    if (! ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])) {
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided
is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

required_opts <- c('fileHto')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}


#Parameters originally Null
methodsForConsensus <- opt\$methodsForConsensus
if(is.null(methodsForConsensus)){
  methodsForConsensus <- "NULL"
}
cellbarcodeWhitelist <- opt\$cellbarcodeWhitelist
if(is.null(cellbarcodeWhitelist)){
  cellbarcodeWhitelist <- "NULL"
}
perCellSaturation <- opt\$perCellSaturation
if(is.null(perCellSaturation)){
  perCellSaturation <- "NULL"
}
majorityConsensusThreshold <- opt\$majorityConsensusThreshold
if(is.null(majorityConsensusThreshold)){
  majorityConsensusThreshold <- "NULL"
}
callerDisagreementThreshold <- opt\$callerDisagreementThreshold
if(is.null(callerDisagreementThreshold)){
  callerDisagreementThreshold <- "NULL"
}
print(opt\$fileHto)
#Transform logical parameters
do_TSNE <- as.logical(opt\$doTSNE)
do_Heatmap <- as.logical(opt\$doHeatmap)

#saving parameters in a dataframe
Argument <- c("HTO-File", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(opt\$fileHto, opt\$methods, methodsForConsensus, cellbarcodeWhitelist, opt\$metricsFile, perCellSaturation, majorityConsensusThreshold, callerDisagreementThreshold, opt\$doTSNE, opt\$doHeatmap,opt\$chemistry)
params <- data.frame(Argument, Value)

if(opt\$preprocess_bff == TRUE){
  print("Preprocessing activated")
  print(opt\$preprocess)
  #get barcodes
  string <- opt\$barcodeWhitelist
  #separate the barcodes by comma
  words <- strsplit(string, ",")[[1]]
  #Remove leading/trailing whitespace from each word
  words <- trimws(words)
  # Step 3: Create a vector from the barcodesl
  vector <- unlist(words)
  print("Preprocessing")
  counts <- ProcessCountMatrix(rawCountData = opt\$fileHto, barcodeBlacklist = vector)
  print("Preprocessing done")
}else{
  print("No preprocessing")
  counts <- Read10X(opt\$fileHto) 
}

substring_vector <- NULL
if (!is.null(opt\$methodsForConsensus)) { 
  consensus_methods = opt\$methodsForConsensus
  substring_vector <- strsplit(consensus_methods, ",")[[1]]
}

perCell_args <- opt\$perCellSaturation
perCell <- ifelse(perCell_args == "null" || perCell_args == "Null", NULL, perCell_args)

if(opt\$methodsForConsensus=="bff_raw" || opt\$methodsForConsensus=="bff_cluster" || opt\$methodsForConsensus=="bff_raw,bff_cluster" || opt\$methodsForConsensus=="bff_cluster,bff_raw"|| is.null(opt\$methodsForConsensus)  )
  #Only Bff in its different variations is available
  if (opt\$methods == "bff_raw") {
    print("Executing BFF raw")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw"), doTSNE = do_TSNE, doHeatmap = do_Heatmap, methodsForConsensus = substring_vector,cellbarcodeWhitelist = opt\$cellbarcodeWhitelist, metricsFile = opt\$metricsFile, perCellSaturation = opt\$perCellSaturation, majorityConsensusThreshold = opt\$majorityConsensusThreshold, chemistry = opt\$chemistry, callerDisagreementThreshold = opt\$callerDisagreementThreshold )
  }else if (opt\$methods == "bff_cluster") {
    print("Executing BFF cluster")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_cluster"), doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector,cellbarcodeWhitelist = opt\$cellbarcodeWhitelist,metricsFile = opt\$metricsFile, perCellSaturation = opt\$perCellSaturation, majorityConsensusThreshold = opt\$majorityConsensusThreshold, chemistry = opt\$chemistry, callerDisagreementThreshold = opt\$callerDisagreementThreshold)
  }else if (opt\$methods == "combined_bff") {
    print("Executing BFF combined")
    cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw", "bff_cluster") , doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector, cellbarcodeWhitelist = opt\$cellbarcodeWhitelist ,metricsFile = opt\$metricsFile, perCellSaturation = NULL, majorityConsensusThreshold = opt\$majorityConsensusThreshold )
    #cell_hash_R_res <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = c("bff_raw", "bff_cluster") , doTSNE = do_TSNE, doHeatmap = do_Heatmap,methodsForConsensus = substring_vector,cellbarcodeWhitelist = opt\$cellbarcodeWhitelist, metricsFile = opt\$metricsFile, perCellSaturation = NULL, majorityConsensusThreshold = opt\$majorityConsensusThreshold, chemistry = opt\$chemistry, callerDisagreementThreshold = opt\$callerDisagreementThreshold )
  
  }else {
    print("Method not available on the pipeline")
}else{
  print("Consensus only available using BFF methods on the pipeline")
}

if(is.null(cell_hash_R_res)){
  print("No results found")
  df <- data.frame()
  write.csv(df, paste0(opt\$outputdir, "/", opt\$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)
}else{
  write.csv(cell_hash_R_res, paste0(opt\$outputdir, "/", opt\$assignmentOutBff, "_assignment_bff.csv"), row.names=FALSE)
}
write.csv(params, paste0(opt\$outputdir, "/params.csv"))