


#' Read stacks tsv.gz files for RootID
#'
#' @description
#' Read stacks *.tsv.gz files for known aboveground samples and convert them into input for the RootID package.
#' 
#' @details 
#' 
#' This function is the first step in the RootID package, the results is the main input for the find.diag function. All aboveground samples should be in a single directory (without belowground samples) and should all have the "matches.tsv.gz" suffix. If min.dep is set too low, it can use a large amount of memory, the default works fine on a standard laptop with our test data but larger datasets may require more memory or a higher min.dep setting. 
#' 
#' @param data.dir path to the directory containing aboveground tsv.gz files
#' @param min.dep minimum depth to include haplotypes
#' @param verbose whether to output progress information
#' 
#' @return a list containing two data frames, containing the depth of each marker (by.marker) or haplotype (by.haplotype) in each individual. 
#' @importFrom data.table as.data.table merge.data.table
#' @importFrom stats aggregate
#' @export
#' 
#' 
read.stacks <- function(data.dir, min.dep = 3, verbose = F){
  # add a final file separator to data.dir if it doesn't have one
  data.dir <- fix.dir.name(data.dir)
  # get list of individuals
  inds <- get.samp.names(data.dir)
  dat.list <- list()
  for(ind in inds){
    if(verbose) message(paste("Loading ", ind, ".matches.tsv.gz", sep = ""))
    # open and filter matches.tsv.gz
    dat <- read.tsv.gz(data.dir = data.dir, ind = ind, min.dep = min.dep, coln = c("marker", "haplotype", paste("depth", ind, sep = "_")))
    # add to list
    dat.list[[ind]] <- data.table::as.data.table(dat)
  }
  if(verbose) message("Merging data")
  # merge the dataframes for each individual
  df.hap <- Reduce(function(x, y) data.table::merge.data.table(x, y, by = c("marker", "haplotype"), all = T, allow.cartesian = T), dat.list)
  df.hap <- as.data.frame(df.hap)
  df.hap[is.na(df.hap)] <- 0
  # make per-marker dataframe with haplotype-specific depths for each locus summed
  if(verbose) message("Summing haplotype depths per marker")
  df.marker <- df.hap[, c(1, 3:ncol(df.hap))]
  df.marker <- stats::aggregate(. ~ marker, data = df.marker, FUN = sum)
  # make LOCHAP column for df.hap
  df.hap$LOCHAP <- paste(df.hap$marker, df.hap$haplotype, sep = ":")
  # return list of the two dataframes
  out.list <- list(by.marker = df.marker, by.haplotype = df.hap)
  if(verbose) message("Done")
  out.list
} 
