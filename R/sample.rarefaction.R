#' Run rarefaction analysis on root samples
#'
#' @description
#' 
#' Subsample root data and calculate the number of markers across different subsample sizes to estimate the sufficiency of sequencing. 
#' 
#' @details 
#' 
#' This function attempts to estimate whether the sequencing depth for each root sample is sufficient to recover all the diversity in the sample at the level of species-diagnostic markers, individual-diagnostic haplotypes, number of species and number of individuals. It uses a rarefaction-based approach similar to those used in metabarcoding studies. It randomly subsamples the data n.reps times across a range of subsample sizes (determined by n.subsamples). It then calculates the number of diagnostic haplotypes, markers, species and individuals in each subsample. The mean of these can then be plotted against subsample size. Steeper final curves indicate lower sequencing sufficiency whereas if the curve goes flat, this indicates that most diversity is likely to have been sampled.
#' 
#' @param data.dir he path to the directory which contains the stacks 'matches.tsv.gz' files for the unknown root samples
#' @param diag output of find.diag
#' @param n.subsamples number of subsamples 
#' @param n.reps number of random replicates per subsample
#' @param verbose whether to output progress information
#' 
#'@return Returns a list of data frames. n.reads is a data frame with one subsample and a column for each root sample containing the size of each subsample. Species, individual, marker and haplotype are lists, each containing dataframes with the same dimensions as n.reads with the mean, standard deviation (sd), upper (upper95CI) and lower (lower95CI) 95% confidence intervals.
#' @importFrom fastmatch fmatch
#' @importFrom stats as.dist quantile sd
#' @export
#' 
sample.rarefaction <- function(data.dir, diag, n.subsamples = 10, n.reps = 10, verbose = TRUE){
  # add a final file separator to data.dir if it doesn't have one
  data.dir <- fix.dir.name(data.dir)
  # get names of samples, species and species with multiple individuals
  samples <- get.samp.names(data.dir)
  species <- names(diag$species.diag.markers)
  multi.ind <- names(diag$individual.diag.haplotypes)
  # make the output list
  out.list <- list()
  out.list[["n.reads"]] <- data.frame(matrix(NA, nrow = n.subsamples, ncol = length(samples), dimnames = list(NULL, samples)))
  out.list[["species"]] <- list() 
  out.list[["species"]][["mean"]] <- out.list[["n.reads"]]
  out.list[["species"]][["sd"]] <- out.list[["n.reads"]]
  out.list[["species"]][["upper95CI"]] <- out.list[["n.reads"]] 
  out.list[["species"]][["lower95CI"]] <- out.list[["n.reads"]] 
  out.list[["individuals"]] <- out.list[["species"]]
  out.list[["haplotypes"]] <- out.list[["species"]]
  out.list[["markers"]] <- out.list[["species"]]
  
  # process each sample
  for(samp in 1:length(samples)){
    #open matches.tsv.gz
    if(verbose) message(paste("Sample:",samples[[samp]]))
    dat <- read.tsv.gz(data.dir = data.dir, ind = samples[[samp]])
    # make LOCHAP column to store combined locus and haplotype
    dat$LOCHAP <- paste(dat$marker, dat$haplotype, sep = ":")
    # get total read number
    total.reads <- sum(dat$depth)
    # get subsample sizes
    subsample.sizes <- round(seq(total.reads / n.subsamples, total.reads, total.reads / n.subsamples))
    # make a long form version of the data where rows with depth > 1 are repeated the number of times indicated in the depth column. This allows sampling without replacement.
    dat.long <- dat[rep(seq_len(nrow(dat)) ,dat$depth),]
    # for each subsample size: 
    for (ss in 1:length(subsample.sizes)){
      # vectors to hold species and individual counts
      my.reps.mar <- rep(NA, n.reps)
      my.reps.hap <- rep(NA, n.reps)
      my.reps.spp <- rep(NA, n.reps)
      my.reps.ind <- rep(NA, n.reps)
      # for each replicate, subsample the reads and get the number of species and individuals in each root sample
      for(N in 1:n.reps){
        # make subsample
        my.ss.dat <- dat
        my.subsamp <- sample(x=dat.long$LOCHAP, size = subsample.sizes[[ss]], replace = F)
        my.ss.dat$depth <- as.vector(table(my.subsamp)[dat$LOCHAP])
        my.ss.dat <- my.ss.dat[which(!is.na(my.ss.dat$depth)),]
        # count N species present
        counter <- 0
        for(sp in species){
          if(length(which(fastmatch::fmatch(my.ss.dat$marker, diag[["species.diag.markers"]][[sp]], nomatch = 0) > 0)) > 0){
            counter <- counter + 1
          }
        }
        my.reps.spp[[N]] <- counter
        # count N individuals present
        counter <- 0
        for(sp in multi.ind){
          for(ind in names(diag$individual.diag.haplotypes[[sp]])){
            if(ind != "non.diagnostic" & length(which(fastmatch::fmatch(my.ss.dat$LOCHAP, diag[["individual.diag.haplotypes"]][[sp]][[ind]], nomatch = 0) > 0)) > 0){
              counter <- counter + 1
            }
          }
        }
        my.reps.ind[[N]] <- counter
        # count markers and haplotypes
        my.reps.mar[[N]] <- length(unique(my.ss.dat$marker))
        my.reps.hap[[N]] <- length(unique(my.ss.dat$LOCHAP))
      }
      #add mean, sd and CIs to output list
      out.list[["species"]][["mean"]][ss, samples[[samp]]] <- mean(my.reps.spp)
      out.list[["species"]][["sd"]][ss, samples[[samp]]] <- stats::sd(my.reps.spp)
      out.list[["species"]][["upper95CI"]][ss, samples[[samp]]] <- quantile(my.reps.spp, 0.975)
      out.list[["species"]][["lower95CI"]][ss, samples[[samp]]] <- quantile(my.reps.spp, 0.025)
      out.list[["individuals"]][["mean"]][ss, samples[[samp]]] <- mean(my.reps.ind)
      out.list[["individuals"]][["sd"]][ss, samples[[samp]]] <- stats::sd(my.reps.ind)
      out.list[["individuals"]][["upper95CI"]][ss, samples[[samp]]] <- quantile(my.reps.ind, 0.975)
      out.list[["individuals"]][["lower95CI"]][ss, samples[[samp]]] <- quantile(my.reps.ind, 0.025)
      out.list[["markers"]] [["mean"]][ss, samples[[samp]]] <- mean(my.reps.mar)
      out.list[["markers"]] [["sd"]][ss, samples[[samp]]] <- stats::sd(my.reps.mar)
      out.list[["markers"]][["upper95CI"]][ss, samples[[samp]]] <- quantile(my.reps.mar, 0.975)
      out.list[["markers"]][["lower95CI"]][ss, samples[[samp]]] <- quantile(my.reps.mar, 0.025)
      out.list[["haplotypes"]][["mean"]][ss, samples[[samp]]] <- mean(my.reps.hap)
      out.list[["haplotypes"]][["sd"]][ss, samples[[samp]]] <- stats::sd(my.reps.hap)
      out.list[["haplotypes"]][["upper95CI"]][ss, samples[[samp]]] <- quantile(my.reps.hap, 0.975)
      out.list[["haplotypes"]][["lower95CI"]][ss, samples[[samp]]] <- quantile(my.reps.hap, 0.025)
      rm(list = c("my.ss.dat", "my.subsamp"))
      gc()
    }
    #add subsample sizes to output list
    out.list[["n.reads"]][, samples[[samp]]] <- subsample.sizes
    rm(list = c("dat", "dat.long"))
    gc()
  }
  out.list
}
