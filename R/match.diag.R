#' Match unknown root samples to diagnostic markers and haplotypes 
#'
#' @description
#' Matches markers from unknown belowground samples to diagnostic markers from aboveground samples produced by find.diag  
#' 
#' @details 
#' 
#' Takes the diagnostic markers and haplotypes identified by find.diag and a directory with matches.tsv.gz files from STACKS for all unknown root samples. It then matches the root samples to the diagnostic markers and haplotypes. Minimum depth filters can be set with min.reads.mar and min.reads.hap.
#' 
#' @param data.dir the path to the directory which contains the stacks 'matches.tsv.gz' files for the unknown root samples
#' @param diag list containing diagnostic markers made by find.diag
#' @param min.reads.mar minimum read depth for species-diagnostic marker matches
#' @param min.reads.hap minimum read depth for individual-diagnostic haplotype matches
#' 
#' @return Returns a list with the following elements:
#' \itemize{
#' \item{species.diag.markers:}{ a data frame with one row for each root sample and a column for each species containing counts of species-diagnostic markers detected in each root sample}
#' \item{individual.diag.haplotypes:}{ a list of data frames, one for each species, with one row for each root sample and a column for each individual, containing the proportion of detected species-diagnostic markers for the focal species which had diagnostic haplotypes for each individal}
#' \item{N.reads}{ a data frame containing the number of catalogue-matching reads in each sample}
#' }.
#' @importFrom fastmatch fmatch
#' @export
#' 

match.diag <- function(data.dir, diag, min.reads.mar = 1, min.reads.hap = 1){
  # add a final file separator to data.dir if it doesn't have one
  data.dir <- fix.dir.name(data.dir)
  # get names of samples, species and species with multiple individuals
  samples <- get.samp.names(data.dir)
  species <- names(diag$species.diag.markers)
  multi.ind <- names(diag$individual.diag.haplotypes)
  # make the output list
  out.list <- list()
  out.list[["species.diag.markers"]] <- data.frame(matrix(NA,ncol=length(species)+1,nrow=length(samples),dimnames = list(NULL,c("sample",paste(species,"N.reads",sep=".")))))
  out.list[["species.diag.markers"]]$sample <- samples
  out.list[["individual.diag.haplotypes"]] <- list()
  out.list[["N.reads"]] <- data.frame(matrix(NA,ncol=2,nrow=length(samples),dimnames = list(NULL,c("sample","N.reads"))))
  out.list[["N.reads"]]$sample <- samples
  #for each species with multiple individuals, add a dataframe dataframe to the list
  for(sp in multi.ind){
    out.list[["individual.diag.haplotypes"]][[sp]] <- data.frame(matrix(NA,ncol=length(diag[["individual.diag.haplotypes"]][[sp]])+1,nrow=length(samples),dimnames = list(NULL,c("sample",paste(names(diag[["individual.diag.haplotypes"]][[sp]]),"N.reads",sep=".")))))
    out.list[["individual.diag.haplotypes"]][[sp]]$sample <- samples
  }
  # do by sample so each .gz file only needs to be opened once
  for(samp in 1:length(samples)){
    dat <- read.tsv.gz(data.dir = data.dir, ind = samples[[samp]])
    # make LOCHAP column
    dat$LOCHAP <- paste(dat$marker,dat$haplotype,sep=":")
    # match markers for single-individual species
    for(sp in setdiff(species, multi.ind)){
      out.list[["species.diag.markers"]][samp, paste(sp, "N.reads", sep = ".")] <- sum(dat[which(fastmatch::fmatch(dat$marker, diag[["species.diag.markers"]][[sp]], nomatch = 0) > 0), "depth"])
    }
    # match haplotypes for individuals
    for(sp in multi.ind){
      my.inds <- names(diag[["individual.diag.haplotypes"]][[sp]])
      for(ind in 1:length(my.inds)){
        # get count
        my.hapcount <- sum(dat[which(fastmatch::fmatch(dat$LOCHAP, diag[["individual.diag.haplotypes"]][[sp]][[my.inds[[ind]]]], nomatch = 0) > 0), "depth"])
        # add to output
        out.list[["individual.diag.haplotypes"]][[sp]][samp, paste(my.inds[[ind]], "N.reads", sep = ".")] <- my.hapcount
      }
    }
    # get N reads per root sample
    out.list[["N.reads"]][samp, "N.reads"] <- sum(dat$depth)
  }
  # add marker read numbers for multi-individual species
  for(sp in multi.ind){
    out.list[["species.diag.markers"]][, paste(sp, "N.reads", sep = ".")] <- rowSums(out.list[["individual.diag.haplotypes"]][[sp]][, 2:ncol(out.list[["individual.diag.haplotypes"]][[sp]])])
  }
  # filter sp markers and haplotypes by min.reads.mar and min.reads.hap
  for(sp in species){
    # sp markers
    out.list[["species.diag.markers"]][which(out.list[["species.diag.markers"]][, paste(sp, "N.reads",sep = ".")] < min.reads.mar), paste(sp, "N.reads", sep = ".")] <- 0
    # ind haplotypes
    if(sp %in% multi.ind){
      my.inds <- names(diag[["individual.diag.haplotypes"]][[sp]])
      for(ind in 1:length(my.inds)){
        # set to zero if less than min.reads.hap 
        out.list[["individual.diag.haplotypes"]][[sp]][which(out.list[["individual.diag.haplotypes"]][[sp]][, paste(my.inds[[ind]], "N.reads", sep = ".")] < min.reads.hap), paste(my.inds[[ind]], "N.reads", sep = ".")] <- 0
      }
    }
  }
  # convert haplotype values to proportion of species reads present
  for(sp in multi.ind){
    for(i in 2:ncol(out.list[["individual.diag.haplotypes"]][[sp]])){
      out.list[["individual.diag.haplotypes"]][[sp]][which(out.list[["individual.diag.haplotypes"]][[sp]][, i] > 0), i] <- out.list[["individual.diag.haplotypes"]][[sp]][which(out.list[["individual.diag.haplotypes"]][[sp]][, i] > 0), i] / out.list[["species.diag.markers"]][which(out.list[["individual.diag.haplotypes"]][[sp]][, i] > 0), paste(sp, "N.reads", sep = ".")]
    }
  }
  out.list
}