#' Find diagnostic markers and haplotypes for known aboveground samples
#'
#' @description
#' Finds markers which distinguish one species in the dataset from all other species (species-diagnostic markers) and haplotypes within these which distinguish one individual of the species from all other individuals (individual-diagnostic haplotypes).
#' 
#' @details 
#' 
#' This function takes the output of read.stacks and identifies diagnostic markers and haplotypes which can be used to identify species and individuals in mixed unknown samples. It identifies two classes of diagnostic markers: 1. Species-diagnostic markers, these are markers that are present in at least the proportion of individuals of a species set by the inverse of max.md.marker but not in any individuals outside the species. 2. Individual-diagnostic haplotypes, which are haplotypes which are only present in one individual, although the marker must be present (with a different haplotype) in other individuals of the species. 
#' 
#' @param ind.list a named list with a vector for each species containing the names of all individuals of that species
#' @param stacks.dat output of read.stacks
#' @param max.md.marker the maximum proportion of missing data among individuals of the focal species to call a species-diagnostic marker
#' @param max.md.hap the maximum proportion of missing data in other individuals of a species to call an individual-diagnostic haplotype
#' @param min.dep the minimum read depth required for a marker to be considered present in an individual
#' @param max.haps the maximum number of haplotypes allowed per-marker (no limit if NA)
#' @param verbose whether to output progress information
#' 
#' @return Returns a list, intended as input for match.diag, with the following elements:
#' \itemize{
#' \item{species.diag.markers:}{ a list containing a vector of diagnostic marker index numbers for each species}
#' \item{individual.diag.haplotypes:}{ a list of lists, one for each species, containing a vector of diagnostic haplotype index numbers for each individual of the species, and one of non-individual-diagnostic haplotypes (non.diagnostic)}
#' \item{summary:}{ a list containing dataframes with counts of diagnostic markers and haplotypes for each species and individual}
#' }
#' @importFrom fastmatch fmatch
#' @export
#' 
find.diag <- function(ind.list, stacks.dat, max.md.marker = 0.2, max.md.hap = 0.2, min.dep = 20, max.haps=2, verbose = FALSE){
  if(verbose) message("Processing markers")
  # get individual names
  inds <- gsub("^depth_", "", grep("^depth_", colnames(stacks.dat$by.marker), value = T))
  # check individual names in ind.list and stacks.dat match
  if(any(sort(unlist(ind.list)) != sort(inds))){
    stop("Individual names in ind.list and stacks.dat must match exactly")
  }
  # get names of species with multiple individuals
  multi.ind <- names(which(lengths(ind.list) > 1))
  # remove markers with > max.haps haplotypes which have > min.dep/max.haps depth if max.haps is set
  if(!is.na(max.haps)){
    haps_gt_min_dep <- stacks.dat$by.haplotype[which(do.call("pmax",stacks.dat$by.haplotype[,grep("^depth_",colnames(stacks.dat$by.haplotype),value=T)]) > min.dep/max.haps),"marker"]
    nhaps <- table(haps_gt_min_dep)
    blacklist <- names(nhaps)[which(nhaps>max.haps)]
    stacks.dat$by.haplotype <- stacks.dat$by.haplotype[which(fastmatch::fmatch(stacks.dat$by.haplotype$marker, blacklist,nomatch = 0)==0),]
    stacks.dat$by.marker <- stacks.dat$by.marker[which(fastmatch::fmatch(stacks.dat$by.marker$marker, blacklist,nomatch = 0)==0),]
  } else {
    # if max.haps isn't set, set it to 2 so we can use it to get haplotype depth later
    max.haps <- 2
  }
  # initiate marker presence lists with and without minimum depth thresholds 
  mar.deep <- list()
  mar.shal <- list()
  mar.tables <- list()
  # fill marker presence lists
  for(ind in inds){
    mar.deep[[ind]] <- stacks.dat$by.marker[which(stacks.dat$by.marker[,paste("depth",ind,sep="_")] >= min.dep),"marker"]
    mar.shal[[ind]] <- stacks.dat$by.marker[which(stacks.dat$by.marker[,paste("depth",ind,sep="_")] > 0),"marker"]
  }
  # find species-diagnostic markers
  if(verbose) message("Identifying species-diagnostic markers")
  diag.mar <- list()
  for(sp in names(ind.list)){
    # get vector of all the markers inside and outside the ingroup (i.e. species), keep duplicates for the ingroup so they can be counted
    in.mar <- Reduce(c,mar.deep[ind.list[[sp]]])
    notin.mar <- unique(Reduce(c,mar.shal[setdiff(inds,ind.list[[sp]])]))
    # find markers that are unique to the ingroup
    uniq.mar <- setdiff(in.mar,notin.mar)
    # count all unique markers
    mar.tables[[sp]] <- table(in.mar)[as.character(uniq.mar)]
    # convert counts to proportion of individuals in ingroup and filter by marker.prop
    mar.tables[[sp]] <- mar.tables[[sp]]/length(ind.list[[sp]])
    my.diag <- as.numeric(names(mar.tables[[sp]][which(mar.tables[[sp]]>=1-max.md.marker)]))
    if(verbose) message(paste("     ",length(my.diag),"diagnostic markers for",sp))
    diag.mar[[sp]] <- my.diag
  }
  # make summary list
  diag.sum <- list()
  diag.sum[["species.diagnostic.markers"]] <-  data.frame(species=names(diag.mar),N.diag.markers=lengths(diag.mar),row.names = NULL)
  diag.sum[["individual.diagnostic.haplotypes"]] <- list()
  # find individual-diagnostic haplotypes among species specific markers
  if(verbose) message("Identifying individual-diagnostic haplotypes")
  diag.hap <- list()
  for(sp in multi.ind){
    diag.hap[[sp]] <- list()
    # get dataset for species
    sp.haps <- stacks.dat$by.haplotype[which(fastmatch::fmatch(stacks.dat$by.haplotype$marker,diag.mar[[sp]],nomatch = 0)>0),c("marker","haplotype",paste("depth",ind.list[[sp]],sep="_"),"LOCHAP")]
    # make haplotype missing data list and haplotype presence lists with and without minimum depth thresholds
    hap.deep <- list()
    hap.shal <- list()
    hap.miss <- list()
    for(ind in ind.list[[sp]]){
      hap.deep[[ind]] <- sp.haps[which(sp.haps[,paste("depth",ind,sep="_")] >= min.dep/max.haps),"LOCHAP"]
      hap.shal[[ind]] <- sp.haps[which(sp.haps[,paste("depth",ind,sep="_")] > 0),"LOCHAP"]
      hap.miss[[ind]] <- sp.haps[which(fastmatch::fmatch(sp.haps$marker,as.numeric(names(mar.tables[[sp]][which(mar.tables[[sp]]<1-max.md.hap)])),nomatch = 0)>0) ,"LOCHAP"]
    }
    # convert hap.miss to single list of haplotypes with missing data for the species
    hap.miss <- unique(Reduce(c,hap.miss))
    # get diag for each ind
    for(ind in ind.list[[sp]]){
      # get vector of all the haplotypes with depth over min.dep/max.haps in the individual
      in.hap <- hap.deep[[ind]]
      # get vector of all the haplotypes of any depth in other individuals of the species
      notin.hap <- unique(Reduce(c,hap.shal[ind.list[[sp]][which(ind.list[[sp]]!=ind)]]))
      # find haplotypes which are unique to the ingroup
      uniq.hap <- setdiff(in.hap,notin.hap)
      # remove haplotypes which have missing data in the ingroup or outgroup
      uniq.hap <- setdiff(uniq.hap,hap.miss)
      if(verbose) message(paste("     ",length(uniq.hap)," diagnostic haplotypes for ",sp,": ",ind,sep=""))
      # add to list
      diag.hap[[sp]][[ind]] <- uniq.hap
    }
    # add vector of haplotypes from species-diagnostic markers which aren't diagnostic for any individuals
    diag.hap[[sp]][["non.diagnostic"]] <- setdiff(sp.haps$LOCHAP,unique(Reduce(c,diag.hap[[sp]])))
    diag.sum[["individual.diagnostic.haplotypes"]][[sp]] <- data.frame(individual=names(diag.hap[[sp]]),N.diag.haplotypes=lengths(diag.hap[[sp]]),row.names = NULL)
  }
  out <- list(species.diag.markers = diag.mar, individual.diag.haplotypes = diag.hap, summary = diag.sum)
  if(verbose) message("Done!")
  # return output
  return(out)
} # function end

