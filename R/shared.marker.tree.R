#' Make dendrogram from shared markers
#'
#' @description
#' Takes the output of read.stacks and produces a rooted UPGMA tree based on shared marker presence.
#' 
#' @details 
#' 
#' For very closely related species, it is possible that marker presence/absence is not a reliable way to distinguish them. This function produces an Unweighted Pair Group Method with Arithmetic Mean dendrogram of all individuals based on shared marker presence, which can be plotted with e.g ape::plot.phylo. If all conspecific individuals cluster monophyletically, then it suggests that marker presence-absence does distinguish the species in the dataset. This is not intended to produce a phylogenetic tree!
#' 
#' @param stacks.dat output of read.stacks
#' @param ind.list a named list with a vector for each species containing the names of all individuals of that species
#' @param root the name of the outgroup species (or a vector of species names of the outgoup clade) by which to root the dendrogram
#' @param collapse.to.species if TRUE, conspecific individuals are collapsed to produce a "species tree" 
#' 
#' @return an object of class phylo
#' @importFrom phangorn upgma
#' @importFrom ape root getMRCA
#' @export
#' 
#' 
shared.marker.tree <- function(stacks.dat, ind.list, root = NULL, collapse.to.species = F){
  # This function takes the output of read.stacks and produces a rooted UPGMA tree based on shared marker presence. 
  # extract presence absence markers
  my.dat <- list()
  if(collapse.to.species){
    for(sp in names(ind.list)){
      my.sp.list <- list()
      for(ind in ind.list[[sp]]){
        my.sp.list[[ind]] <- stacks.dat[["by.marker"]][which(stacks.dat[["by.marker"]][, paste("depth", ind, sep = "_")] > 0), "marker"]
      }
      my.dat[[sp]] <-  unique(Reduce(c, my.sp.list))
    }
  }
  else{
    for(sp in names(ind.list)){
      for(ind in ind.list[[sp]]){
        my.dat[[paste(sp, ind, sep = "_")]] <- stacks.dat[["by.marker"]][which(stacks.dat[["by.marker"]][, paste("depth", ind, sep = "_")] > 0), "marker"]
      }
    }
  }
  # initialise matrix
  my.mat <- matrix(nrow = length(my.dat), ncol = length(my.dat))
  # fill matrix with the inverse of the proportion of markers shared between each pair of individuals (i.e. the proportion of the markers in the individual with fewest markers which are shared with the individual with more markers).
  for (i in 1:length(my.dat)){
    for (j in 1:length(my.dat)){
      my.mat[i, j] <- 1 - (length(intersect(my.dat[[i]], my.dat[[j]])) / min(length(my.dat[[i]]), length(my.dat[[j]])))
    }
  }
  #name matrix
  rownames(my.mat) <- names(my.dat)
  colnames(my.mat) <- names(my.dat)
  # make distance matrix
  my.mat <- as.dist(my.mat)
  # make tree
  my.tree <- phangorn::upgma(my.mat)
  # root tree if root is set
  if(!is.null(root)){
    if(collapse.to.species){
      my.tree <- ape::root(my.tree, root, resolve.root = T)
    } else {
      if(length(ind.list[[root]]) > 1){
        my.tree <- ape::root(my.tree, ape::getMRCA(my.tree, paste(root, ind.list[[root]], sep = "_")), resolve.root = T)
      } else {
        my.tree <- ape::root(my.tree, paste(root, ind.list[[root]], sep = "_"), resolve.root = T)
      }
    }
  }
  return(my.tree)
}