
fix.dir.name <- function(data.dir){
  # if data.dir doesn't end with file separator
  if(substr(data.dir, nchar(data.dir), nchar(data.dir)) != .Platform$file.sep){
    # then add it
    data.dir <- paste(data.dir, .Platform$file.sep, sep="")
  }
  data.dir
}

read.tsv.gz <- function(data.dir, ind, min.dep = 0, coln = c("marker", "haplotype", "depth")){
  # gz file
  zz <- gzfile(paste(data.dir, ind, ".matches.tsv.gz", sep = ""), 'rt')
  # read from gz file
  dat <- utils::read.table(zz, header = F, sep = "\t", comment.char = "#")
  # close gz file
  close(zz)
  # keep only necessary columns
  dat <- dat[, c(1, 4, 5)]
  # name columns
  colnames(dat) <- coln
  # filter out haplotypes with depth lower than min_dep
  dat <- dat[which(dat[, coln[[3]]] >= min.dep), ]
  # return table
  dat
}
  
get.samp.names <- function(data.dir){
  # file names
  files <- list.files(data.dir, "*matches.tsv.gz")
  # remove suffix
  samples <- sub(files, pattern = ".matches.tsv.gz", replacement = "")
  # output
  samples
}

