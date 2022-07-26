#' @importFrom rgl subdivision3d cube3d scale3d translate3d shade3d cylinder3d wire3d legend3d par3d mtext3d grid3d spheres3d
#' @importFrom reshape2 melt
#' @importFrom Rfast rowMaxs
#' @importFrom scales rescale

get.root.dat <- function(M, matches, ind.list, tax.type, taxa, scale.type = "colmax"){
  # check scale type
  if(!(scale.type %in% c("max", "colmax"))){
    stop("abundance.scale.type must be either 'max' or 'colmax'")
  }
  if(tax.type == "sp"){
    # for each species in taxa
    for(sp in taxa){
      if(max(matches$species.diag.markers[[paste(sp,".N.reads",sep="")]])>0){
        if(scale.type == "colmax"){
          # extract N species diagnostic markers for each root sample scaled by max markers detected for that species in any root sample
          M[,sp] <- matches$species.diag.markers[, paste(sp, ".N.reads", sep = "")][match(M$value, matches$species.diag.markers$sample)] / max(matches$species.diag.markers[, paste(sp, ".N.reads", sep = "")])
        } else if(scale.type == "max"){
          # extract N species diagnostic markers for each root sample scaled by max markers detected for any species in taxa
          M[,sp] <- matches$species.diag.markers[, paste(sp, ".N.reads", sep = "")][match(M$value, matches$species.diag.markers$sample)]/max(matches$species.diag.markers[, paste(taxa, ".N.reads", sep = "")])
        }
      } else {
        # skip species which weren't detected in any root samples
        warning(paste("Skipping", sp, "as the species wasn't detected in any root samples"))
      }
    }
  } else if(tax.type == "ind"){
    # for each individual
    for(ind in taxa){
      # identify correct species 
      sp <- names(ind.list)[[which(sapply(ind.list, function(x) ind %in% x))]]
      if(max(matches$individual.diag.haplotypes[[sp]][, paste(ind, ".N.reads", sep = "")]) > 0){
        if(scale.type == "colmax"){
          # extract proportion of individual diagnostic haplotypes for each root sample scaled by max proportion of haplotypes detected for that individual in any root sample
          M[,ind] <- matches$individual.diag.haplotypes[[sp]][, paste(ind, ".N.reads", sep = "")][match(M$value, matches$individual.diag.haplotypes[[sp]]$sample)] / max(matches$individual.diag.haplotypes[[sp]][, paste(ind, ".N.reads", sep = "")])
        } else if(scale.type == "max"){
          # extract proportion of individual diagnostic haplotypes for each root sample scaled by max proportion of haplotypes detected for any individual in taxa
          M[, ind] <- matches$individual.diag.haplotypes[[sp]][, paste(ind, ".N.reads", sep = "")][match(M$value, matches$individual.diag.haplotypes[[sp]]$sample)] / max(matches$individual.diag.haplotypes[[sp]][, paste(taxa, ".N.reads", sep = "")])
        }
      } else {
        # skip individuals which weren't detected in any root samples
        warning(paste("Skipping", ind, "as the individual wasn't detected in any root samples"))
      }
    }
  } else if(tax.type == "sp.ind"){
    # for each species in taxa
    for(sp in taxa){
      if(length(ind.list[[sp]])>1){
        # for each individual in the species
        for(ind in ind.list[[sp]]){
          if(max(matches$individual.diag.haplotypes[[sp]][,paste(ind,".N.reads",sep="")]) > 0){
            if(scale.type == "colmax"){
              # extract proportion of individual diagnostic haplotypes for each root sample scaled by max proportion of haplotypes detected for that individual in any root sample
              M[,ind] <- matches$individual.diag.haplotypes[[sp]][,paste(ind,".N.reads",sep="")][match(M$value, matches$individual.diag.haplotypes[[sp]]$sample)]/max(matches$individual.diag.haplotypes[[sp]][,paste(ind,".N.reads",sep="")])
            } else if(scale.type == "max"){
              # extract proportion of individual diagnostic haplotypes for each root sample scaled by max proportion of haplotypes detected for any individual of that species
              M[,ind] <- matches$individual.diag.haplotypes[[sp]][,paste(ind,".N.reads",sep="")][match(M$value, matches$individual.diag.haplotypes[[sp]]$sample)]/max(matches$individual.diag.haplotypes[[sp]][,paste(ind.list[[sp]],".N.reads",sep="")])
            }
          } else {
            # skip individuals which weren't detected in any root samples
            warning(paste("Skipping", ind, "as the individual wasn't detected in any root samples"))
          }
        }
      } else {
        # skip species without multiple individuals
        warning(paste("Skipping", sp, "as it contains only one individual"))
      }
    }
  } else if(tax.type == "vec"){
    if(scale.type == "max"){
      warning("scale.type = 'max' is only implemented for tax.type = 'sp', 'ind' and 'sp.ind'. Using scale.type = 'colmax'")
    }
    for(n in taxa){
      # if the element is a species name
      if(n %in% names(ind.list)){
        sp <- n
        if(max(matches$species.diag.markers[[paste(sp,".N.reads",sep="")]]) > 0){
          M[,sp] <- matches$species.diag.markers[[paste(sp,".N.reads",sep="")]][match(M$value,matches$species.diag.markers$sample)]/max(matches$species.diag.markers[[paste(sp,".N.reads",sep="")]])
        } else {
          # skip species which weren't detected in any root samples
          warning(paste("Skipping", sp, "as the species wasn't detected in any root samples"))
        }
        # if the element is an individual name
      } else {
        ind <- n
        # identify correct species 
        sp <- names(ind.list)[[which(sapply(ind.list,function(x) ind %in% x))]]
        if(max(matches$individual.diag.haplotypes[[sp]][[paste(ind,".N.reads",sep="")]]) > 0){
          M[,ind] <- matches$individual.diag.haplotypes[[sp]][[paste(ind,".N.reads",sep="")]][match(M$value, matches$individual.diag.haplotypes[[sp]]$sample)]/max(matches$individual.diag.haplotypes[[sp]][[paste(ind,".N.reads",sep="")]])
        } else {
          # skip individuals which weren't detected in any root samples
          warning(paste("Skipping", ind, "as the individual wasn't detected in any root samples"))
        }
      }
    }
  } else {
    stop("tax.type must be either 'sp', 'ind', 'sp.ind' or 'vec'")
  }
  M
}

rgl.ellipsoid <- function (x=0, y=0, z=0, a=1, b=1, c=1,subdivide = 5){
  # adapted from cda::rgl.ellipsoid, downloaded from https://github.com/nano-optics/cda/blob/master/R/rgl.r
  sphere <- rgl::subdivision3d(rgl::cube3d(), subdivide)
  class(sphere) <- c("mesh3d","shape3d")
  
  norm <- sqrt(sphere$vb[1, ]^2 + 
                 sphere$vb[2, ]^2 + 
                 sphere$vb[3, ]^2 )
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
  sphere$vb[4, ] <- 1
  sphere$normals <- sphere$vb
  result <- rgl::scale3d(sphere, a,b,c)
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}

tree.cols <- function(ind.list, tax.type, cols){
  if(tax.type %in% c("sp.ind","vec","ind")){
    out <- cols
  } else if(tax.type == "sp"){
    out <- c()
    for(s in names(cols)){
      for(i in ind.list[[s]]){
        out[[i]] <- cols[[s]]
      }
    }
  }
  out
}

tree.3D <- function(tree.dat,col,trunk.col="saddlebrown",z.base=0,trunk.sides=100,trunk.rad=0.1,alpha=1,subdivide=5){
  # check column names
  if(any(sort(colnames(tree.dat)) != c("crown.height","height","rad","sample","x","y"))) {
    stop("tree.dat must be a data frame with columns: 'crown.height','height','rad','sample','x','y'")
  }
  # make colour vectors correct length
  if(length(col) < nrow(tree.dat)){
    col <- rep(col,ceiling(nrow(tree.dat)/length(col)))
  }
  if(length(trunk.rad) < nrow(tree.dat)){
    trunk.rad <- rep(trunk.rad,ceiling(nrow(tree.dat)/length(trunk.rad)))
  }
  if(length(trunk.col) < nrow(tree.dat)){
    trunk.col <- rep(trunk.col,ceiling(nrow(tree.dat)/length(trunk.col)))
  }
  # replace missing values for crown height or height
  tree.dat[which(is.na(tree.dat$height)),"height"] <- 0
  tree.dat[which(is.na(tree.dat$crown.height)),"crown.height"] <- 0
  for(n in 1:nrow(tree.dat)){
    # draw crown
    rgl::shade3d(rgl.ellipsoid(x=tree.dat[n,"x"],
                               y=tree.dat[n,"y"],
                               z=z.base+mean(c(tree.dat[n,"crown.height"],tree.dat[n,"height"])),
                               a=tree.dat[n,"rad"],
                               b=tree.dat[n,"rad"],
                               c=(tree.dat[n,"height"]-tree.dat[n,"crown.height"])/2,
                               subdivide = subdivide),
                 col=col[[n]],
                 alpha=alpha,
                 specular = col[[n]])
    # draw trunk
    if(tree.dat[n,"crown.height"] > 0){
      rgl::shade3d(rgl::cylinder3d(center=cbind(rep(tree.dat[n,"x"],2),rep(tree.dat[n,"y"],2),c(z.base,z.base+tree.dat[n,"crown.height"])),
                                   sides = trunk.sides,
                                   radius = trunk.rad[[n]],
                                   closed = -2),
                   col=trunk.col[[n]],
                   alpha=alpha)
    }
  }
}

root.cubes <- function(M, which.draw = "a",threshold=0,scale=c(1,1,1),which.vec = NA,type="c",cube.col="cornflowerblue",cube.alpha=0.1){
  if(which.draw != "n"){
    if(which.draw == "a"){
      n.v <- 1:nrow(M)
    } else if (which.draw == "p") {
      n.v <- which(M$score > threshold)
    } else if (which.draw == "s") {
      if(!is.na(which.vec)){
        n.v <- which.vec
      } else {
        stop("which.vec (a vector specifying which rows to draw cubes from) must be provided if which.draw is set to 's'")
      }
    } else {
      stop("which.draw must be 'a' (all), 'p' (present), 's' (specified) or 'n' (none)")
    }
    for (n in n.v){
      if(type=="l"){
        rgl::wire3d(rgl::translate3d(rgl::scale3d(rgl::cube3d(),scale[1]/2,scale[2]/2,scale[3]/2),x=M[n,"Var1"]*scale[1],y=M[n,"Var2"]*scale[2],z=M[n,"Var3"]*scale[3]))
      } else if(type=="c"){
        rgl::shade3d(rgl::translate3d(rgl::scale3d(rgl::cube3d(),scale[1]/2,scale[2]/2,scale[3]/2),x=M[n,"Var1"]*scale[1],y=M[n,"Var2"]*scale[2],z=M[n,"Var3"]*scale[3]),col=cube.col,alpha=cube.alpha)
      } else {
        stop("type must be 'l' (lines) or 'c' (cube)")
      }
      
    }
  }
}

root.bbox <- function(M,scale=c(1,1,1),draw=F){
  if(draw){
    col="black" 
  } else {
    col="transparent"
  }
  rgl::wire3d(
    rgl::translate3d(rgl::scale3d(rgl::cube3d(),
                                  scale[1]/2*max(M$Var1),
                                  scale[2]/2*max(M$Var2),
                                  scale[3]/2*max(M$Var3)),
                      x=mean(c(max(M$Var1),min(M$Var1)))*scale[1],
                      y=mean(c(max(M$Var2),min(M$Var2)))*scale[2],
                      z=mean(c(max(M$Var3),min(M$Var3)))*scale[3]),
    color=col,
    meshColor="edges")
}

#' @importFrom stats runif 
root.spheres <- function(M, threshold = 0, col="black",max.n.spheres,scale,radius,alpha){
  for (n in which(M$score > threshold)){
    # scale N particles to root abundance
    n.spheres <- round(M[n,"score"]*max.n.spheres)
    coords <- list()
    # generate particle coordinates
    for(C in 1:3)(
      coords[[C]] <- runif(n.spheres,M[n,paste("Var",C,sep="")]*scale[[C]]-scale[[C]]/2, M[n,paste("Var",C,sep="")]*scale[[C]]+scale[[C]]/2)
    )
    # draw particles
    rgl::spheres3d(x=coords[[1]], y=coords[[2]], z=coords[[3]], radius = radius, color = col, alpha = alpha) 
  }
}

root.alpha.scale.cubes <- function(M, col, scale, min.alpha, max.alpha){
  # check if more than one species
  if(ncol(M) > 5){
    warning("Coloured cube mode only works with one species or individual. Only using the first data column.")
  }
  if(max(M[,5]) > 0){
    # scale to min and max alpha
    M[which(M[,5] > 0),5] <- scales::rescale(M[which(M[,5] > 0),5], to=c(min.alpha, max.alpha))
    #M[,5] <- M[,5]/max(M[,5])*max.alpha
    for(n in which(M[,5] > 0)){
      rgl::shade3d(rgl::translate3d(rgl::scale3d(rgl::cube3d(),scale[1]/2,scale[2]/2,scale[3]/2),x=M[n,"Var1"]*scale[1],y=M[n,"Var2"]*scale[2],z=M[n,"Var3"]*scale[3]),col=col,alpha=M[n,5])
    }
  } else {
    warning("No detected roots")
  }
}

root.grid <- function(M,scale,which.root.grid){
  if("x+" %in% which.root.grid) rgl::grid3d("x+",at=list(y=(1:(max(M$Var2)-1)+0.5)*scale[[2]],z=(1:(max(M$Var3)-1)+0.5)*scale[[3]]),lty =2)
  if("x-" %in% which.root.grid) rgl::grid3d("x-",at=list(y=(1:(max(M$Var2)-1)+0.5)*scale[[2]],z=(1:(max(M$Var3)-1)+0.5)*scale[[3]]),lty =2)
  if("y+" %in% which.root.grid) rgl::grid3d("y+",at=list(x=(1:(max(M$Var1)-1)+0.5)*scale[[1]],z=(1:(max(M$Var3)-1)+0.5)*scale[[3]]),lty =2)
  if("y-" %in% which.root.grid) rgl::grid3d("y-",at=list(x=(1:(max(M$Var1)-1)+0.5)*scale[[1]],z=(1:(max(M$Var3)-1)+0.5)*scale[[3]]),lty =2)
  if("z+" %in% which.root.grid) rgl::grid3d("z+",at=list(x=(1:(max(M$Var1)-1)+0.5)*scale[[1]],y=(1:(max(M$Var2)-1)+0.5)*scale[[2]]),lty =2)
  if("z-" %in% which.root.grid) rgl::grid3d("z-",at=list(x=(1:(max(M$Var1)-1)+0.5)*scale[[1]],y=(1:(max(M$Var2)-1)+0.5)*scale[[2]]),lty =2)
}

root.axes <- function(M, scale, root.x.names=NULL, root.y.names=NULL, root.z.names=NULL,which.axes, which.sides){
  # positions
  x.at <- sort(unique(M$Var1)) * scale[[1]]
  y.at <- sort(unique(M$Var2)) * scale[[2]]
  z.at <- sort(unique(M$Var3)) * scale[[3]]
  # names
  if(is.null(root.x.names)) root.x.names <- 1:length(x.at)
  if(is.null(root.y.names)) root.y.names <- 1:length(y.at)
  if(is.null(root.z.names)) root.z.names <- 1:length(z.at)
  #
  
  # add labels
  if(which.axes[[1]] == TRUE){
    rgl::mtext3d(root.x.names, edge=paste("x",which.sides[[1]],sep = ""), at=x.at, line = 2) # pos=c(min(y.at),min(z.at))
  }
  if(which.axes[[2]] == TRUE){
    rgl::mtext3d(root.y.names, edge=paste("y",which.sides[[2]],sep = ""), at=y.at, line = 2)
  }
  if(which.axes[[3]] == TRUE){
    rgl::mtext3d(rev(root.z.names), edge=paste("z",which.sides[[3]],sep = ""), at=z.at, line = 2, cex = 0.8)
  }
}

#' Plot a 3d root map
#'
#' @description
#' 
#' This function uses the output of `match.diag` and information on the position and dimensions of root samples and trees to show the root sampling layout as a three-dimensional grid. Each grid square represents one root sample, and visually displays the abundance of the focal tree or species (either in the form of colour intensity or density of randomly distributed particles within each root sample). Optional three-dimensional models of the trees show their position, height, crown base height, and crown diameter.
#' 
#' @details 
#' 
#' The root.pos option shows the position of the root samples relative to each other. It is a 3-dimensional array with dimensions (x, y, z), where x and y are the number of samples along the x and y horizontal axes (corresponding to those in tree.dat) and z is the number of samples along the depth axis (with deeper samples coming first). For example, for samples arranged in a 2 x 2 horizontal grid at 2 depths named "x1_y1_d1", "x2_y1_d1" ... etc, the correct root.pos input would be:
#' 
#'  `array(c("x1_y2_d2", "x1_y1_d2", "x2_y2_d2", "x2_y1_d2", "x1_y2_d1", "x1_y1_d1", "x2_y2_d1", "x2_y1_d1"), dim = c(2,2,2))`. 
#'  
#'  The function currently only works in cases where root samples are evenly spread in each direction. 
#' 
#' The taxa and tax.type options control which species or individuals are plotted. For `tax.type = 'sp'`, one or more species' roots are plotted and `taxa` should contain the names of these species. For `tax.type = 'ind'`, one or more individuals' roots are plotted and `taxa` should contain the names of these individuals. For `tax.type = 'sp.ind'`, all individuals of one or more species are plotted and `taxa` should contain the names of these species. For `tax.type = 'vec'` a mixture of individuals and species can be plotted (listed in `taxa`).
#' 
#' For `tax.type = 'sp'`, root presence and density are determined from species-diagnostic markers. For `tax.type = 'sp.ind'` and `tax.type = 'ind'`, root presence and density are determined from individual-diagnostic haplotypes. For `tax.type = 'vec'`, root presence and density are determined by individual-diagnostic haplotypes for individuals and species-diagnostic markers for species.
#' 
#' How the root abundance is scaled can be controlled with `abundance.scale.type`. For `abundance.scale.type = colmax`, abundance (i.e. number of diagnostic marker or haplotype reads detected) is scaled by the maximum detected in any root sample for each species or individual separately. For `abundance.scale.type = max`, abundance is scaled by the maximum detected in any root sample across all species or individuals in the plot.
#' 
#' To draw trees, a data frame - `tree.dat` - must be provided. It should have one row for each individual tree and the following columns: 'sample' (the tree name), 'x', 'y' (the tree position on the x and y horizontal axes), 'rad' (tree canopy radius), 'height' (tree height) and 'crown.height' (the height of the bottom of the canopy). These measurements should be in units of the horizontal distance between root samples (i.e. if root samples are one metre apart on the x and y axes, they should be in metres, but if root samples are 2 metres apart, they should be in metres / 2).
#' 
#' @param matches the output of `match.diag`.
#' @param ind.list a named list with a vector for each species containing the names of all individuals of that species.
#' @param root.pos a 3D array containing the names of each root sample. See details. 
#' @param taxa a character vector containing one or more names of species or individuals in `matches`
#' @param tax.type either 'sp', 'ind', 'sp.ind', or 'vec' (see details).
#' @param cols a named character vector containing colors for each species and/or individual to be plotted.
#' @param draw.legend logical, whether to draw a legend.
#' @param abundance.scale.type either 'colmax' or 'max' (see details)
#' @param draw.particles logical, whether to draw particles to indicate root abundance. 
#' @param max.n.particles integer, the maximum number of root particles to draw per root sample
#' @param particle.radius the radius of root particles
#' @param particle.alpha float between 0 and 1 determining the opacity of root particles (0 = transparent, 1 = opaque).
#' @param particle.threshold minimum abundance above which to draw root particles
#' @param alpha.scale.cubes logical, whether to draw scaled transparent cubes to indicate root abundance.
#' @param max.alpha.scale.cubes float between 0 and 1, determining the maximum alpha (transparency) of root cubes.
#' @param min.alpha.scale.cubes the minimum alpha of root cubes. 
#' @param alpha.scale.cube.col scaled root cube colour
#' @param draw.root.cubes determines whether uniform (unscaled) cubes should be drawn around root samples. Either  'a' (all samples), 'p' (samples where roots are present), 's' (samples specified in `cube.which.vec`) or 'n' (none).
#' @param cube.threshold the minimum abundance above which to draw uniform cubes if `draw.root.cubes = 'p'`
#' @param cube.which.vec vector of root sample names to draw uniform cubes round if `draw.root.cubes = 's'`.
#' @param cube.type determines what type of uniform cubes should be drawn. Either 'l' (lines, i.e. a wire frame) or 'c' (cube, a semi-transparent cube).
#' @param cube.alpha float between 0 and 1 determining the opacity of uniform root cubes
#' @param cube.col uniform root cube colour
#' @param draw.root.bbox logical, whether to draw a bounding box around the root samples
#' @param draw.root.grid logical, whether to draw grids on the edge of the root area for x, y and depth axes
#' @param which.root.grid character vector determining which sides the grids should be on (see `?rgl::grid3d` for details)
#' @param root.scale a numeric vector of length = 3, determining how the root axes should be scaled relative to each other, e.g. should the depth axis be flattened. 
#' @param draw.trees logical, whether 3D scale models of the trees should be drawn.
#' @param tree.dat a data frame with (scaled) tree dimension info (see details).
#' @param tree.trunk.col colour of the tree trunks.
#' @param tree.trunk.rad tree trunk radius (currently the same for all trees in the plot).
#' @param tree.alpha transparency of the trees 
#' @param tree.trunk.sides the number of sides of the tree trunk polygon (higher = smoother but more computationally intensive).
#' @param tree.crown.subdivide the number of subdivisions of the tree crown (higher = smoother but more computationally intensive).
#' @param draw.axes logical, whether to draw axis labels
#' @param root.x.names,root.y.names,root.z.names names for the x, y and z (depth) axes.
#' @param which.axes logical vector of length = 3, which axes (x, y and z) to plot axis labels for.
#' @param which.axis.sides vector of length = 3 of '+', '-', '++' or '--', determining which side of the axes (x, y and z), axis labels should be plotted on (see 'edge' option in `?rgl::mtext3d` for details)
#' 
#' @return returns NULL and plots an `rgl` object
#' @importFrom grDevices rainbow 
#' @export
#' 
plot_roots_3d <- function(# essential options
                          matches, 
                          ind.list, 
                          root.pos, 
                          taxa, 
                          tax.type,
                          # general plotting options
                          cols=NULL,
                          draw.legend = NULL,
                          abundance.scale.type = "colmax",
                          # root particles
                          draw.particles = TRUE,
                          max.n.particles = 500,
                          particle.radius = 0.02,
                          particle.alpha = 1,
                          particle.threshold = 0,
                          # transparent scaled root cubes
                          alpha.scale.cubes = FALSE,
                          max.alpha.scale.cubes = 0.5,
                          min.alpha.scale.cubes = 0.05,
                          alpha.scale.cube.col = "grey80",
                          # uniform root cubes
                          draw.root.cubes = "n",
                          cube.threshold = 0,
                          cube.which.vec = NA,
                          cube.type = "c",
                          cube.alpha = 0.1,
                          cube.col = "cornflowerblue",
                          # additional root decorations
                          draw.root.bbox = TRUE,
                          draw.root.grid = FALSE,
                          which.root.grid = c("x-", "y+", "z-"),
                          root.scale = c(1,1,0.2),
                          # tree drawing options
                          draw.trees = TRUE, 
                          tree.dat = NULL,
                          tree.trunk.col = "saddlebrown",
                          tree.trunk.rad = 0.1,
                          tree.alpha=0.5,
                          tree.trunk.sides = 100,
                          tree.crown.subdivide = 5,
                          # axis options
                          draw.axes = TRUE,
                          root.x.names=NULL, 
                          root.y.names=NULL, 
                          root.z.names=NULL,
                          which.axes = c(TRUE,TRUE,TRUE),
                          which.axis.sides = c("-","+","++")
                          ){
  if(is.null(tree.dat) & draw.trees == TRUE){
    warning("To draw trees, tree.dat must be provided. Setting draw.trees to FALSE")
    draw.trees <- FALSE
  }
  ####### ADD REST OF OPTION CHECKS
  on.exit(rgl::par3d(skipRedraw=FALSE))
  # initialise data frame with root position
  M <- reshape2::melt(root.pos)
  # add root data to data frame
  M <- get.root.dat(M = M, matches = matches, ind.list = ind.list, tax.type = tax.type, taxa = taxa, scale.type = abundance.scale.type)
  # set up colours
  if(is.null(cols)){
    cols <- rainbow(ncol(M)-4)
    names(cols) <- colnames(M)[5:ncol(M)]
  }
  if(any(sort(names(cols)) != sort(colnames(M)[5:ncol(M)]))){
    stop(paste("names of cols must match the selected species or individuals:", paste(colnames(M)[5:ncol(M)],collapse = ", ")))
  }
  root.cols <- cols
  tree.cols <- tree.cols(ind.list = ind.list, tax.type = tax.type, cols = cols)
  # prevent elements being drawn until the function ends, this speeds up drawing, especially for complex plots.
  rgl::par3d(skipRedraw=TRUE)
  # 3D plots
  # draw bbox
  if(draw.root.bbox){
    root.bbox(M = M, scale = root.scale,draw = TRUE)
  } else {
    root.bbox(M = M, scale = root.scale,draw = FALSE)
  }
  # draw legend
  if(is.null(draw.legend)){
    if(ncol(M) > 5) draw.legend <- TRUE else draw.legend <- FALSE
  }
  if(draw.legend){
    rgl::legend3d("top", fill=root.cols, legend = names(root.cols))
  }
  # draw grid
  if(draw.root.grid){
    root.grid(M = M, scale = root.scale, which.root.grid = which.root.grid)
  }
  # draw axes
  if(draw.axes){
    root.axes(M = M, scale = root.scale, root.x.names = root.x.names, root.y.names = root.y.names, root.z.names = root.z.names, which.sides = which.axis.sides, which.axes = which.axes)
  }
  # draw root cubes
  if(draw.root.cubes != "n"){
    my.M <- M[,1:4]
    my.M$score <- Rfast::rowMaxs(as.matrix(M[,5:ncol(M)]),value=T)
    root.cubes(M = my.M, which.draw = draw.root.cubes, scale = root.scale, type = cube.type, cube.col = cube.col, cube.alpha = cube.alpha, which.vec =  cube.which.vec, threshold = cube.threshold)
  }
  # draw transparent scaled root cubes
  if(alpha.scale.cubes){
    if(tax.type == "sp.ind"){
      Mcubes <- get.root.dat(M = reshape2::melt(root.pos), matches = matches, ind.list = ind.list, tax.type = "sp", taxa = taxa, scale.type = "colmax")
      root.alpha.scale.cubes(M = Mcubes, col = alpha.scale.cube.col, scale = root.scale, max.alpha = max.alpha.scale.cubes, min.alpha = min.alpha.scale.cubes)
    } else {
      alpha.scale.cube.col <- root.cols[[1]]
      root.alpha.scale.cubes(M = M[,1:5], col = alpha.scale.cube.col, scale = root.scale, max.alpha = max.alpha.scale.cubes, min.alpha = min.alpha.scale.cubes)
    }
  }
  # draw spheres
  if(draw.particles){
    for(n in 5:ncol(M)){
      my.M <- M[,c(1:4,n)]
      colnames(my.M) <- c("Var1","Var2","Var3","value","score")
      root.spheres(M = my.M, col = root.cols[[colnames(M)[[n]]]], max.n.spheres = max.n.particles, radius = particle.radius, scale = root.scale, threshold = particle.threshold, alpha = particle.alpha)
    }
  }
  # draw trees
  if(draw.trees){
    z.base <- max(M$Var3)*root.scale[[3]]+root.scale[[3]]/2
    tree.dat <- tree.dat[which(tree.dat$sample %in% names(tree.cols)),]
    tree.3D(tree.dat = tree.dat,col=tree.cols[tree.dat$sample],z.base = z.base, trunk.rad = tree.trunk.rad, trunk.sides = tree.trunk.sides, trunk.col = tree.trunk.col, alpha = tree.alpha, subdivide = tree.crown.subdivide)
  }
}
