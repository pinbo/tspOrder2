library(qtl)

#' use LKH2 to reorder the markers within a R/qtl object
#' @usage tspOrde2(cross, method="concorde", execPath, chr=chrnames(cross))
#' @param cross - a R/qtl cross object
#' @param method - which software to use: concorde OR lkh
#' @param execPath - path to LKH executable
#' @param chr - a vector of chromosome names to process, default is all
#' @return A new cross object with re-ordered markers for each chromosome
#' @import qtl
#' @export
tspOrder2 <- function(cross, method="concorde", execPath, chr=chrnames(cross))
{
  require(qtl)
  cross = est.rf(cross)
  tmpdir = "tsptemp"
  if(!dir.exists(tmpdir)) dir.create(tmpdir)
  neworder = lapply(chr, function(x){
    cat("Processing chromosome", x, "\n")
    filestem = paste0("./", tmpdir, "/temp_", x)
    
    # pull rfmatrix from the cross
    rf0 = pull.rf(cross, "rf", x)
    mn = rownames(rf0) # marker names
    nmarker = nrow(rf0)
    
    # Create the TSP matrix file.
    cat("\nCreating tsp problem file\n")
    subtspmat = createTSPFile(rf0 , filestem)
    
    # Call concorde
    if (method == "concorde"){
      cat("\nRun concorde\n")
      callConcorde(execPath, paste0(filestem, ".tsp"))
      cat("\nProcess concorde output\n")
      good.order = processConcordeOutput(paste0(filestem, ".sol"))
    }
    
    # Call LKH
    if (method == "lkh"){
      callLKH(execPath, paste0(filestem, ".tsp"))
      cat("\nProcess LKH output\n")
      good.order = processLKHOutput(paste0(filestem, ".LKH"), nmarker)
    }
    
    
    # return ordered marker names
    return(mn[good.order])
  })
  
  names(neworder) = chr
  newcross = reorder.marker(cross, neworder)
  # Remove the temp files.
  unlink(tmpdir, recursive = T)
  return(newcross)
}


# function to reorder a map manually
reorder.marker = function(cross, new.marker.order) {# new.marker.order: a list of chromosomes with marker numbers or marker names
  # require(ASMap)
  cc = names(new.marker.order)
  for (chr in cc){
    genodata = cross$geno[[chr]]$data
    mapdata = cross$geno[[chr]]$map
    neworder = new.marker.order[[chr]]
    new.geno = genodata[,neworder]
    new.map = mapdata[neworder]
    cross$geno[[chr]]$data = new.geno
    cross$geno[[chr]]$map = new.map
  }
  cross2 = quickEst(cross)
  # print(pull.map(cross2, chr, as.table = T))
  return(cross2)
}


##### functions
# check orientation compared to a good order
# return a new cross with good order
matchOrientation = function(cross, cross0){ # cross0 is the one with good order
  chroms = chrnames(cross)
  toflip = c()
  for (i in chroms){
    oldmap = pull.map(cross0, chr = i, as.table = T)
    newmap = pull.map(cross, chr = i, as.table = T)
    mm = match(rownames(oldmap),rownames(newmap))
    cc = cor(1:nrow(oldmap), mm)
    if (cc < 0) toflip = c(toflip, i)
  }
  if(length(toflip)==0) return(cross)
  else {
    newcross = flip.order(cross,toflip)
    plotMap(newcross, cross0)
    return(newcross)
  }
}

# add marker physical chromosomes from markers like "IWB39261_2A_635.6" 
addmarkerchr = function(maps){
  aa = strsplit(rownames(maps), "_")
  bb = sapply(aa, function(x) as.numeric(x[3])) # Mb
  cc = sapply(aa, function(x) x[2]) # chrom names
  maps$phy.chr = cc
  maps$phy.pos = bb
  return(maps)
}

# map region compare, return a boolean whether need2switch
mapcmp = function(map1, map2, start=1, end=NULL){
  if(is.null(end)) end = nrow(map1)
  map.part1 = map1[start:end,]
  map.part2 = map2[start:end,]
  need2switch = F
  if(sum(map.part1$phy.chr==map.part1$chr)>1){
    pos1 = map.part1$phy.pos[map.part1$phy.chr==map.part1$chr]
    pos2 = map.part2$phy.pos[map.part2$phy.chr==map.part2$chr]
    if (all(pos1==pos2)) return(F)
    ss = sort(pos1)
    cat("cor1 is",cor(ss, pos1), "cor2 is",cor(ss, pos2),"\n")
    if (cor(ss, pos1)<cor(ss,pos2)) need2switch = T
  }
  return(need2switch)
}

# compare two orders and get a mixed order that is more consistent to the physical map
better.order = function(map1, map2){
  nmarker = nrow(map1)
  map1 = addmarkerchr(map1)
  map2 = addmarkerchr(map2)
  newmap = map1
  dif = F
  start = 1
  end = 1
  nm1 = rownames(map1)
  nm2 = rownames(map2)
  if (!all(sort(nm1) == sort(nm2))) stop("map1 and map2 do not have the same markers")
  for (i in 1:(nrow(map1)+1)){
    chrom = ifelse(i>nrow(map1), "dum", as.character(map1$chr[i]))
    lastchrom = ifelse(i==1, as.character(map1$chr[1]), as.character(map1$chr[i-1]))
    cat("i is",i,"start is",start, "end is",end, lastchrom, chrom, "\n")
    if (chrom == lastchrom & nm1[i] != nm2[i]){
      # if (dif) end = i
      # else {
      #   dif = T
      #   start = i
      #   end = i
      # }
      dif = T
      end = i
      
    } else {# marker name nm1 == nm2 or enter a new chromosome
      if (end > start){# there are difference before this
        part1 = nm1[start:end]
        part2 = nm2[start:end]
        if (all(sort(part1) == sort(part2))){# can compare phy order
          if (mapcmp(map1,map2, start, end)){
            cat("replace map segment when i is",i,"\n")
            newmap[start:end,] = map2[start:end,]
            row.names(newmap)[start:end] = row.names(map2)[start:end]
          }
          # reset
          dif = F
          start = i
          end = i
        } else end = i
      } else if (dif) end = i
    }
  }
  return(newmap)
}


## modified functions from TSPmap package and Rqtltools
callLKH <- function(execpath, filename)
{
  
  # Get full paths for each file path
  execpath = path.expand(execpath)
  filename = path.expand(filename)

  # Need to create LKH input file (.par), which references the TSP file created by the call to createTSPFile.
  parfilename = sub(".tsp", ".par", filename, fixed = T)
  LKHoutputfilename = sub(".tsp", ".LKH", filename, fixed = T)
  zz <- file(parfilename, "w")
  cat(paste("PROBLEM_FILE = ", filename, sep=""),
      "RUNS = 1",
      paste0("SEED = ",sample(1:1000000,1)),
      paste("OUTPUT_TOUR_FILE = ", LKHoutputfilename, sep = ""),
      "TRACE_LEVEL = 0",
      file = zz, sep = "\n")
  close(zz)

  # Replace all spaces with "\\\\ " before performing the system calls.
  execpath = gsub(" ", "\\\\ ", execpath)
  filename = gsub(" ", "\\\\ ", filename)
  parfilename = gsub(" ", "\\\\ ", parfilename)
  
  system(paste(execpath, file=parfilename), wait=TRUE)
  return(LKHoutputfilename)
}

#' Create TSP problem file.
#'
#'  Create TSP problem file (for use with both Concorde and LKH).
#'  This converts the rf matrix to an integer matrix that the TSP solvers need.
#'  HAMILTONIAN PATH CORRECTION: This also adds a dummy node which has 0 weight to all other nodes so that the TSP solver will find a Hamiltonian path instead of a circuit.
#' @usage createTSPFile(rfmatrix,strainName, dirname)
#' @param rfmatrix recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param strainName name of the strain as a string (no spaces), used to name the file
#' @param dirname directory in which to create the file
#' @return Returns the path and filename of the TSP file (.tsp extension).
#' @export
createTSPFile <- function(rfmatrix, strainName)#, dirname)
{
  # Disable scientific notation so the .tsp files does not have character strings in it.
  options(scipen=999)
  
  # Update numcols with the new value, since we have removed duplicate markers.
  numcols = dim(rfmatrix)[2]
  
  # All values need to be turned into ints so that LKH can read them in, so multiply by 100000 so that we can preserve 5 significant digits (i.e. 0.12345 becomes 12345).  Any factor larger than this seems to cause problems in LKH.
  tspmatrix <- round(rfmatrix*100000)
  
  # Create data file that will be sent to Concorde.
  filename = paste0(strainName, ".tsp")
  zz <- file(filename, "w")
  
  cat("NAME: TSP",
      "COMMENT: Generated by write_TSPLIB (R-package TSP)",
      "TYPE: TSP",
      paste("DIMENSION:", dim(rfmatrix)[1]+1),
      "EDGE_WEIGHT_TYPE: EXPLICIT",
      "EDGE_WEIGHT_FORMAT: UPPER_ROW",
      "EDGE_WEIGHT_SECTION",
      file = zz, sep = "\n")
  
  for(i in 1:(dim(tspmatrix)[1]-1))
  {
    # Start at i+1 to exclude the diagonal elements.
    # Also append a zero at the end to represent the dummy node.
    outputthis = c(tspmatrix[(i+1):dim(tspmatrix)[1],i],0)
    cat(outputthis, "\n", file = zz)
  }
  # Add a final zero for the dummy node.
  # Apply Hamilton path modification.
  # Instead of having the TSP solver find a Hamiltonian circuit, we will add a dummy node with rf = 0 at all elements so that the TSP solver will produce a path that does not have to connect head-to-tail.
  # NOTE: this marker is removed from the TSP solution in the processLKHOutput/processConcordeOutput methods.
  cat(0, "EOF", file = zz, sep = "\n")
  close(zz)

  # Reset scientific notation to the default value.
  options(scipen=0)
  
  return(filename)
}

#' Process LKH solution file.
#'
#'    Process LKH solution file.
#' @usage processLKHOutput(filename, markerdata, rfmat)
#' @param filename - path and name of LKH solution file, which is returned from the function callLKH (this file has a .LKH extension)
#' @return
#'    the best order:

#' @export
processLKHOutput <- function(filename, nmarker)
{
  # The LKH output file has a header that needs to be removed after we read it in.
  # The path terminates with a "-1".
  
  # Skip the first 6 lines since these are the header.
  # Only read in length(markerdata) values since this is the length of the path.  This avoids the problem of reading in the terminal -1 and "EOF".
  # Add 1 to the value of n because of the dummy marker.
  path <- scan(file=filename, skip=6, n = nmarker+1)
  
  # We want to remove the dummy marker from path.  This marker always has the largest index value, which is 1 greater than the length of markerdata.
  maxpos = which(path==max(path))
  if (maxpos == 1 | maxpos == length(path)) {# at either end of the vector
    newpath = path[-maxpos]
  } else {
    newpath = path[c((maxpos+1):length(path), 1:(maxpos-1))]
  }
  return(newpath)
}



#' @title Method to improve a genetic map.
#'
#' @description
#' \code{dropSimilarMarkers} finds markers that have a small recombination fraction and
#' drops the one with combined greater segregation distortion and/or missing data. ***Note:
#' if the cross object has many markers (>1000), avoid running this function on more than one
#' chromosome at a time. Also make sure to run est.rf first and use re.est.map = FALSE***
#'
#' @param cross The qtl cross object to search
#' @param chr numeric, Should the analysis be restricted to a set of chromosomes.
#' Defaults to all chromosomes in the cross object.
#' @param rf.threshold The recombination fraction threshold to drop a marker. If est.rf has
#' not been run on cross, it will be done so automatically. See qtl::est.rf for details
#' @param sd.weight The weighting of segregation distortion rank in dropping a marker.
#' Higher values relative to na.weight increase the weight of the sd rank. Setting a value
#' of 0 removes sd as a factor in choosing the best marker.
#' @param na.weight Same as sd.weight, but for the number of NAs.
#' @param keepEnds Logical, should markers on the ends of the chromosomes always be retained?
#' @param doNotDrop Character vector of markers to retain no matter their rfs.
#' @param verbose Logical, should updates be printed?
#' @param blockSize If not NULL, do an initial culling by splitting markers into
#' blocks of this size. Smaller blocks run more quickly than large blocks, but when the total
#' number of blocks excedes ~ 2000, it can take a very long time to parse the cross object
#' into blocks.
#' @param byChr Should the procedure be run chromosome-by-chromosome. If blocksize != NULL,
#' this procedure is run following block-wise culling. If there are many thousands of markers
#' it is recommended to run multiple block-wise calls prior to whole-chromosome procedures. In
#' general, chromosomes with > 1k markers should first be culled using blockSize != NULL.
#' @param runFullMatrix should the full matrix ever be assessed?
#'
#' @param ... if recombination fractions are not included in the cross object,
#' pass on additional arguments to est.rf.
#'
#' @return A new cross object with the similar markers dropped.
#'
#' @examples
#' set.seed(42)
#' map<-sim.map(len = c(50,20), n.mar = c(20,30), include.x=FALSE)
#' cross0<-sim.cross(map, n.ind=50, type="f2", map.function="kosambi",
#'    error.prob=.01, missing.prob = .05)
#' cross0<-est.rf(cross0)
#' cross1<-dropSimilarMarkers(cross0)
#' cross2<-dropSimilarMarkers(cross0, keepEnds=TRUE)
#' par(mfrow=c(2,1))
#' plot.map(cross0, cross1, main = "comparison of full and culled maps")
#' plot.map(cross0, cross2, main = "comparison of full and culled maps")
#'
#' @import qtl
#' @export

dropSimilarMarkers<-function(cross,
                             chr = NULL,
                             rf.threshold=0.01,
                             sd.weight=1,
                             na.weight=1,
                             keepEnds = FALSE,
                             doNotDrop = NULL,
                             verbose=TRUE,
                             blockSize = 100,
                             byChr = TRUE,
                             runFullMatrix = FALSE,
                             ...){
  loadNamespace("qtl")
  dsm<-function(cross,
                rf.threshold=0.01,
                sd.weight=1,
                na.weight=1,
                keepEnds = FALSE,
                doNotDrop = NULL,
                verbose=TRUE){
    # 1. Get the rfs, geno table and chromosomes in order
    gt<-geno.table(cross)
    if(!"rf" %in% names(cross)){
      cross<-est.rf(cross)
    }
    rf<-pull.rf(cross, what = "rf")
    rf[!upper.tri(rf)]<-1
    diag(rf)<-1
    
    # 1.1 Fancy calculation of sd.weight
    if(class(cross)[1]=="4way"){
      gt.names<-colnames(gt)[3:6]
      gt$P.value<-apply(gt[,gt.names],1,function(x)
        min(x, na.rm = T)/nind(cross))
    }
    gt$rank.p<-with(gt, rank(rank(-P.value, ties.method = "min")*sd.weight))
    gt$rank.sd<-with(gt, rank(rank(missing, ties.method = "min")*na.weight))
    gt$rank<-with(gt, rank(rank.p+rank.sd))
    
    # 2. drop the markers to retain from the matrix
    if(!is.null(doNotDrop)){
      dnd.index<-which(colnames(rf) %in% doNotDrop)
      rf<-rf[-dnd.index, -dnd.index]
    }
    if(keepEnds){
      tokeep<-as.character(
        unlist(
          lapply(pull.map(cross),function(x)
            c(names(x)[1], names(x)[length(x)]))))
      ends.index<-which(colnames(rf) %in% tokeep)
      rf<-rf[-ends.index, -ends.index]
    }
    
    # 3. Loop through the rf matrix, dropping one of the two markers with the lowest
    # recombination fraction.
    nmarstart<-sum(nmar(cross))
    while(min(rf)<rf.threshold){
      worst<-colnames(rf)[which(rf == min(rf, na.rm=TRUE), arr.ind=T)[1,]]
      gtm<-gt[worst,]
      badmars<-rownames(gtm)[which.max(gtm$rank)[1]]
      
      which.todrop<-which(colnames(rf) == badmars)
      rf<-rf[-which.todrop,-which.todrop]
      
      cross<-drop.markers(cross, markers = badmars)
    }
    nmarend<-sum(nmar(cross))
    return(cross)
  }
  
  if(!is.null(chr)) {
    cross = subset(cross, chr = chr)
  }
  if(!is.null(blockSize)){
    if(verbose) cat("initial n markers =", totmar(cross),"\n")
    spl<-split(markernames(cross), ceiling(seq_along(markernames(cross))/blockSize))
    if(length(spl)>2000) warning("breaking cross into > 2k blocks can be very slow\n")
    if(verbose) cat("parsing cross object into blocks\n")
    # temp.cross<-newLG(cross, markerList=spl)
    temp.cross<-reorder.marker(cross, spl)
    if(verbose) cat("running on", length(spl), "blocks of", blockSize, "markers ... \nblock: ")
    goodMars<-lapply(chrnames(temp.cross), function(x){
      ctp<-ifelse(length(spl)>1000, 100, ifelse(length(spl)>100,10, ifelse(length(spl)>50,5,1)))
      if(which(chrnames(temp.cross) == x) %% ctp == 0) cat(x,"")
      cr<-subset(temp.cross, chr = x)
      if(totmar(cr)>1){
        cr<-dsm(cr, rf.threshold = rf.threshold,
                sd.weight = sd.weight,verbose = FALSE,
                keepEnds = keepEnds,
                doNotDrop = doNotDrop)
      }
      return(markernames(cr))
    })
    if(verbose) cat("\n")
    toKeep<-unlist(goodMars)
    toDrop<-markernames(cross)[!markernames(cross) %in% toKeep]
    cross<-drop.markers(cross, markers = toDrop)
    if(verbose) cat("n markers after block-wise culling:", totmar(cross),"\n")
  }
  
  if(byChr){
    if(verbose) cat("running on each chromosome\ninitial n markers:", nmar(cross),"\n")
    if(verbose) cat("Chromosome:")
    goodMars<-lapply(chrnames(cross), function(x){
      if(verbose) cat(x,"")
      cr<-subset(cross, chr = x)
      if(totmar(cr)>1){
        cr<-dsm(cr, rf.threshold = rf.threshold,
                sd.weight = sd.weight,verbose = FALSE,
                keepEnds = keepEnds,
                doNotDrop = doNotDrop)
      }
      return(markernames(cr))
    })
    toKeep<-unlist(goodMars)
    toDrop<-markernames(cross)[!markernames(cross) %in% toKeep]
    cross<-drop.markers(cross, markers = toDrop)
    if(verbose) cat("\nn markers after chromosome-wise culling:", nmar(cross),"\n")
  }
  if(runFullMatrix){
    if(verbose) cat("running for the whole matrix of",totmar(cross),"markers\n")
    cross<-dsm(cross, rf.threshold = rf.threshold,
               sd.weight = sd.weight,verbose = FALSE,
               keepEnds = keepEnds,
               doNotDrop = doNotDrop)
    if(verbose) cat("n markers after whole-map culling:",totmar(cross),"\n")
  }
  return(cross)
}



#' Run Concorde TSP solver.
#'
#'    Run Concorde TSP solver.
#' @usage callConcorde(execpath,filename, outputdir)
#' @param execpath - path and name of Concorde executable.
#' @param filename path and input filename as a string, a .tsp file
#' @param outputdir - directory in which Concorde will save the solution file
#' @return Returns the path and filename of the Concorde solution file (.sol extension).
#' @export
callConcorde <- function(execpath, filename)
{
  execpath = path.expand(execpath)
  filename = path.expand(filename)
  # Replace all spaces with "\\\\ "
  execpath = gsub(" ", "\\\\ ", execpath)
  filename = gsub(" ", "\\\\ ", filename)
  
  # Create the output file name (.sol file extension).
  outfilename = sub(".tsp", ".sol", filename, fixed = T)
  
  # Need to add these parameters to the concorde call:
  #   -x    delete files on completion (sav pul mas)
  #   -o f  output file name (for optimal tour)
  exec = paste(execpath, '-x -o', outfilename)
  cat("CMD is:", exec, "\n")
  
  system(paste(exec, file=filename), wait=TRUE)
  
  return(outfilename)
}


#'Process Concorde solution file
#'
#'    Process Concorde solution file.
#' @usage processConcordeOutput(filename, markerdata, rfmat)
#' @param filename path and name of Concode solution file, which is returned from the function callConcorde (this file has a .sol extension)
#' @return
#'    new the marker order
#' @export
processConcordeOutput <- function(filename)
{
  # IMPORTANT: Concorde considers the first node to be #0, so we will need to adjust the indices to match those of the input file.
  
  # Read in Concorde's solution file.
  path <- scan(file=filename)
  
  # Remove first element of list, since this is just the number of nodes.
  path = path[-1]
  
  # Increment all values in path by 1 because Concorde starts indexing at 0 but R starts at 1.
  path <- path + 1
  
  # We want to remove the dummy marker from path.  This marker always has the largest index value, which is 1 greater than the length of markerdata.
  maxpos = which(path==max(path))
  if (maxpos == 1 | maxpos == length(path)) {# at either end of the vector
    newpath = path[-maxpos]
  } else {
    newpath = path[c((maxpos+1):length(path), 1:(maxpos-1))]
  }
  return(newpath)
}

#' Fast estimate the genetic map for an order of markers
#' from package ASMap
quickEst <- function(object, chr, map.function = "kosambi", ...){
  if (!any(class(object) == "cross"))
    stop("Input should have class \"cross\".")
  if (missing(chr))
    chr <- names(nmar(object))
  imf <- switch(map.function, kosambi = imf.k, haldane = imf.h,
                morgan = imf.m, cf = imf.cf)
  nm <- nmar(object)
  for(i in chr){
    temp <- subset(object, chr = i)
    if(nmar(temp) != 1){
      est <- est.rf(temp)$rf
      nc <- dim(est)[1]
      er <- est[cbind(2:nc,1:(nc - 1))]
      temp$geno[[i]]$map <- c(0,cumsum(imf(er)))
      names(temp$geno[[i]]$map) <- dimnames(temp$geno[[i]]$data)[[2]]
      tempa <- argmax.geno(temp, step = 0, map.function = map.function, ...)
      tempa$geno[[i]]$data <- tempa$geno[[i]]$argmax
      tempa$geno[[i]] <- tempa$geno[[i]][-3]
      esta <- est.rf(tempa)$rf
      era <- esta[cbind(2:nc,1:(nc - 1))]
      if(class(object)[1] == "riself")
        era <- (era/2)/(1 - era)
      object$geno[[i]]$map <- c(0,cumsum(imf(era)))
      names(object$geno[[i]]$map) <- dimnames(object$geno[[i]]$data)[[2]]
    }
  }
  object
}

