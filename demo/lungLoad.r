
load.gsea <- function(pathwayFile="c2.cp.kegg.v3.0.symbols.gmt",cohort = "Boston", path, npars=2, ncats=3) {

  lines <- readLines(pathwayFile)
  n <- length(lines)
  if(is.character(path)) { 
    for(ll in lines) {
      ls <- strsplit(ll,"\t")[[1]]
      if(ls[1] == path)
        break
    }
    if(length(ls) < 1 || ls[1] != path)
      stop("no valid pathway")
  }
  else if(is.numeric(path)) {
    path <- as.integer(path)
    if(path < 1 || path > n)
      stop("n out of range")
    ls <- strsplit(lines[[path]],"\t")[[1]]
    cat(ls[1],"\n")
  } 
 if(length(ls) < 3)
    stop("no valid pathway")
  pathname <- ls[1]
  path <- ls[3:length(ls)]
    
  if(cohort == "Boston") {
    data <- read.table("Lung_Boston_collapsed_symbols.gct", header=TRUE, sep='\t', skip=2)
    cls <- c(rep(1, 31), rep(2,31))
  }
  else if(cohort == "Michigan") {
    data <- read.table("Lung_Michigan_collapsed_symbols.gct", header=TRUE, sep='\t', skip=2)
    cls <- c(rep(1, 24), rep(2, 62))
  }
  else {
    stop("no valid cohort")
  }
  
  pathid <- NULL
  levs <- levels(data[,1])
  for(i in 1:length(path)) {
    id <- which(data[,1] == path[i])
    if(length(id) != 1)
      next
    str <- levs[data[id[1], 1]]
    if(nchar(str) > 16 || length(which(pathid==id[1]))>0)
      next
    pathid <- c(pathid, id[1])
  }
  path <- as.character(data[pathid,1])

  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("-",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("@",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub(" /// ",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))
  path <- as.character(sapply(path, function(ss) sub("_",".",ss)))

  cdata <- data[pathid, 3:ncol(data)]
  rownames(cdata) <- path
  numnodes <- nrow(cdata)
  numsamples <- ncol(cdata)
  cdata <- as.matrix(cdata, nrow=numnodes)
    
  cat(numnodes, " nodes and ", numsamples, " samples\n")
  
  nodePars <- rep(npars, numnodes)
  nodeOrder <- 1:numnodes
  nodeCats <- lapply(1:nrow(cdata), function(i) return(1:ncats))
  names(nodeCats) <- path
  
  return(list(pathname=pathname, cdata=cdata, cls=cls, pathway=path, nodeOrder=nodeOrder, nodePars=nodePars, nodeCats=nodeCats))
}

