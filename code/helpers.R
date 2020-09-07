

# convert a list of dataframe into a long dataframe
list_to_df <- function(mylist){
  
  # make a vector of row ids that correspond to the list names
  rowid.indx <- lapply(mylist, function(x) dim(x)[1])
  sourceVec.list <- list()
  for(i in 1:length(rowid.indx)){
    sourceName <- names(rowid.indx)[i]
    numRows <- rowid.indx[[i]]
    sourceVec.list[[i]] <- rep(sourceName, numRows)
  }
  rowVec <- unlist(sourceVec.list)
  
  # combine into df
  df <- data.frame(do.call(rbind, mylist), row.names = NULL)
  df$source <- rowVec
  
  return(df)
}

# source a bunch of files from the same folder; code taken from R help -> source()
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# create an informative unique name for each OTU from the RDP taxonomy table
add_niceOTUname <- function(curr.tab.taxa){
  
  tmp <- curr.tab.taxa
  tmp$best.name <- NA
  taxlevels <- c("kingdom","phylum","class","order","family", "genus")
  for(i in 1:dim(tmp)[1]){ #loop through each row in the df
    curr.row <- tmp[i, taxlevels]
    highest.pos <- max(which(curr.row != "unidentified")) # highest classified position
    highest.pos
    best.name <- paste("unclassified", curr.row[,taxlevels[highest.pos]], sep ="_")
    best.name
    tmp[i,"best.name"] <- best.name
  }
  
  # if genus level confidence is acceptable, use the species field as the best name
  tmp %>%
    mutate(best.name = ifelse(genus != "unidentified", 
                              genusSpecies, 
                              best.name)) -> tmp.new
  
  return(tmp.new)
  
}

