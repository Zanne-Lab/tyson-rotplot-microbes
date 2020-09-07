# all microbes

#------------------------------------------------------------------#
# Tidy the raw OTU tables

load_raw_otu <- function(otu.name){
  
  fileName.part1<-"./data/microbes/otuTables/otuTable_"
  fileName.part2<- ".tab"
  fileName<-paste(fileName.part1, otu.name, fileName.part2, sep = "")
  raw.otu<- read.delim(fileName,
                       header=TRUE,                 #make the barcodes the headers
                       stringsAsFactors = FALSE,
                       strip.white = TRUE,          #column names have trailing white space
                       row.names = 1)               #make the OTU ids the row names
  raw.otu <- as.matrix(raw.otu)

  return(raw.otu)
}

remove_zeroOTUcols <- function(curr.otuTab){
  
  #check if there are entirely empty OTU columns and prune them out
  num0colSums <- sum(colSums(curr.otuTab) == 0)
  if(num0colSums > 0){
    print(paste("needed to prune ", num0colSums, " columns", sep=""))
    curr.otuTab <- curr.otuTab[, colSums(curr.otuTab) != 0]
  }
  return(curr.otuTab)
}

remove_zeroSamprows <- function(curr.otuTab){
  
  #check if there are entirely empty Sample rows and prune them out
  num0rowSums <- sum(rowSums(curr.otuTab) == 0)
  if(num0rowSums > 0){
    print(paste("needed to prune ", num0rowSums, " rows", sep=""))
    curr.otuTab <- curr.otuTab[rowSums(curr.otuTab) != 0, ]
  }
  return(curr.otuTab)
}


#------------------------------------------------------------------#
# Subset OTU table information

SelectSamps <- function(select.samps, otuTab.list, taxon.list){
  
  #identify the selected OTU table
  tableName <- unique(select.samps$table)
  curr.otuTab <- otuTab.list[[as.character(tableName)]]
  
  #identify selected samples within the OTU table
  select.samps$tidySeq_sampleName
  row.names(curr.otuTab)
  x <- rownames(curr.otuTab) %in% select.samps$tidySeq_sampleName
  curr.otuTab <- curr.otuTab[x,]

  length(dim(curr.otuTab)) == 2
  if(length(dim(curr.otuTab)) == 2){
    #remove OTUs that do not appear at all in the subset
    curr.otuTab <- remove_zeroOTUcols(curr.otuTab = curr.otuTab)
    #transform OTU table, make long
    curr.otuTab.t <- data.frame(t(curr.otuTab))
    curr.otuTab.t$OTUid <- row.names(curr.otuTab.t)
    
  }else{
    #remove OTUs that do not appear at all in the subset
    curr.otuTab <- curr.otuTab[curr.otuTab != 0]
    # make OTU table into a long df
    curr.otuTab.t <- data.frame(OTUid = names(curr.otuTab), curr.otuTab)
    colnames(curr.otuTab.t)[2] <- select.samps$tidySeq_sampleName
  }
  
  #add taxon info
  curr.taxonTab <- taxon.list[[as.character(tableName)]]
  curr.otuTab.t %>%
    gather(key = "tidySeq_sampleName", value = "numReads", -OTUid) %>%
    filter(numReads != 0) %>%
    left_join(curr.taxonTab) -> longtab
  
  #remove extra columns
  #hier<-c("kingdom","phylum","class","order","family","genus","species")
  gene <- unique(select.samps$gene)
  if(gene == "16S" | gene == "16SA" | gene == "16SB"){
    keepcols<-c("tidySeq_sampleName","numReads","OTUid","genusSpecies", "phylum")
  }else{
    keepcols<-c("tidySeq_sampleName","numReads", "OTUid","genusSpecies", "phylum",
                "Trophic.Mode","Trait")
  }
  longtab %>%
    select(keepcols) %>%
    arrange(tidySeq_sampleName, OTUid, desc(numReads)) -> longtab
  
  return(longtab)
}

SelectSamps.OTUtab <- function(select.samps, otuTab.list){
  
  #identify the OTU tables in select.samps
  TABLES <- names(otuTab.list)[names(otuTab.list) %in% unique(select.samps$table)]
  new.otuTab.list <- list()
  for(i in 1:length(TABLES)){ #loop through each OTU table
    
    #identify the current OTU table and samples
    select.samps %>%
      filter(table == TABLES[i]) -> curr.select.samps
    curr.otuTab <- otuTab.list[[TABLES[i]]]
    
    #subset the OTU table by the select samples
    x <- row.names(curr.otuTab) %in% curr.select.samps$tidySeq_sampleName
    new.otuTab <- curr.otuTab[x,]
    
    #check if there are entirely empty OTU columns and prune them out
    new.otuTab <- remove_zeroOTUcols(curr.otuTab = new.otuTab)
    
    #save the new OTU table
    new.otuTab.list[[i]] <- new.otuTab
  }
  names(new.otuTab.list) <- TABLES
  
  return(new.otuTab.list)
}


#------------------------------------------------------------------#
# Subtract negative controls from OTU tables (and examine select samples, e.g. positive controls)

RemoveNegatives <- function(negTab, otuTab){
  
  for (i in 1:dim(negTab)[1]){ #loop through each row in negTab
    
    #identify the otu and number of reads to subtract
    curr.OTUid<-negTab[i,"OTUid"]
    curr.subt<-negTab[i,"numReads"]
    
    #do the subtraction and update the column values
    col.tmp <- otuTab[, colnames(otuTab) == curr.OTUid] - curr.subt
    col.tmp[col.tmp < 0] <- 0 # change negative read counts to 0
    otuTab[, colnames(otuTab) == curr.OTUid] <- col.tmp
    
  }
  
  #check if there are entirely empty OTU columns and prune them out
  otuTab <- remove_zeroOTUcols(curr.otuTab = otuTab)
  
  return(otuTab)
}


#------------------------------------------------------------------#
# Trim failed samples from OTU tables

HistoSeqsPerSamp <- function(curr.otuTab.list, maxXval){
  
  hist.list <- list()
  for(i in 1:length(curr.otuTab.list)){
    otuTab <- curr.otuTab.list[[i]]
    sampReads <- data.frame(OTUid = row.names(otuTab), numReads = rowSums(otuTab))
    p <- ggplot(sampReads, aes(x = numReads)) + 
      geom_freqpoly() + 
      scale_x_continuous(limits = c(0, maxXval)) + 
      geom_vline(xintercept=20, linetype=2) + 
      geom_vline(xintercept=100, linetype=2, color=3) +
      geom_vline(xintercept=500, linetype=2, color=4) +
      geom_vline(xintercept=1000, linetype=2, color=5) +
      ggtitle(names(otuTab.list)[i])
    hist.list[[i]] <- p
  }
  names(hist.list) <- names(curr.otuTab.list)
  
  
  return(hist.list)
}

RemoveSamps.OTUtab <- function(minReads, curr.otuTab.list){
  
  #loop through each table in the curr.otuTab.list
  new.otuTab.list<-list()
  for(i in 1:length(curr.otuTab.list)){
    
    #identify the OTU table
    curr.otuTab<-curr.otuTab.list[[i]]
    
    #remove rows with rowsSums less than minReads
    new.otuTab <- curr.otuTab[rowSums(curr.otuTab) > minReads,]
    
    #check if there are entirely empty OTU columns and prune them out
    new.otuTab.list[[i]] <- remove_zeroOTUcols(new.otuTab)
  }
  names(new.otuTab.list) <- names(curr.otuTab.list)
  
  return(new.otuTab.list)
}

#------------------------------------------------------------------#
# Rarify, normalize

MakeRare.OTUtab <- function(curr.otuTab.list, sampleSize){
  
  require("phyloseq")
  
  curr.otuTab.list.rare <- list()
  for(i in 1:length(curr.otuTab.list)){
    
    #identify the OTU table
    curr.otuTab <- curr.otuTab.list[[i]]
    tableName <- names(curr.otuTab.list)[i]
    
    #put the OTU table into phyloseq format and rarify
    otumat <- t(curr.otuTab)
    OTU <- otu_table(otumat, taxa_are_rows = TRUE)
    OTU.rare <- rarefy_even_depth(OTU, trim=TRUE, rngseed=711, sample.size = sampleSize)
    new.otuTab <- t(OTU.rare@.Data)
    
    #remove empty otu cols and save
    curr.otuTab.list.rare[[i]] <- remove_zeroOTUcols(new.otuTab)
    
  }
  names(curr.otuTab.list.rare)<-names(curr.otuTab.list)
  
  return(curr.otuTab.list.rare)
}

MakeCSSNorm.OTUtab <- function(curr.otuTab.list){
  
  require("metagenomeSeq")
  
  curr.otuTab.list.norm <- list()
  for(i in 1:length(curr.otuTab.list)){
    
    curr.otuTab <- curr.otuTab.list[[i]]
    
    #turn OTU table into a MRexperiment object
    counts <- t(curr.otuTab) # columns need to be samples
    mrobj <- newMRexperiment(counts = counts)
    
    #do the abund normalization
    #p <- cumNormStatFast(mrobj, pFlag = TRUE, rel = 0.1) #defaultsp
    p <- cumNormStat(mrobj, pFlag = TRUE, rel = 0.1)
    
    new.mrobj <- cumNorm(mrobj, p = p)
    otuTab.norm <- MRcounts(new.mrobj, norm = TRUE, log = TRUE)
    curr.otuTab.list.norm[[i]] <- t(otuTab.norm)
  }
  names(curr.otuTab.list.norm) <- names(curr.otuTab.list)
  
  return(curr.otuTab.list.norm)
  
}

#remove OTUs that are present in less than 20% of samples
trimLowSampleOTUs <- function(otuTab){
  x <- 1 * (otuTab > 0) # transform into pres/abs
  minSamps <- floor(dim(x)[1] * 0.2)
  x.trim <- otuTab[,colSums(x) >= minSamps]
  
  #remove samples with zero OTU observations
  x.trim <- x.trim[rowSums(x.trim) > 0,]
  
  return(x.trim)
}



#------------------------------------------------------------------#
# Prep OTU tables
# ... see prep_otutabs_rotplot.R and prep_otutabs_cwd.R 

#------------------------------------------------------------------#
# Read and write tidy otuTabs

write_tidy_otuTabs <- function(folderName, prepName, otuTab.list){
  
  path <- "data_Rsynth/otuTables_tidy"
  path <- paste0(path, "/", folderName)
  
  for(i in 1:length(names(otuTab.list))){
    fileName <- paste0(path,"/", prepName, "_", names(otuTab.list)[i], ".csv")
    write.csv(otuTab.list[[i]], file = fileName)
  }
}

load_tidy_otuTabs <- function(folderName, prepName){
  
  path <- "data_Rsynth/otuTables_tidy"
  path <- paste0(path, "/", folderName)
  
  files <- data.frame(fileName=list.files(path))
  
  if(folderName == "cwd"){
    files %>%
      separate(fileName, into = c("prepName_col","tableName","drop","drop2"), remove = F) %>%
      transform(tableName = paste0(tableName, "_", drop)) %>%
      select(-c(drop, drop2)) %>%
      filter(prepName_col == prepName) -> file.indx
  }else{
    files %>%
      separate(fileName, into = c("prepName_col","tableName","drop"), remove = F) %>%
      select(-drop) %>%
      filter(prepName_col == prepName) -> file.indx
  }
  
  otuTab.list <- list()
  for(i in 1:dim(file.indx)[1]){
    fileName <- paste0(path, "/", file.indx[i, "fileName"])
    otuTab.list[[i]] <- read.csv(file=fileName, row.names = 1)
  }
  names(otuTab.list) <- file.indx$tableName
  
  return(otuTab.list)
  
}

