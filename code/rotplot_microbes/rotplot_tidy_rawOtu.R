# rotplot-specific: requires prep_otutabs_all.R, load_otuTaxa.R, load_rotplotMetaData.R

#------------------------------------------------------------------#
# Tidy the raw OTU tables

MakeSampleLookup_12 <- function(){
  
  #load OTU tables to pull out the sample names
  table.names<-c("2012ITS","2012LSU","201216S")
  otuTab.list <- lapply(table.names, load_raw_otu)
  names(otuTab.list) <- table.names
  barcode.list <- lapply(otuTab.list, colnames)
  
  # rotplot study design info
  logmicrobs <- load_microbes() # load_rotplotMetaData.R
  logmicrobs %>%
    filter(yrharv == 2012) -> logmicrobs.12
  
  #MiSeq 2012 barcode to sample lookup table
  barcodes<-read.delim("./data/studyMetadata/miSeqBarcodes/Identified_Barcodes_samples_Corrected.txt", stringsAsFactors = FALSE)        
  barcodes %>%
    select(Barcode, Sample) %>%
    separate(col=Sample, into=c("sampleName","sampleType"), sep="z", fill = "right", remove=FALSE) %>%
    mutate(sampleType = ifelse(sampleType == "" & !is.na(sampleType), "rotplot", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("Pos_", sampleName), "positiveControl", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("Neg_", sampleName), "negativeControl", sampleType)) %>%
    mutate(sampleType = ifelse(is.na(sampleType), "WMstudy", sampleType)) %>%
    left_join(logmicrobs.12) %>%
    rename("seq_sampleName"="Barcode") %>%
    mutate(yrseq = 2012) %>%
    select(seq_sampleName, sampleName, sampleType, logtrt, plot, logLoc, bag, yrseq, yrharv) -> samplookup.2012
  
  # elongate the table with full list of sequenced samples (e.g. by gene)
  samplookup.2012.ann <- list()
  for(i in 1:length(barcode.list)){
    df <- data.frame(seq_sampleName = barcode.list[[i]], table = names(barcode.list)[i])
    samplookup.2012.ann[[i]] <- left_join(df, samplookup.2012)
  }
  names(samplookup.2012.ann)<-names(barcode.list)
  ann.samplookup.2012 <- list_to_df(samplookup.2012.ann)
  
  return(ann.samplookup.2012)
  
}
MakeSampleLookup_14 <-function(){
  
  #load OTU tables to pull out the sample names
  table.names<-c("2014ITS","2014LSU","201416S")
  otuTab.list <- lapply(table.names, load_raw_otu)
  names(otuTab.list) <- table.names
  samp.list <- lapply(otuTab.list, colnames)
  
  # rotplot study design info
  logmicrobs <- load_microbes() 
  logmicrobs %>%
    filter(yrharv == 2014) -> logmicrobs.14
  
  # make the format of "sampleName" match the sampleNames used in the OTU table
  # for bag samples...
  logmicrobs.14 %>%
    filter(bag == "b") %>%
    separate(col="sampleName", into=c("mainSampleName","bag1"), sep="-", remove = FALSE) %>%
    mutate(sampleName1 = paste(mainSampleName, bag1, sep = "_")) %>%
    select(logtrt, plot, logLoc, bag, rep_code, yrharv, sampleName, sampleName1) -> bagsamps
  # add non-bag samples
  logmicrobs.14 %>%
    filter(is.na(bag)) %>%
    mutate(sampleName1 = sampleName) -> non.bagsamps
  logmicrobs.14.ann <- rbind(bagsamps, non.bagsamps)
  
  #re-format
  logmicrobs.14.ann %>%
    select(-sampleName) %>%
    rename("sampleName"="sampleName1") %>%
    mutate(sampleType = "rotplot") %>%
    mutate(yrseq = 2014) %>%
    select(sampleName, sampleType, logtrt, plot, logLoc, bag, yrseq, yrharv) -> samplookup.2014
  
  # complete the sample lookup table with full list of sequenced samples
  samplookup.2014.ann <- list()
  for(i in 1:length(samp.list)){
    df <- data.frame(sampcode = samp.list[[i]], table = names(samp.list)[i])
    df %>%
      separate(sampcode, into = c("R","gene","sampid"), remove = FALSE, extra = "merge") %>%
      rename("seq_sampleName" = "sampcode",
             "sampleName" = "sampid") %>%
      select(-c(R, gene)) -> df
    samplookup.2014.ann[[i]] <- left_join(df, samplookup.2014)
  }
  names(samplookup.2014.ann)<-names(samp.list)
  ann.samplookup.2014 <- list_to_df(samplookup.2014.ann)
  
  #----------------------------#
  # annotate the metadata for extra samples
  
  # identify sampleType of control samples
  ann.samplookup.2014 %>%
    mutate(sampleType = ifelse(grepl("soil", seq_sampleName), "rotplot_soil", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("Pth", seq_sampleName), "pathogens", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("NTC", seq_sampleName), "negativeControl", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("GTaq", seq_sampleName), "2012Taq_on_2014samples", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("Q5", seq_sampleName), "2014Taq_on_2012samples", sampleType)) %>%
    mutate(sampleType = ifelse(grepl("mockcomm", seq_sampleName), "mockcomm", sampleType)) %>%
    mutate(yrseq = ifelse(!is.na(sampleType), 2014, yrseq)) -> pause
  
  # fill-in rotplot_soil missing metadata
  plotLoc <- load_plotLoc()
  pause %>%
    filter(sampleType == "rotplot_soil", !is.na(sampleType)) %>%
    select(seq_sampleName) %>%
    separate(seq_sampleName, into=c("thing1","thing2"), sep="soil", remove = F) %>%
    separate(thing2, into=c("wshed","topo"), sep = 1, remove = F) %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    select(seq_sampleName, plot) -> soilindx
  soilindx <- unique(soilindx)
  for(i in 1:dim(soilindx)[1]){
    pause[pause$seq_sampleName == soilindx[i,1], "plot"] <- soilindx[i,2]
  }
  pause %>%
    mutate(yrharv =  ifelse(sampleType == "rotplot_soil", 2014, yrharv)) -> pause
  
  # fill-in 2012Taq_on_2014samples missing metadata
  logmicrobs %>%
    filter(bag == "b") %>%
    select(sampleName, logtrt, plot, logLoc) %>%
    mutate(sampleName = gsub("-","_", sampleName)) -> logmicrobs.indx
  pause %>%
    filter(sampleType == "2012Taq_on_2014samples") %>%
    select(seq_sampleName) %>%
    separate(seq_sampleName, into = c("thing1","thing2"), sep = "GTaq", remove = FALSE) %>%
    separate(thing1, into=c("sampleName","thing4"), sep = "bag") %>% #all are bag samples
    mutate(bag = "b") %>%
    mutate(sampleName = paste(sampleName, "bag", sep="_")) %>%
    select(seq_sampleName, sampleName, bag) %>%
    left_join(logmicrobs.indx) -> taq12.indx
  taq12.indx<-unique(taq12.indx)
  for(i in 1:dim(taq12.indx)[1]){
    pause[pause$seq_sampleName == taq12.indx[i,1], c("bag","logtrt","plot","logLoc")] <- 
      taq12.indx[i, c("bag","logtrt","plot","logLoc")]
  }
  pause %>%
    mutate(yrharv =  ifelse(sampleType == "2012Taq_on_2014samples", 2014, yrharv)) -> pause
  
  # fill-in 2014Taq_on_2012samples missing metadata
  pause %>%
    filter(sampleType == "2014Taq_on_2012samples") %>%
    filter(!sampleName %in% c("Q5b","Q5t")) %>%
    separate(seq_sampleName, into=c("thing1","thing2"), sep = "12", remove = FALSE) %>%
    separate(thing2, into=c("thing3","thing4"), sep = "Q5") %>%
    select(-c(thing1, thing4, logtrt, plot, logLoc)) %>%
    separate(thing3, into=c("logtrt","plot","logLoc","rep_code"), sep = c(1,2,3)) %>%
    mutate(rep_code = tolower(rep_code)) %>%
    select(seq_sampleName, logtrt, plot, logLoc, sampleType, yrharv) %>%
    mutate(sampleName = paste0(logtrt, plot, logLoc)) -> tmpindx
  tmpindx <- unique(tmpindx)
  for(i in 1:dim(tmpindx)[1]){
    pause[pause$seq_sampleName == tmpindx[i,1], c("logtrt","plot","logLoc")] <-
      tmpindx[i, c("logtrt","plot","logLoc")]
  }
  pause %>%
    mutate(yrharv =  ifelse(sampleType == "2014Taq_on_2012samples", 2012, yrharv)) -> ann.samplookup.2014
  
  return(ann.samplookup.2014)
  
}
make_miSeq_sampleList <- function(){
  
  ann.samplookup.2012 <- MakeSampleLookup_12()
  ann.samplookup.2014 <- MakeSampleLookup_14()
  
  # combine
  samplookup <- rbind(ann.samplookup.2012, ann.samplookup.2014)
  samplookup %>%
    mutate(tidySeq_sampleName = paste(sampleName, yrseq, sep = "__")) %>%
    separate(table, into = c("drop","gene"), sep = 4, remove = FALSE) %>%
    select(table, gene, tidySeq_sampleName, seq_sampleName, sampleName, sampleType, logtrt, plot, logLoc, bag, yrseq, yrharv) %>%
    arrange(table, sampleType, sampleName) -> samps
  
  # check that all the tidySeq_sampleNames are unique
  # samps %>%
  #   group_by(table) %>%
  #   summarize(total = length(tidySeq_sampleName),
  #             uniq = length(unique(tidySeq_sampleName)))
  
  return(samps)
  
}
make_taqRedo_sampleList <- function(){
  
  # make a sample list and standardize sample names in otu table
  path <- "./data/studyMetaData/miSeqBarcodes/ZanneAug2017_SubmissionSheet_Amplicon_v1.94.csv"
  barcodes17 <- read.csv(file = path)
  barcodes17 %>%
    filter(Project_ID == "TysonReDo") %>%
    select(Sample_ID, Target, Barcode_Name) %>%
    rename("seq_sampleName"="Barcode_Name",
           "gene"="Target",
           "sampleName"="Sample_ID") %>%
    mutate(sampleType = "redos") -> samps.tmp
  
  unique(samps.tmp$gene) #only for ITS
  
  samps.tmp %>%
    mutate(table = "2017ITS") %>%
    mutate(yrseq = 2017) %>%
    separate(sampleName, into = c("drop1","taqType","sampName", "yrharv"), remove = F) %>%
    select(-drop1) %>%
    mutate(yrharv = ifelse(yrharv == 12, 2012, yrharv)) %>%
    mutate(yrharv = ifelse(yrharv == 14, 2014, yrharv)) %>%
    separate(sampName, into = c("logtrt","plot","logLoc"), sep=c(1,2)) %>%
    mutate(bag = NA) %>%
    mutate(tidySeq_sampleName = paste0(taqType, "-", logtrt, plot, logLoc, "-", yrharv, "__", yrseq)) %>%
    select(table, gene, tidySeq_sampleName, seq_sampleName, sampleName, sampleType, logtrt, plot, logLoc, bag, yrseq, yrharv, taqType) -> samps
  
  return(samps)
  
}

standardize_otu_samples <- function(otuTab, tableName, samps){
  
  # identify the samples to keep in the OTU table
  samps %>%
    filter(table == tableName) %>%
    filter(sampleType != "WMstudy") -> samps.select
  
  # subset the OTU table
  x <- colnames(otuTab) %in% samps.select$seq_sampleName
  otuTab.trim <- otuTab[,x]
  
  # replace OTU sample names with tidy ones
  x <- match(colnames(otuTab.trim), samps.select$seq_sampleName)
  o.samps.select<-samps.select[x,] #reorganize the lookup table
  colnames(otuTab.trim) <- o.samps.select$tidySeq_sampleName
  
  # transform table
  otuTab.trim.t<-t(otuTab.trim)
  
  # remove OTUs that no longer show up in the dataset
  otuTab.trim.tt <- remove_zeroOTUcols(otuTab.trim.t)
  
  return(otuTab.trim.tt)
}
batch_standarize_otu_samples <- function(otuTab.list, samps){
  
  otuTab.list.standard <- list()
  for(i in 1:length(otuTab.list)){
    otuTab.list.standard[[i]] <- standardize_otu_samples(otuTab = otuTab.list[[i]], tableName = names(otuTab.list)[i], samps = samps)
  }
  names(otuTab.list.standard) <- names(otuTab.list)
  
  return(otuTab.list.standard)
  
}

#------------------------------------------------------------------#
# Subtract negative controls from OTU tables (and examine select samples, e.g. positive controls)

negTab_2012ITS <- function(samps, otuTab.list, taxon.list){
  
  # 2012 ITS
  samps %>%
    filter(table == "2012ITS") %>%
    filter(sampleType == "negativeControl") -> select.samps
  df <- SelectSamps(select.samps = select.samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
  df %>%
    group_by(tidySeq_sampleName) %>%
    summarize(OTUrich = length(genusSpecies))
  #identify the negative control sequences to subtract out
  df %>%
    filter(tidySeq_sampleName != "Neg_1__2012") %>%
    select(tidySeq_sampleName, numReads, OTUid) -> neg2012ITS
  
  return(neg2012ITS)
  
}
negTab_2012LSU <- function(samps, otuTab.list, taxon.list){
  
  # 2012 LSU
  samps %>%
    filter(table == "2012LSU") %>%
    filter(sampleType == "negativeControl") -> select.samps
  df <- SelectSamps(select.samps = select.samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
  df %>%
    group_by(tidySeq_sampleName) %>%
    summarize(OTUrich = length(genusSpecies))
  #identify the negative control sequences to subtract out
  df %>%
    filter(tidySeq_sampleName != "Neg_1__2012") %>%
    select(tidySeq_sampleName, numReads, OTUid) -> neg2012LSU
  
  return(neg2012LSU)
  
}
negTab_2014ITS <- function(samps, otuTab.list, taxon.list){
  
  #2014 ITS
  samps %>%
    filter(table == "2014ITS") %>%
    filter(sampleType == "negativeControl") -> select.samps
  df <- SelectSamps(select.samps = select.samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
  df %>%
    group_by(tidySeq_sampleName) %>%
    summarize(OTUrich = length(genusSpecies))
  #identify the negative control sequences to subtract out
  df %>%
    select(tidySeq_sampleName, numReads, OTUid) -> neg2014ITS
  
  return(neg2014ITS)
}
negTab_201416S <- function(samps, otuTab.list, taxon.list){
  
  #2014 16S
  samps %>%
    filter(table == "201416S") %>%
    filter(sampleType == "negativeControl") -> select.samps
  df <- SelectSamps(select.samps = select.samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
  df %>%
    group_by(tidySeq_sampleName) %>%
    summarize(OTUrich = length(genusSpecies))
  #identify the negative control sequences to subtract out
  df %>%
    select(tidySeq_sampleName, numReads, OTUid) -> neg201416S
  
  return(neg201416S)
  
}
subtract_negativeControls <- function(otuTab.list, taxon.list, samps){
  
  #initialize new OTU list to hold OTU tables that have had negatives removed
  otuTab.list.negRem <- list()
  
  if(sum(grepl('2012', names(otuTab.list))) != 0){
    # 2012 ITS
    tmp <- negTab_2012ITS(samps = samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
    tmp
    otuTab.list.negRem[["2012ITS"]] <- RemoveNegatives(negTab = tmp, otuTab = otuTab.list[["2012ITS"]])
    
    # 2012 16S
    otuTab.list.negRem[["201216S"]] <- otuTab.list[["201216S"]] #no negative controls for 16S in 2012
  }
  
  if(sum(grepl('2014', names(otuTab.list))) != 0){
    #2014 ITS
    tmp <- negTab_2014ITS(samps = samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
    otuTab.list.negRem[["2014ITS"]] <- RemoveNegatives(negTab = tmp, otuTab = otuTab.list[["2014ITS"]])
    
    #2014 16S
    tmp <- negTab_201416S(samps = samps, otuTab.list = otuTab.list, taxon.list = taxon.list)
    
    # check for Bradyrhizo in negative controls
    # tax <- taxon.list$`201416S`
    # tmp %>%
    #   left_join(tax) -> tmp.tax
    # tmp.tax$genus == "Bradyrhizobium"
    
    otuTab.list.negRem[["201416S"]] <- RemoveNegatives(negTab = tmp, otuTab = otuTab.list[["201416S"]])
  } 
  
  return(otuTab.list.negRem)
  
}

#------------------------------------------------------------------#
# Trim failed samples from OTU tables

plot_libSize_distributions <- function(otuTab.list, select.samps, maxXval, ammendfileName){
  
  require(gridExtra)
  
  curr.otuTab.list <- SelectSamps.OTUtab(select.samps = select.samps, otuTab.list = otuTab.list)
  
  #plot distribution of the number of reads per sample
  hist.list <- HistoSeqsPerSamp(curr.otuTab.list = curr.otuTab.list, maxXval=maxXval)
  
  #save plot
  path <- "output/prep_otuTabs/min_readsPerSample_"
  fileName <- paste0(path, ammendfileName, ".pdf")
  pdf(file = fileName)
  grid.arrange(hist.list[["2012ITS"]],
               hist.list[["2012LSU"]],
               hist.list[["201216S"]],
               hist.list[["2014ITS"]],
               hist.list[["2014LSU"]],
               hist.list[["201416S"]],
               nrow=2)
  dev.off()
  
}

TooFewReads <- function(minReads, curr.otuTab.list, select.samps){
  
  # first, subset the OTUs by the select.samps
  curr.otuTab.list <- SelectSamps.OTUtab(select.samps = select.samps, otuTab.list = curr.otuTab.list)
  
  # next, loop through each OTU tab and compile a loss list
  loss.list<-list()
  for(i in 1:length(curr.otuTab.list)){
    
    #select the OTU table
    curr.otuTab <- curr.otuTab.list[[i]]
    
    #id samples with fewer than minReads
    sampReads <- data.frame(tidySeq_sampleName = rownames(curr.otuTab), numReads = rowSums(curr.otuTab))
    sampReads %>%
      filter(numReads < minReads) -> fails
    
    #annotate samples
    select.samps %>%
      filter(table == names(curr.otuTab.list)[i]) -> curr.select.samps
    fails %>%
      left_join(curr.select.samps) -> fails
    
    loss.list[[i]] <- fails
    
  }
  names(loss.list)<-names(curr.otuTab.list)
  loss.df <- list_to_df(loss.list)
  logtrts <- load_logTrt()
  loss.df %>%
    left_join(logtrts) -> loss.df
  
  #summarize by OTU table
  loss.df %>%
    group_by(gene, yrseq) %>%
    summarize(numSampsDropped = length(tidySeq_sampleName)) -> summ.loss
  summ.loss
  
  #summarize by treatment
  loss.df %>%
    group_by(gene, yrseq, species, bag) %>%
    summarize(numSampsDropped = length(tidySeq_sampleName)) %>%
    mutate(bag = ifelse(is.na(bag), "no", bag)) -> plot.df
  
  ylabname <- paste0("Number of samples with less than ", minReads, " reads")
  p <- ggplot(plot.df, aes(y = numSampsDropped, x = reorder(species, numSampsDropped), shape = bag, color = species)) +
    geom_point() + facet_grid(yrseq ~ gene) +
    scale_shape_manual(values = c(2, 16)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab(ylabname) + xlab("Wood species") + guides(color = FALSE)
  p
  
  result<-list(summ=summ.loss, df=loss.df, plot=p)
  
  return(result)
}

#------------------------------------------------------------------#
# Prep OTU tables -- requires prep_otutabs_all.R, load_otuTaxa.R, load_rotplotMetaData.R

prep_rotplotNobags_otuTabs <- function(){
  
  # requires prep_otutabs_helpers.R, load_otuTaxa.R, load_rotplotMetaData.R
  
  ## Read in the tidy taxonomy and FUNGuild tables
  table.names <- c("2012ITS","201216S","2014ITS","201416S")
  taxon.list <- read_tidy_taxon(table.names = table.names)
  
  ## Tidy the raw OTU tables
  samps <- make_miSeq_sampleList()
  
  table.names <- names(taxon.list)
  otuTab.list <- lapply(table.names, load_raw_otu) 
  names(otuTab.list) <- table.names
  otuTab.list <- batch_standarize_otu_samples(otuTab.list = otuTab.list, samps = samps)
  
  ## Subtract negative controls from OTU tables
  otuTab.list <- subtract_negativeControls(otuTab.list = otuTab.list, taxon.list = taxon.list, samps = samps)
  
  ## Trim failed samples from OTU tables
  #plot_libSize_distributions(otuTab.list = otuTab.list, select.samps = samps, maxXval = 2000, ammendfileName = "all")
  otuTab.list.t100 <- RemoveSamps.OTUtab(minReads = 100, curr.otuTab.list = otuTab.list)
  
  ## Subset OTU tables by rotplot, no bags
  samps %>%
    filter(sampleType == "rotplot") %>%
    filter(is.na(bag)) -> select.samps
  select.otuTab.list.t100 <- SelectSamps.OTUtab(select.samps = select.samps, otuTab.list = otuTab.list.t100)
  
  ## Remove unnecessary OTU columns
  tab.t100 <- lapply(select.otuTab.list.t100, remove_zeroOTUcols)
  
  ## Remove ITS OTUs that are not "Fungi" and 16S OTUs that are not "Bacteria"
  # fungi
  tab.taxa.its <- read_tidy_taxon(table.names = c("2012ITS","2014ITS"))
  taxa.its <- list_to_df(tab.taxa.its)
  taxa.its %>%
    filter(kingdom != "Fungi") # good all are fungi
  # bacteria
  tab.taxa.16s <- read_tidy_taxon(table.names = c("201216S","201416S"))
  taxa.16s <- list_to_df(tab.taxa.16s)
  taxa.16s %>%
    filter(kingdom != "Bacteria") %>%
    group_by(kingdom, source) %>%
    summarize(n = sum(!is.na(percConfid))) # there are quite a few that aren't identified as Bacteria
  #2012
  taxa.16s %>%
    filter(kingdom != "Bacteria" & source == "201216S") -> drop.otus
  tab.t100[['201216S']] <- tab.t100[['201216S']][,!colnames(tab.t100[['201216S']]) %in% drop.otus$OTUid]
  #2014
  taxa.16s %>%
    filter(kingdom != "Bacteria" & source == "201416S") -> drop.otus
  tab.t100[['201416S']] <- tab.t100[['201416S']][,!colnames(tab.t100[['201416S']]) %in% drop.otus$OTUid]
  
  return(tab.t100)
  
}

prep_taqRedo_2017ITS <- function(){
  
  # read in raw otu tables from the 2017 run
  otu.name <- "2017ITS"
  fileName.part1<-"./data/microbes/otuTables/otuTable_"
  fileName.part2<- ".tab"
  fileName<-paste0(fileName.part1, otu.name, fileName.part2)
  raw.otu<- read.delim(fileName,
                       header=TRUE,                 #make the barcodes the headers
                       stringsAsFactors = FALSE,
                       strip.white = TRUE,          #column names have trailing white space
                       row.names = 1)               #make the OTU ids the row names
  raw.otu <- as.matrix(raw.otu)
  
  # make a sample list and standardize sample names in otu table
  redo.samps <- make_taqRedo_sampleList()
  otu.redos <- standardize_otu_samples(otuTab = raw.otu, tableName = "2017ITS", samps = redo.samps)
  
  ## Subtract negative controls from OTU tables
  #NA
  
  ## Trim failed samples from OTU tables -- leave these in for now...
  #hist(rowSums(otu.redos))
  #hist(rowSums(otu.redos), breaks = 10000, xlim = c(0,500))
  #otuTab.list <- list(oturedos = otu.redos)
  #select.otuTab.list.t100 <- RemoveSamps.OTUtab(minReads = 100, curr.otuTab.list = otuTab.list)
  
  ## Remove unnecessary OTU columns
  tab.t100 <- remove_zeroOTUcols(otu.redos)
  
  ## Remove ITS OTUs that are not "Fungi" and 16S OTUs that are not "Bacteria"
  # fungi
  #tab.taxa.its <- read_tidy_taxon(table.names = c("2017ITS")) #need to add this table
  #taxa.its <- list_to_df(tab.taxa.its)
  #taxa.its %>%
  #  filter(kingdom != "Fungi")
  
  tab.t100 <- list('2017ITS' = tab.t100)
  
  return(tab.t100)
  
}
