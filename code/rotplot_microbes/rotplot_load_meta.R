
# -------------------------------------------------------------------#
# codes for plot location and log treatments

load_decoder <- function(){
  decoder<-read.delim("./data/studyMetadata/decoders/decoder.txt", stringsAsFactors = FALSE,
                      na.strings = c(" ",""))          #some elements have odd spaces so just turn these into NAs
  return(decoder)
}

load_plotLoc <- function(){
  
  decoder <- load_decoder()
  
  decoder.tmp<-decoder[,grepl("plot", colnames(decoder))] # select relevent columns
  decoder.tmp<-decoder.tmp[!is.na(decoder.tmp$plot),] # get rid of rows of NA
  colnames(decoder.tmp)<-c("plot_code","plot","wshed","topo")
  plotLoc<-decoder.tmp[,-1]
  plotLoc$plot <- as.character(plotLoc$plot)
  
  return(plotLoc)
}

load_logTrt <- function(){
  
  require('tidyverse')
  
  decoder <- load_decoder()
  
  # create a real species column (no experiment info in it)
  decoder$species_real<-unlist(strsplit(decoder$species, "2"))
  
  # create a yrdeploy column (this is the experiment (1, 2) info)
  decoder$yrdeploy<-"2009"
  decoder[decoder$numDecayYrs_2012==1,"yrdeploy"]<-"2011"
  
  # pull out log treatment levels
  logtrt.df<-decoder[,c("species_code","species","species_real","yrdeploy")]
  colnames(logtrt.df)<-c("logtrt","species","species4","yrdeploy")
  
  # there are plank species too (pull them out of the harvest.2010 df)
  harvest.2010 <- read.csv("./data/harvest/1st Harvest Data Original.csv", stringsAsFactors = FALSE)        
  
  logtrt.df %>%
    filter(yrdeploy == 2009) -> logTrt.09
  harvest.2010 %>%
    left_join(logTrt.09) %>%
    filter(is.na(species4)) -> tmp
  logtrt.planks<-data.frame(logtrt=paste(rep("plank", length(unique(tmp$species))),
                                         seq.int(from=1, to=length(unique(tmp$species))), sep="_"),
                            species=unique(tmp$species),
                            species4=unique(tmp$species),
                            yrdeploy=rep(2009, length(unique(tmp$species))))
  logTrt <- rbind(logtrt.df, logtrt.planks)
  
  # add binomials
  logBinomials <- read.csv("data/studyMetadata/logBinomials.csv", stringsAsFactors = FALSE)
  logBinomials %>%
    select(species4, family, binomial, angio.gymno) %>%
    transform(binomial = ifelse(binomial == "Vitus_vulpina", "Vitis_vulpina", binomial)) -> logBinomials
  logTrt %>%
    left_join(logBinomials) -> logTrt
  
  return(logTrt)
}

# -------------------------------------------------------------------#
# log-level information

load_harvest2010 <- function(){
  
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  harvest.2010 <- read.csv("./data/harvest/1st Harvest Data Original.csv", stringsAsFactors = FALSE)        
  
  harvest.2010 %>%
    separate(plot, into=c("wshed","topo.old"), sep="-") %>%
    mutate(topo = dplyr::recode(topo.old,
                           `LOW` = "L",
                           `HIGH` = "H")) %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    select(-topo.old) %>%
    left_join(logTrt) %>%
    mutate(bag = NA) %>%
    select(plot, symbol, logtrt, yrdeploy, yrharv, pullDate, bag) -> harvest.2010
  
  return(harvest.2010)
  
}

load_harvest2012 <- function(){
  
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  harvest.2012<-read.csv("./data/harvest/2nd Harvest Data Original AEZ 29 June 2012.csv", stringsAsFactors = FALSE)        
  
  harvest.2012 %>%
    select(Species, Symbol, Plot, Pull.date) %>%
    separate(Plot, into=c("wshed","topo"), sep="-") %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    rename(species = Species, pullDate = Pull.date, symbol = Symbol) %>%
    left_join(logTrt) %>%
    mutate(yrharv = 2012) %>%
    mutate(bag = NA) %>%
    select(plot, symbol, logtrt, yrdeploy, yrharv, pullDate, bag) -> harvest.2012
  
  return(harvest.2012)
  
}

load_harvest2014 <- function(){
  
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  harvest.2014<-read.csv("./data/harvest/2014 Harvest Spreadsheet June18b.csv", stringsAsFactors = FALSE)
  
  harvest.2014 %>%
    select(Species, Symbol, Plot, Pull.date) %>%
    separate(Plot, into=c("wshed","topo"), sep="-") %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    separate(Species, into=c("species","bag"), sep="-", fill="right") %>%
    left_join(logTrt) %>%
    rename(pullDate = Pull.date, symbol = Symbol) %>%
    mutate(yrharv = 2014) %>%
    select(plot, symbol, logtrt, yrdeploy, yrharv, pullDate, bag) -> harvest.2014
  
  return(harvest.2014)
}

load_logs <- function(){
  
  require('tidyverse')
  
  #Identify unique log samples using harvest datasheets
  harvest.2010 <- load_harvest2010()
  harvest.2012 <- load_harvest2012()
  harvest.2014 <- load_harvest2014()
  
  ## Put everything into the same dataframe
  harvest.samps <- rbind(harvest.2010, harvest.2012, harvest.2014)
  
  # standardize
  harvest.samps %>%
    mutate(pullDate = ifelse(pullDate=="", "unknown", pullDate)) %>%
    mutate(symbol = ifelse(symbol=="", "none", symbol)) -> harvest.samps
  
  # make sure each row is unique
  harvest.samps %>%
    group_by(plot, symbol, logtrt, yrdeploy, yrharv, bag) %>%
    summarize(num = length(plot),
              symb=paste(unique(symbol), collapse="_")) %>% 
  filter(num != 1) # if this is empty, then that means each row is unique
  
  # sort the rows
  harvest.samps %>%
    arrange(plot, logtrt, yrdeploy, yrharv, bag, symbol) %>%
    mutate(yrharv = as.character(yrharv)) -> logs
  
  return(logs)
}

summarize_logs <- function(logs){
  
  logTrt <- load_logTrt()
  
  #identify plank and non-plank logtrts
  sp4.planks<-logTrt[grepl("plank", logTrt$logtrt),"species4"]
  sp4.nonplanks<-logTrt[!grepl("plank", logTrt$logtrt),"species4"]
  
  #add species4
  logs %>%
    left_join(logTrt) -> logs
  
  # MAIN SAMPLES (non-plank and non-bagged)
  logs %>%
    filter(species4 %in% sp4.nonplanks) %>%
    filter(is.na(bag)) %>%
    group_by(species4, logtrt) %>%
    summarize(numSamples=length(plot),
              yrsDeployed=paste(unique(yrdeploy), collapse="_"),
              yrsHarv=paste(unique(yrharv),collapse="_"),
              numPlots=paste(unique(plot), collapse="_")) %>%
    mutate(sampType = "main") -> summ.main
  
  # PLANK SAMPLES
  logs %>%
    filter(species4 %in% sp4.planks) %>%
    filter(is.na(bag)) %>%
    group_by(species4, logtrt) %>%
    summarize(numSamples=length(plot),
              yrsDeployed=paste(unique(yrdeploy), collapse="_"),
              yrsHarv=paste(unique(yrharv),collapse="_"),
              numPlots=paste(unique(plot), collapse="_")) %>%
    mutate(sampType = "planks") -> summ.planks
  
  # BAGGED SAMPLES
  logs %>%
    filter(!is.na(bag)) %>%
    group_by(species4, logtrt) %>%
    summarize(numSamples=length(plot),
              yrsDeployed=paste(unique(yrdeploy), collapse="_"),
              yrsHarv=paste(unique(yrharv),collapse="_"),
              numPlots=paste(unique(plot), collapse="_")) %>%
    mutate(sampType = "bagged") -> summ.bagged
  
  summ <- rbind(summ.main, summ.planks, summ.bagged)
  write.csv(summ, "data_Rsynth/studyMetadata/rotplot_logSamples_summary.csv")
}

# -------------------------------------------------------------------#
# microbial sample -level information

load_2012microbs <- function(){
  
  decoder <- load_decoder()
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  
  # main study in 2012: create sampling design matrix for microbes sampled from logs
  logtrt.nonplanks <- logTrt[!grepl("plank", logTrt$logtrt),"logtrt"]
  logLoc <- decoder[!is.na(decoder$logLocation_code),"logLocation_code"]
  logmicrobs.2012.main <- expand.grid(logtrt.nonplanks, #log treatment
                                      plotLoc$plot, # plot treatment
                                      logLoc) # top and bottom of the log
  colnames(logmicrobs.2012.main)<-c("logtrt","plot","logLoc")
  logmicrobs.2012.main %>%
    mutate(rep_code = "z") %>%
    mutate(yrharv = 2012) %>%
    mutate(bag = NA) %>%
    mutate(sampleName = paste0(logtrt, plot, logLoc)) %>%
    select(logtrt, plot, logLoc, bag, rep_code, yrharv, sampleName) -> logmicrobs.2012
  
  return(logmicrobs.2012)
  
}

load_2014microbs.main <- function(){
  
  decoder <- load_decoder()
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  
  # main study in 2014
  logmicrobs.2014.main <- read.csv("./data/studyMetadata/sampleLists/2014SampleList_main.csv", stringsAsFactors = FALSE)
  logmicrobs.2014.main %>%
    select(1:11) %>%
    filter(!is.na(Sp) & Sp != "") %>%
    gather(plot, sampleName, -c(Sp, Yrs, Side)) %>%
    filter(!is.na(sampleName)) %>%
    separate(plot, into=c("X","wshed","topo"), sep=c(1,2)) %>%
    select(-X) %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    rename("species4" = "Sp") %>%
    mutate(species4 = dplyr::recode(species4, `PRSE` = "PRSE/PRVI")) %>% # need to replace PRSE with PRSE/PRVI so that logTrt df matches
    mutate(yrdeploy = dplyr::recode(Yrs, `5` = "2009", `3` = "2011")) %>% # need to add this so that logTrt df matches
    left_join(logTrt) -> logmicrobs.2014.main.pause
  
  # why didn't these get assigned a logtrt?
  logmicrobs.2014.main.pause %>%
    filter(is.na(logtrt)) %>%
    arrange(species, Yrs) -> tmp
  tmp
  
  # this looks like a data entry mistake....
  # 2014SampleList_main indicates that QUAL samples had been decaying for 5 yrs but...
  logTrt[logTrt$species=="QUAL",] #log treatment table indicates that QUAL was only deployed in 2011 (3 yrs)
  
  # change 2014SampleList_main to reflect logTrt
  logmicrobs.2014.main.pause %>%
    mutate(Yrs = ifelse(species4 == "QUAL", 3, Yrs)) %>%
    mutate(yrdeploy = dplyr::recode(Yrs, `5` = "2009", `3` = "2011")) %>% # update for QUAL
    select(-c(species, family, binomial, angio.gymno, logtrt)) -> new.2014SampleList_main
  
  new.2014SampleList_main %>%
    left_join(logTrt) %>% # try this again, now it works
    separate(sampleName, into=c("first","second","logLoc"), sep=c(1,2), remove=FALSE) %>%
    mutate(logLoc = ifelse(logLoc == "_", "mush", logLoc)) %>%
    mutate(yrharv = "2014") %>%
    mutate(rep_code = NA) %>%
    mutate(bag = NA) %>%
    select(logtrt, plot, logLoc, bag, rep_code, yrharv, sampleName) -> logmicrobs.2014.main
  
  return(logmicrobs.2014.main)
  
}

load_2014microbes.bagged <- function(){
  
  decoder <- load_decoder()
  plotLoc <- load_plotLoc()
  logTrt <- load_logTrt()
  
  bag.2014<-read.csv("./data/studyMetadata/sampleLists/2014SampleList_bagged.csv", stringsAsFactors = FALSE)
  
  bag.2014 %>%
    gather(plot, sampleName, -c(Sp, Yrs, Side)) %>%
    filter(!is.na(sampleName)) %>%
    separate(plot, into=c("X","wshed","topo"), sep=c(1,2)) %>%
    mutate(wshed = as.integer(wshed)) %>%
    left_join(plotLoc) %>%
    rename("species4" = "Sp") %>%
    mutate(yrdeploy = dplyr::recode(Yrs, `5` = "2009", `3` = "2011")) %>% # need to add this so that logTrt df matches
    left_join(logTrt) %>%
    separate(sampleName, into=c("first","second"), sep="-", remove=FALSE) %>% #get rid of the "bag" part
    select(-c(X, second)) %>%
    separate(first, into=c("first","second","logLoc"), sep=c(1,2)) -> pause #split the remaining sample id to capture logLoc (this should be the same as side)
    
  # check that Side matches logLoc... and it does
  sum(pause$Side != pause$logLoc) 
  
  pause %>%
    select(-c(Side, first, second)) %>%
    mutate(logLoc = ifelse(logLoc == "_", "mush", logLoc)) %>%
    mutate(yrharv = "2014") %>%
    mutate(rep_code = NA) %>%
    mutate(bag = "b") %>%
    select(logtrt, plot, logLoc, bag, rep_code, yrharv, sampleName) -> bag.2014
  
  return(bag.2014)
}

load_microbes <- function(){
  
  require('tidyverse')
  
  logmicrobs.2012 <- load_2012microbs() # 2012 main study, df based on experimental design matrix
  logmicrobs.2014.main <- load_2014microbs.main() # 2014 main study, df based on sample list
  logmicrobs.2014.bagged <- load_2014microbes.bagged() # 2014 bagged logs, df based on sample list
  logmicrobs <- rbind(logmicrobs.2012, logmicrobs.2014.main, logmicrobs.2014.bagged)
  
  # make sure each row is unique
  logmicrobs %>%
    group_by(plot, logtrt, logLoc, bag, yrharv) %>%
    summarize(num=length(plot),
              sampName=paste(unique(sampleName), collapse="_")) %>%
    filter(num !=1) # good they are all unique
  
  #sort and add unique identifier
  logmicrobs %>%
    arrange(plot, logtrt, logLoc, yrharv, bag, sampleName) -> logmicrobs
  
  return(logmicrobs)
}

summarize_logmicrobs <- function(logmicrobs){
  
  logTrt <- load_logTrt()
  
  #identify plank and non-plank logtrts
  sp4.planks<-logTrt[grepl("plank", logTrt$logtrt),"species4"]
  sp4.nonplanks<-logTrt[!grepl("plank", logTrt$logtrt),"species4"]
  
  #add species4
  logmicrobs %>%
    left_join(logTrt) -> logmicrobs
  
  # MAIN SAMPLES (non-plank and non-bagged)
  logmicrobs %>%
    filter(species4 %in% sp4.nonplanks) %>%
    filter(is.na(bag)) %>%
    group_by(species4, logtrt) %>%
    summarize(numSamples=length(plot),
              yrsDeployed=paste(unique(yrdeploy), collapse="_"),
              yrsHarv=paste(unique(yrharv),collapse="_"),
              Plots=paste(unique(plot), collapse="_"),
              LogLocs=paste(unique(logLoc), collapse="_")) %>%
    mutate(sampType = "main") -> summ.main
  
  # # PLANK SAMPLES
  # logmicrobs %>%
  #   filter(species4 %in% sp4.planks) %>%
  #   filter(is.na(bag)) %>%
  #   group_by(species4, logtrt) %>%
  #   summarize(numSamples=length(plot),
  #             yrsDeployed=paste(unique(yrdeploy), collapse="_"),
  #             yrsHarv=paste(unique(yrharv),collapse="_"),
  #             numPlots=paste(unique(plot), collapse="_"),
  #             LogLocs=paste(unique(logLoc), collapse="_")) %>%
  #   mutate(sampType = "planks") -> summ.planks
  
  # BAGGED SAMPLES
  logmicrobs %>%
    filter(!is.na(bag)) %>%
    group_by(species4, logtrt) %>%
    summarize(numSamples=length(plot),
              yrsDeployed=paste(unique(yrdeploy), collapse="_"),
              yrsHarv=paste(unique(yrharv),collapse="_"),
              numPlots=paste(unique(plot), collapse="_"),
              LogLocs=paste(unique(logLoc), collapse="_")) %>%
    mutate(sampType = "bagged") -> summ.bagged
  
  summ <- rbind(summ.main, summ.bagged)
  write.csv(summ, "data_Rsynth/studyMetadata/rotplot_microbSamples_summary.csv")
  
  
}


# -------------------------------------------------------------------#
# make data subsets 

# function to identitfy key fungal groups (sapros and basidios) in the taxa tables
id_group <- function(tab.taxa, group){
  fung.tab.names <- c('2012ITS', '2014ITS')
  tab.taxa.list <- list()
  for(i in 1:length(fung.tab.names)){
    fung.taxa <- tab.taxa[[fung.tab.names[i]]]
    if(group == "sapros"){
      sapro.levels <- unique(fung.taxa$Trophic.Mode)[grepl("Sapro", unique(fung.taxa$Trophic.Mode))]
      fung.taxa %>%
        filter(Trophic.Mode %in% sapro.levels) %>%
        filter(Confidence.Ranking %in% c("Probable", "Highly Probable")) -> group.result
    }
    if(group == "basidios"){
      fung.taxa %>%
        filter(phylum == "Basidiomycota") -> group.result
    }
    tab.taxa.list[[i]]<- group.result
  }
  names(tab.taxa.list) <- fung.tab.names
  return(tab.taxa.list)
}


#### testing
# overlap.otuPreplist <- lapply(tab.list.pa, function(x){
#tmp <- make_overlap_subsets(tab.otu = x, samp.meta = samp.meta, add_group_tf = T, tab.taxa = tab.taxa)
# })
#tab.otu <- tab.list.pa[[2]]
#samp.meta = samp.meta
#add_group_tf = T
#tab.taxa = tab.taxa

# add group fxn
#group = "sapros"
add_group <- function(tab.dataset, tab.taxa, group){
  
  tab.taxa.group <- id_group(tab.taxa = tab.taxa, group = group)  
  tab.dataset.fung <- tab.dataset[names(tab.dataset) %in% names(tab.taxa.group)]
  tab.dataset.fung.group <- list()
  
  #loop through the otu tabs for different years (2012, 2014)
  for(i in 1:length(tab.dataset.fung)){
    
    # id the current otu and meta tables
    curr.otu <- tab.dataset.fung[[i]]$otu # 1158 OTUs (all fungi from the overlap samples)
    curr.meta <- tab.dataset.fung[[i]]$meta
    
    # id the current sapro OTUs
    curr.group.otus <- tab.taxa.group[[i]] # 523 OTUs (only sapros from all the samples)
    condition <- colnames(curr.otu) %in% curr.group.otus$OTUid # of all the otus from overlap samples, which are sapros?
    sum(condition) # 309 of 1158 are sapros
    curr.otu.onlygroup <- curr.otu[,condition]
    
    # check for empty rows
    emptyrows <- rowSums(curr.otu.onlygroup) == 0
    sum(emptyrows) > 0
    if(sum(emptyrows) > 0){
      # remove the empty row
      curr.otu.onlygroup <- curr.otu.onlygroup[!emptyrows,]
      dim(curr.otu.onlygroup)
      
      # check for empty cols
      emptycols <- colSums(curr.otu.onlygroup) == 0
      sum(emptycols) != 0
      if(sum(emptycols) != 0){
        # remove the empty cols
        curr.otu.onlygroup <- curr.otu.onlygroup[,!emptycols]
      }
      
      # make the new otu table match curr.meta
      o <- match(row.names(curr.otu.onlygroup), curr.meta$tidySeq_sampleName)
      curr.meta <- curr.meta[o,]
    }
    
    # update the otu table
    tab.dataset.fung.group[[i]] <- list(otu = curr.otu.onlygroup, 
                                        meta = curr.meta)
    
  }
  new.names <- paste0(names(tab.dataset.fung),"_",group)
  names(tab.dataset.fung.group) <- new.names
  
  return(tab.dataset.fung.group)
}

# function to compile data subsets w/ groups
#samp.data = samp.overlap
compile_datasubsets_wgroups <- function(tab.otu, samp.data, tab.taxa){
  
  # define the yrharv, gene, and otu table
  yrseq.vec <- substring(names(tab.otu), 1, 4)
  gene.vec <- substring(names(tab.otu), 5, 8)
  
  # loop through each otu table
  tab.dataset <- list()
  for(i in 1:length(tab.otu)){
    
    otu.mat <- tab.otu[[i]]
    
    # make a sample list
    samp.data %>%
      filter(yrseq == yrseq.vec[i]) %>%
      filter(gene == gene.vec[i]) %>%
      filter(tidySeq_sampleName %in% row.names(otu.mat)) %>%
      select(-c(seq_sampleName, sampleType, bag)) -> samp.select
    
    # subset the otu mat and reorder the rows so they match samp.select
    otu.mat.select <- otu.mat[row.names(otu.mat) %in% samp.select$tidySeq_sampleName,]
    o <- match(samp.select$tidySeq_sampleName, row.names(otu.mat.select))
    otu.mat.select.o <- otu.mat.select[o,]
    
    tab.dataset[[i]] <- list(otu = otu.mat.select.o, meta = samp.select)
    
  }
  names(tab.dataset) <- names(tab.otu)
  lapply(tab.dataset, function(x){dim(x$otu)})
  
  # make sapros
  tab.dataset.fung.sapros <- add_group(tab.dataset = tab.dataset, 
                                       tab.taxa = tab.taxa, 
                                       group = "sapros")
  lapply(tab.dataset.fung.sapros, function(x){dim(x$otu)})
  
  # make basidios
  tab.dataset.fung.basidios <- add_group(tab.dataset = tab.dataset, 
                                         tab.taxa = tab.taxa, group = "basidios")
  lapply(tab.dataset.fung.basidios, function(x){dim(x$otu)})
  
  # combine into data object
  names1 <- names(tab.dataset)
  names2 <- names(tab.dataset.fung.sapros)
  names3 <- names(tab.dataset.fung.basidios)
  
  tab.dataset <- list(tab.dataset[[names1[1]]], 
                      tab.dataset[[names1[2]]],
                      tab.dataset[[names1[3]]],
                      tab.dataset[[names1[4]]],
                      
                      tab.dataset.fung.sapros[[names2[1]]], 
                      tab.dataset.fung.sapros[[names2[2]]],
                      
                      tab.dataset.fung.basidios[[names3[1]]], 
                      tab.dataset.fung.basidios[[names3[2]]]
  )
  
  names(tab.dataset) <- c(names1, names2, names3)
  #lapply(tab.dataset, function(x){dim(x$otu)})
  #lapply(tab.dataset, function(x){dim(x$meta)})
  
  return(tab.dataset)
}

compile_datasubsets <- function(tab.otu, samp.data){
  
  # define the yrharv, gene, and otu table
  yrseq.vec <- substring(names(tab.otu), 1, 4)
  gene.vec <- substring(names(tab.otu), 5, 8)

  # loop through each otu table
  tab.dataset <- list()
  for(i in 1:length(tab.otu)){
    
    otu.mat <- tab.otu[[i]]
    
    # make a sample list
    samp.data %>%
      filter(yrseq == yrseq.vec[i]) %>%
      filter(gene == gene.vec[i]) %>%
      filter(tidySeq_sampleName %in% row.names(otu.mat)) %>%
      select(-c(seq_sampleName, sampleType, bag)) -> samp.select
    
    # subset the otu mat and reorder the rows so they match samp.select
    otu.mat.select <- otu.mat[row.names(otu.mat) %in% samp.select$tidySeq_sampleName,]
    o <- match(samp.select$tidySeq_sampleName, row.names(otu.mat.select))
    otu.mat.select.o <- otu.mat.select[o,]

    tab.dataset[[i]] <- list(otu = otu.mat.select.o, meta = samp.select)
    
  }
  names(tab.dataset) <- names(tab.otu)
  lapply(tab.dataset, function(x){dim(x$otu)})
  
  return(tab.dataset)
  
}

# this fxn uses plottingParams_rotplotNobags()
make_overlap_subsets <- function(tab.otu, samp.meta, add_group_tf, tab.taxa){
  
  ### Overlapping wood species in 2012, 2014
  pars <- plottingParams_rotplotNobags()
  samp.meta %>%
    filter(species4 %in% pars$overlapSp$levels) -> samp.overlap
  if(add_group_tf == T){
    data.overlap <- compile_datasubsets_wgroups(tab.otu = tab.otu, samp.data = samp.overlap, 
                                                tab.taxa = tab.taxa) # fxn in rotplot_load_meta.R
    data.overlap <- list(`2012ITS_o` = data.overlap[['2012ITS']], `2014ITS_o` = data.overlap[['2014ITS']],
                         `201216S_o` = data.overlap[['201216S']], `201416S_o` = data.overlap[['201416S']],
                         `2012ITS_sapros_o` = data.overlap[['2012ITS_sapros']], `2014ITS_sapros_o` = data.overlap[['2014ITS_sapros']],
                         `2012ITS_basidios_o` = data.overlap[['2012ITS_basidios']], `2014ITS_basidios_o` = data.overlap[['2014ITS_basidios']])
  }else{
    data.overlap <- compile_datasubsets(tab.otu = tab.otu, samp.data = samp.overlap) # fxn in rotplot_load_meta.R
    data.overlap <- list(`2012ITS_o` = data.overlap[['2012ITS']], `2014ITS_o` = data.overlap[['2014ITS']],
                         `201216S_o` = data.overlap[['201216S']], `201416S_o` = data.overlap[['201416S']])
  }
  
  
  return(data.overlap)
  
}

make_deploy_subsets <- function(tab.otu, samp.meta, add_group_tf, tab.taxa){
  
  if(add_group_tf == T){
    samp.meta %>%
      filter(yrdeploy == 2009) -> samp.deploy09
    data.deploy09 <- compile_datasubsets_wgroups(tab.otu = tab.otu, samp.data = samp.deploy09, tab.taxa = tab.taxa)
    samp.meta %>%
      filter(yrdeploy == 2011) -> samp.deploy11
    data.deploy11 <- compile_datasubsets_wgroups(tab.otu = tab.otu, samp.data = samp.deploy11, tab.taxa = tab.taxa)
    data.deploy <- list(`2012ITS_d11` = data.deploy11[['2012ITS']], 
                        `2012ITS_d09` = data.deploy09[['2012ITS']], 
                        `2014ITS_d11` = data.deploy11[['2014ITS']], 
                        `2014ITS_d09` = data.deploy09[['2014ITS']],
                        
                        `201216S_d11` = data.deploy11[['201216S']], 
                        `201216S_d09` = data.deploy09[['201216S']], 
                        `201416S_d11` = data.deploy11[['201416S']], 
                        `201416S_d09` = data.deploy09[['201416S']],
                        
                        `2012ITS_d11_sapros` = data.deploy11[['2012ITS_sapros']], 
                        `2012ITS_d09_sapros` = data.deploy09[['2012ITS_sapros']], 
                        `2014ITS_d11_sapros` = data.deploy11[['2014ITS_sapros']], 
                        `2014ITS_d09_sapros` = data.deploy09[['2014ITS_sapros']],
                        
                        `2012ITS_d11_basidios` = data.deploy11[['2012ITS_basidios']], 
                        `2012ITS_d09_basidios` = data.deploy09[['2012ITS_basidios']], 
                        `2014ITS_d11_basidios` = data.deploy11[['2014ITS_basidios']], 
                        `2014ITS_d09_basidios` = data.deploy09[['2014ITS_basidios']])
    
  }else{
    samp.meta %>%
      filter(yrdeploy == 2009) -> samp.deploy09
    data.deploy09 <- compile_datasubsets(tab.otu = tab.otu, samp.data = samp.deploy09)
    samp.meta %>%
      filter(yrdeploy == 2011) -> samp.deploy11
    data.deploy11 <- compile_datasubsets(tab.otu = tab.otu, samp.data = samp.deploy11)
    data.deploy <- list(`2012ITS_d11` = data.deploy11[['2012ITS']], 
                        `2012ITS_d09` = data.deploy09[['2012ITS']], 
                        `2014ITS_d11` = data.deploy11[['2014ITS']], 
                        `2014ITS_d09` = data.deploy09[['2014ITS']],
                        
                        `201216S_d11` = data.deploy11[['201216S']], 
                        `201216S_d09` = data.deploy09[['201216S']], 
                        `201416S_d11` = data.deploy11[['201416S']], 
                        `201416S_d09` = data.deploy09[['201416S']])
  }
  
  
  return(data.deploy)
  
}


# -------------------------------------------------------------------#
# identify missing microbial sequence samples

identify_missing_rotplotsamps <- function(samples){
  
  # create a full experimental matrix
  exp.mat <- expand.grid(gene = c("16S","ITS"),
                         topo = c("H","L"),
                         logLoc = c("t","b"),
                         dep.harv.sp = unique(samples$dep.harv.sp),
                         wshed = c(1,2,3,4))
  
  # remove mush samples
  samples %>%
    filter(logLoc != "mush") -> samples.nomush
  # identify samples missing from the full experimental matrix
  exp.mat.missing <- exp.mat %>%
    left_join(samples.nomush) %>%
    filter(is.na(table)) %>%
    select(gene, topo, logLoc, dep.harv.sp, wshed)
  # which of these match with "mush" samples?
  samples %>%
    filter(logLoc == "mush") %>%
    select(gene, topo, dep.harv.sp, wshed) %>%
    mutate(mushfail = TRUE) -> samples.mush
  exp.mat.missing %>%
    left_join(samples.mush) %>%
    mutate(mushfail = ifelse(is.na(mushfail), FALSE, TRUE)) -> exp.mat.missing
  # which of these differ between genes?
  exp.mat.missing %>%
    filter(mushfail == FALSE) %>%
    group_by(topo, logLoc, dep.harv.sp, wshed) %>%
    summarize(n = length(gene)) %>%
    filter(n == 1) %>%
    select(topo, logLoc, dep.harv.sp, wshed) %>%
    mutate(seqfail = TRUE) -> samples.seqfail
  exp.mat.missing %>%
    left_join(samples.seqfail) %>%
    mutate(seqfail = ifelse(is.na(seqfail), FALSE, TRUE)) %>%
    separate(dep.harv.sp, into = c("yrdeploy","yrharv","species4"), sep = "_", remove = F) -> exp.mat.missing
  
  return(exp.mat.missing)
}

# -------------------------------------------------------------------#
# make a set of taq test objects : list of otu tables, dataframe of sample ids

make_taqtest_objs <- function(samp.taqredos, samp.meta, tab.t100, taqredo.t100){
  
  # subset the same samples from different sequencing runs
  samp.taqredos %>%
    select(gene, plot, logtrt, logLoc) %>%
    group_by(gene, plot, logtrt, logLoc) %>%
    summarize(n = length(gene)) -> summ.samp.indx
  # 4 samples per samptype (plot+logtrt+logLoc) because yrharv x taqType
  
  # find these samptypes in the 2012 and 2014 otu tables
  samp.meta %>%
    left_join(summ.samp.indx) %>%
    filter(!is.na(n)) %>%
    group_by(gene, plot, logtrt, logLoc) %>%
    summarize(n2 = length(gene)) -> tmp
  # 2 samples per samptype because yrseq = 2012, 2014
  samp.meta %>%
    left_join(summ.samp.indx) %>%
    filter(!is.na(n)) -> select.samps
  tab.select <- SelectSamps.OTUtab(select.samps = select.samps, otuTab.list = tab.t100)
  
  # make a list of otu tables
  tab.taqcomp <- list(tab.select[['2012ITS']],
                      tab.select[['2014ITS']],
                      taqredo.t100[['2017ITS']])
  names(tab.taqcomp) <- c('2012ITS','2014ITS','2017ITS')
  
  # make a sample index
  select.samps %>%
    select(-c(species, species4, yrdeploy, family, binomial, angio.gymno, wshed, topo, n)) %>%
    mutate(taqType =  ifelse(yrseq == "2012", "GoTaq","Q5")) -> tmp
  samp.indx <- rbind(tmp, samp.taqredos)
  samp.indx %>%
    mutate(sampid = paste(plot, logtrt, logLoc, yrharv, sep ="_")) %>%
    mutate(taqType.yrseq = paste(taqType, yrseq, sep = "_")) -> samp.indx
  
  result <- list(otulist = tab.taqcomp, sampdf = samp.indx)
  return(result)
}
