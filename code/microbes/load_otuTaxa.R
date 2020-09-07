
#------------------------------------------------------------------#
# Tidy the taxonomy and FUNGuild tables

TidyFungTax <- function(tableName){
  
  path <- "./data_Rsynth/FUNGtables/FUNGtables_"
  
  #read in FUNGuild-formated otu file and pick off the OTUid and taxonomy columns
  fileName<-paste(path, tableName,".tab", sep="")
  curr.otut<-read.delim(fileName, header=TRUE, stringsAsFactors = FALSE, row.names = 1)
  curr.taxon<-data.frame(OTUid = row.names(curr.otut), taxonomy = curr.otut[,"taxonomy"])
  
  #simplify the taxonomy formating
  curr.taxon %>%
    separate(col=taxonomy, into=c("percConfid","species","EU","SH","reps","taxon_hier"), sep = "\\|") %>%
    select(OTUid, percConfid, species, taxon_hier) %>%
    rename("genusSpecies"="species") %>%
    separate(taxon_hier, into = c("k","p","c","o","f","g","s"), sep = ";") -> pause
  #get rid of hier codes in front of data in each cell
  pause %>%
    select(-c(OTUid, percConfid, genusSpecies)) -> hier.expand.df
  hier<-c("kingdom","phylum","class","order","family","genus","species")
  hier.letter <- substring(hier, 1, 1)
  hier.code<-paste(hier.letter, "__", sep="")
  for(t in 1:length(colnames(hier.expand.df))){
    hier.expand.df %>%
      separate_(col=hier.letter[t], into=c("drop", hier[t]), sep=hier.code[t]) %>%
      select(-drop) -> hier.expand.df
  }
  curr.taxon <- data.frame(pause[,c("OTUid","percConfid","genusSpecies")], hier.expand.df)
  curr.taxon %>%
    separate(col="species", into=c("species","drop"), sep=" ") %>%
    select(-drop) -> curr.taxon
  
  return(curr.taxon)
  
}

AddFUNGuild <- function(tidyTab, tableName){
  
  path <- "./data/microbes/matchedFUNGuild/FUNGtables_"
  
  #read in matched FUNGuild table and slice out the OTU abundance columns
  #fileName<-paste(path, tableName,".guilds_matched.txt", sep="")
  fileName<-paste(path, tableName,".guilds.txt", sep="")
  curr.fungt<-read.delim(fileName, sep="\t", stringsAsFactors = FALSE, header=FALSE, row.names=1) 
  colnames(curr.fungt)<-curr.fungt[1,] #fix format
  curr.fungt<-curr.fungt[-1,]
  colNum<-which(colnames(curr.fungt) == 'taxonomy')
  curr.fungt.simp<-data.frame(OTUid=row.names(curr.fungt), curr.fungt[,c(colNum:length(colnames(curr.fungt)))]) 
  
  #fix format issue with the column names
  # if(tableName != "2017ITS_CWD"){
  #   curr.fungt.simp %>%
  #     rename("duplicate"="Notes",
  #            "Notes"="Citation.Source",
  #            "Citation.Source"="Var.12") %>%
  #     select(-duplicate) -> curr.fungt.simp
  # }
  
  #add FUNGuild info to the tidy taxon file if there isn't FUNGuild info there already
  #sum(colnames(tidyTab) %in% "Trophic.Mode") == 0
  if(sum(colnames(tidyTab) %in% "Trophic.Mode") == 0){ # if this is 0, then there are no columns named Trophic.Mode
    tidyTab %>%
      left_join(curr.fungt.simp) -> tidyTab
    funguildCols <- colnames(curr.fungt.simp)[-c(1,2)]
    tidyTab[tidyTab$Taxon == "-", funguildCols] <- NA
  }
  
  return(tidyTab)
}

annotate_fungTax <- function(table.names){
  
  # read in the taxonomy tables and tidy them up
  taxon.indx.list.fung <- list()
  for(i in 1:length(table.names)){
    taxon.indx.list.fung[[i]] <- TidyFungTax(tableName = table.names[i])
  }
  names(taxon.indx.list.fung) <- table.names
  
  # FUNGuild data and add it to the tidy taxon tables
  taxon.indx.list.fungG <- list()
  for(i in 1:length(table.names)){
    taxon.indx.list.fungG[[i]] <- AddFUNGuild(tidyTab = taxon.indx.list.fung[[i]], tableName = table.names[i])
  }
  names(taxon.indx.list.fungG) <- table.names
  
  return(taxon.indx.list.fungG)
}

TidyBactTax <- function(tableName){
  
  path <- "./data/microbes/classifiedRDP/RDP_"
  
  #read in RDP output file and pick off the OTUid and taxonomy columns
  fileName <- paste(path, tableName,".txt", sep="")
  rdp.table <- read.delim(fileName, header=FALSE, stringsAsFactors = FALSE)
  
  #these next steps are really similar to the function in'formatRDPforFUNGuild.Rmd'
  
  #a. simplify taxonomy
  df.tmp <- rdp.table[, c("V1","V3","V5","V6","V8","V9","V11","V12","V14","V15","V17","V18","V20")]
  colnames(df.tmp) <- c("OTUId","kingdom","k_perc","phylum","p_perc","class","c_perc","order","o_perc","family","f_perc","genus","g_perc")
  
  #strip out taxa identifiers with confidence estimate less than 50%
  for(i in 1:dim(df.tmp)[1]){
    curr.otu <- df.tmp[i,]
    
    #curr.otu$g_perc < .5 # if g perc is less than 50%, label it as unidentified
    if(curr.otu$g_perc < .5){
      curr.otu$genus <- "unidentified"
      
      #curr.otu$f_perc < .5 # if f perc is less than 50%, label it as unidentified
      if(curr.otu$f_perc < .5){
        curr.otu$family <- "unidentified"
        
        #curr.otu$o_perc < .5 # if o perc is less than 50%, label it as unidentified
        if(curr.otu$o_perc < .5 ){
          curr.otu$order <- "unidentified"
          
          #curr.otu$c_perc < .5 # if c perc is less than 50%, label it as unidentified
          if(curr.otu$c_perc < .5){
            curr.otu$class <- "unidentified"
            
            # curr.otu$p_perc < .5 # if p perc is less than 50%, label it as unidentified
            if(curr.otu$p_perc < .5){
              curr.otu$phylum <- "unidentified"
            }
          }
        }
      }
    }
    
    df.tmp[i,] <- curr.otu
  }
  
  #b. add missing columns
  df.tmp$species <- "sp" 
  
  #c. replace columns that include "unidentified with just unidentified"
  COLS <- c("kingdom","phylum","class","order","family","genus")
  for (s in 1:length(COLS)){
    #pull out the current column
    columnThing <- df.tmp[,COLS[s]]
    #find "unidentified"s and replace cell with simple "unidentified"
    columnThing[grepl("unidentified", columnThing)] <- "unidentified"
    #update the column in the original dataframe
    df.tmp[,COLS[s]] <- columnThing
  }
  
  #here's where this starts to differ from the ITS fxn referenced above...
  
  #d. pull stuff together to make the simplified columns
  perc <- df.tmp$g_perc
  kingdom <- df.tmp$kingdom
  phylum <- df.tmp$phylum
  class <- df.tmp$class
  order <- df.tmp$order
  family <- df.tmp$family
  genus <- df.tmp$genus
  species <- df.tmp$species
  # taxon_hier <- paste(paste("k__",kingdom, sep =""),
  #                     paste("p__",phylum, sep =""),
  #                     paste("c__",class, sep =""),
  #                     paste("o__",order, sep =""),
  #                     paste("f__",family, sep =""),
  #                     paste("g__",genus, sep =""),
  #                     paste("s__",species, " ...", sep =""), sep = ";")
  # curr.taxon.indx<-data.frame(OTUid=df.tmp$OTUId,
  #                               percConfid=perc,
  #                               genusSpecies=paste(df.tmp$genus, df.tmp$species, sep="_"),
  #                               taxon_hier=taxon_hier)
  
  # pull taxon_hier apart like for fungal taxonomy tables
  curr.taxon.indx <- data.frame(OTUid=df.tmp$OTUId,
                                percConfid=perc,
                                genusSpecies=paste(df.tmp$genus, df.tmp$species, sep="_"),
                                kingdom=kingdom,
                                phylum=phylum,
                                class=class,
                                order=order,
                                family=family, 
                                genus=genus,
                                species=species)
  
  return(curr.taxon.indx)
  
}

annotate_bactTax <- function(table.names){
  
  # read in the taxonomy tables and tidy them up
  taxon.indx.list.bact <- list()
  for(i in 1:length(table.names)){
    taxon.indx.list.bact[[i]] <- TidyBactTax(tableName = table.names[i])
  }
  names(taxon.indx.list.bact) <- table.names
  
  return(taxon.indx.list.bact)
}

write_tidy_taxon<-function(tidyTab, tableName){
  
  path <- "./data_Rsynth/taxonTables_tidy/taxonTable_"
  fileName <- paste(path, tableName,".csv", sep="")
  write.csv(tidyTab, file = fileName, row.names = FALSE)
  
}

read_tidy_taxon <- function(table.names){
  
  path <- "./data_Rsynth/taxonTables_tidy"
  files <- paste0(path, "/taxonTable_",table.names, ".csv")
  
  taxon.list <- list()
  for(i in 1:length(files)){
    taxon.list[[i]] <- read.csv(file = files[i], stringsAsFactors = F)
  }
  names(taxon.list) <- table.names
  
  return(taxon.list)
  
}

trim_tab.taxa <- function(tab.otus, tab.taxa){
  
  # otu prep may remove OTUIds, so need to also remove them from tab.taxa
  check <- sum(names(tab.otus) != names(tab.taxa)) # make sure the df names match
  if(check == 0){
    for(i in 1:length(tab.otus)){
      
      curr.otu <- tab.otus[[i]]
      curr.taxa <- tab.taxa[[i]]
      new.taxa <- curr.taxa[match(colnames(curr.otu), curr.taxa$OTUid),]
      
      if(dim(curr.otu)[2] != dim(new.taxa)[1]){
        print("number of unique OTUs don't match")
      }
      
      tab.taxa[[i]] <- new.taxa
    }
  }
  
  return(tab.taxa)
}

