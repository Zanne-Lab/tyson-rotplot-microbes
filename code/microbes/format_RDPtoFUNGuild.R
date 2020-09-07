# format_RDPtoFUNGuild.R

RDPtoFUNG_ITS <- function(otu.table, rdp.table){
  
  #1. reformat rdp.table
  
  #a. simplify taxonomy
  df.tmp <- rdp.table[, c("V1","V3","V5","V6","V8","V9","V11","V12","V14","V15","V17","V18","V20","V21","V23")]
  colnames(df.tmp) <- c("OTUId","kingdom","k_perc","phylum","p_perc","class","c_perc","order","o_perc","family","f_perc","genus","g_perc","genus_species_code","s_perc")
  
  #b. split the column with multiple types of info
  df.tmp1 <- separate(data = df.tmp, col = genus_species_code, into = c("genus_species", "shcode"), sep = "\\|")
  df.tmp2 <- separate(data = df.tmp1, col = genus_species, into = c("Genus","species"), sep = "_", extra = "merge")
  #get rid of Genus
  df.tmp <- df.tmp2[,!colnames(df.tmp2) %in% c("Genus")]
  
  #c. strip out taxa identifiers with confidence estimate less than 50% - starting from species and going higher
  for(i in 1:dim(df.tmp)[1]){
    curr.otu <- df.tmp[i,]
  
    #curr.otu$s_perc < .5 # if s perc is less than 50%, label it as unidentified
    if(curr.otu$s_perc < .5){ 
      curr.otu$species <- "sp"
      
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
    }
    
    df.tmp[i,] <- curr.otu
  }
  
  #d. replace columns that include "unidentified with just unidentified"
  i<-0
  COLS <- c("kingdom","phylum","class","order","family","genus","species")
  for (i in 1:length(COLS)){
    #pull out the current column
    columnThing <- df.tmp[,COLS[i]]
    #find "unidentified"s and replace cell with simple "unidentified"
    columnThing[grepl("unidentified", columnThing)] <- "unidentified"
    #update the column in the original dataframe
    df.tmp[,COLS[i]] <- columnThing
  }
  
  #e. make a genus_species column
  # create genus_species
  df.tmp$genus_species <- paste(df.tmp$genus, df.tmp$species, sep="_")
  df.tmp[df.tmp$genus_species=="unidentified_sp","genus_species"]<-"unidentified"
  df.tmp3<-df.tmp
  
  #d. pull stuff together and make a 'taxonomy' column
  perc <- df.tmp3$s_perc
  genus_species <- df.tmp3$genus_species
  eucode <- "EU"
  shcode <- df.tmp3$shcode
  repinfo <- "reps"
  kingdom <- df.tmp3$kingdom
  phylum <- df.tmp3$phylum
  class <- df.tmp3$class
  order <- df.tmp3$order
  family <- df.tmp3$family
  genus <- df.tmp3$genus
  species <- df.tmp3$species
  df.tmp3[,"taxonomy"]<-paste(perc, 
                              genus_species, eucode, shcode, repinfo,
                              paste(paste("k__",kingdom, sep =""),
                                    paste("p__",phylum, sep =""),
                                    paste("c__",class, sep =""),
                                    paste("o__",order, sep =""),
                                    paste("f__",family, sep =""),
                                    paste("g__",genus, sep =""),
                                    paste("s__",species, " ...", sep =""), sep = ";"), 
                              sep = "|")
  
  #2. attach the 'taxonomy' column to the otu table
  ind.tax <- df.tmp3[,c("OTUId","taxonomy")]
  fung.table <- merge(otu.table, ind.tax)
  
  return(fung.table)
}

#need to update to match ITS?
RDPtoFUNG_LSU <- function(otu.table, rdp.table){
  
  #1. reformat rdp.table
  
  #a. simplify taxonomy
  df.tmp <- rdp.table[, c("V1","V3","V5","V6","V8","V9","V11","V12","V14","V15","V17","V18","V20")]
  colnames(df.tmp) <- c("OTUId","kingdom","k_perc","phylum","p_perc","class","c_perc","order","o_perc","family","f_perc","genus","g_perc")
  
  #b. add missing columns
  df.tmp$species <- "sp"
  
  #strip out taxa identifiers with confidence estimate less than 50%
  df.tmp[df.tmp$k_perc<0.5,"kingdom"]<-"unidentified"
  df.tmp[df.tmp$p_perc<0.5,"phylum"]<-"unidentified"
  df.tmp[df.tmp$c_perc<0.5,"class"]<-"unidentified"
  df.tmp[df.tmp$o_perc<0.5,"order"]<-"unidentified"
  df.tmp[df.tmp$f_perc<0.5,"family"]<-"unidentified"
  df.tmp[df.tmp$g_perc<0.5,"genus"]<-"unidentified"
  #df.tmp[df.tmp$s_perc<0.5,"genus_species_code"]<-"unidentified"
  
  #c. replace columns that include "unidentified with just unidentified"
  i<-0
  COLS <- c("kingdom","phylum","class","order","family","genus","species")
  for (i in 1:length(COLS)){
    #pull out the current column
    columnThing <- df.tmp[,COLS[i]]
    #find "unidentified"s and replace cell with simple "unidentified"
    columnThing[grepl("unidentified", columnThing)] <- "unidentified"
    #update the column in the original dataframe
    df.tmp[,COLS[i]] <- columnThing
  }
  
  #e. make a genus_species column
  df.tmp$genus_species <- paste(df.tmp$genus, df.tmp$species, sep="_")
  df.tmp[df.tmp$genus_species=="unidentified_sp","genus_species"]<-"Fungi_sp"
  
  #d. pull stuff together and make a 'taxonomy' column
  perc <- df.tmp$g_perc
  genus_species <- df.tmp$genus_species
  eucode <- "EU"
  shcode <- "SH"
  repinfo <- "reps"
  kingdom <- df.tmp$kingdom
  phylum <- df.tmp$phylum
  class <- df.tmp$class
  order <- df.tmp$order
  family <- df.tmp$family
  genus <- df.tmp$genus
  species <- df.tmp$species
  df.tmp[,"taxonomy"]<-paste(perc, 
                             genus_species, eucode, shcode, repinfo,
                             paste(paste("k__",kingdom, sep =""),
                                   paste("p__",phylum, sep =""),
                                   paste("c__",class, sep =""),
                                   paste("o__",order, sep =""),
                                   paste("f__",family, sep =""),
                                   paste("g__",genus, sep =""),
                                   paste("s__",species, " ...", sep =""), sep = ";"), 
                             sep = "|")
  
  #2. attach the 'taxonomy' column to the otu table
  ind.tax <- df.tmp[,c("OTUId","taxonomy")]
  fung.table <- merge(otu.table, ind.tax)
  
  
  return(fung.table)
}

FormatforFUNGuild <- function(select.table.name){
  
  # current df from RDP output
  fileName1 <- paste("./data/microbes/classifiedRDP/RDP_", select.table.name, ".txt", sep = "")
  rdp.table <- read.table(fileName1, sep = "\t", stringsAsFactors = FALSE)
  
  # current otu table
  fileName2 <- paste("./data/microbes/otuTables/otuTable_", select.table.name, ".tab", sep = "")
  otu.table <- read.table(fileName2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # format fxn based on the gene amplifed
  #ITS
  if( grepl("ITS",select.table.name) ){
    fung.table <- RDPtoFUNG_ITS(otu.table, rdp.table)
  }
  #LSU
  if( grepl("LSU",select.table.name) ){
    fung.table <- RDPtoFUNG_LSU(otu.table, rdp.table)
  }
  
  # write table
  fileName3 <- paste("./data_Rsynth/FUNGtables/FUNGtables_", select.table.name, ".tab", sep = "")
  write.table(fung.table, file=fileName3, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(print("done"))
  
}
