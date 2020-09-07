

#------------------------------------------------------------------#
# Find an appropriate OTU table normalization

quick_ord <- function(otuTab, presAbs, title){
  
  if(presAbs == TRUE){
    mod.obj <- capscale(otuTab ~ 1, data = NULL, distance = "jaccard", binary = TRUE)
  }else{
    otuTab <- decostand(otuTab, 'hellinger') # transform otu tab
    mod.obj <- capscale(otuTab ~ 1, data = NULL, distance = "bray")
  }
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  p <- ggplot(data = scrs, aes(x = MDS1, y = MDS2)) +
    geom_point(alpha = .5) +
    coord_fixed() +
    xlab("PCoA 1") + ylab("PCoA 2") + theme_bw() +
    ggtitle(title)
  
  return(p)
}

collect_plots <- function(otu.list, presAbs, title){
  
  p.list <- list()
  for(i in 1:length(otu.list)){
    mat.otu <- otu.list[[i]]
    p.list[[i]] <- quick_ord(otuTab = mat.otu, presAbs = presAbs, title = title)
  }
  names(p.list) <- names(otu.list)
  return(p.list)
}

#------------------------------------------------------------------#
# Calculate OTU table summaries

calc_otu_summStats<-function(comm.otu){
  
  require(vegan)
  
  #total number of OTUs
  totalOTUs<-dim(comm.otu)[2]
  
  #mean sample richness
  richness<-apply(comm.otu, MARGIN = 1,function(x) sum(x>0))
  meanRichness<-mean(richness)
  seRichness<-sd(richness)/sqrt(length(richness))
  
  # #mean sample evenness
  # H <- diversity(comm.otu) 
  # S <- specnumber(comm.otu) ## rowSums(BCI > 0) does the same...
  # J <- H/log(S) #Pielou's evenness (J)
  # meanJ<-mean(J)
  # seJ<-sd(J)/sqrt(length(J))
  
  #mean number of reads per sample
  meanReads<-mean(rowSums(comm.otu))
  seReads<-sd(rowSums(comm.otu))/sqrt(length(rowSums(comm.otu)))
  
  summaryStats<-data.frame(label=c("totalOTUs","meanRichness","seRichness","meanReads","seReads"),
                           value=c(totalOTUs, meanRichness, seRichness, meanReads, seReads))
  summaryStats$value <- round(summaryStats$value, digits=4)
  
  return(summaryStats)
}

calc_otu_richnessCoverage <-function(curr.otu, curr.tab.taxa, fung){
  
  #total number of OTUs
  totalOTUs <- dim(curr.otu)[2]
  
  #sample richness
  richness <- apply(curr.otu, MARGIN = 1,function(x) sum(x>0))
  nRichness <- sum(!is.na(richness))
  meanRichness<-mean(richness)
  seRichness<-sd(richness)/sqrt(nRichness)
  
  # taxonomic and functional coverage
  curr.tab.taxa %>%
    filter(OTUid %in% colnames(curr.otu)) -> tax.indx
  tax.indx %>%
    filter(genus == "unidentified") -> genus.unclass
  num.genus.unclass <- dim(genus.unclass)[1]
  
  if(fung == TRUE){
    tax.indx %>%
      filter(genus != "unidentified") %>%
      filter(is.na(Trophic.Mode)) -> guild.unclass
    num.guild.unclass <- dim(guild.unclass)[1]
    summaryStats<-data.frame(label=c("totalOTUs","meanRichness","seRichness","nRichness","num.genus.unclass","num.guild.unclass"),
                             value=c(totalOTUs, meanRichness, seRichness, nRichness, num.genus.unclass, num.guild.unclass))
  }else{
    summaryStats<-data.frame(label=c("totalOTUs","meanRichness","seRichness","nRichness","num.genus.unclass","num.guild.unclass"),
                             value=c(totalOTUs, meanRichness, seRichness, nRichness, num.genus.unclass, NA))
  }

  summaryStats$value <- round(summaryStats$value, digits=4)
  
  return(summaryStats)
}


#------------------------------------------------------------------#
# Formatting

pair_otu_and_taxtab <- function(datalist, tab.taxa){
  
  otutaxa.list <- list()
  for(i in 1:length(datalist)){
    
    curr.otu <- datalist[[i]]$otu
    curr.otu.name <- names(datalist)[i]
    
    # pull out the correct taxon table
    # for rotplot
    if(grepl("ITS", curr.otu.name) == TRUE & grepl("2012", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "2012ITS")]]
      fung <- TRUE
    }
    if(grepl("ITS", curr.otu.name) == TRUE & grepl("2014", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "2014ITS")]]
      fung <- TRUE
    }
    if(grepl("16S", curr.otu.name) == TRUE & grepl("2012", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "201216S")]]
      fung <- FALSE
    }
    if(grepl("16S", curr.otu.name) == TRUE & grepl("2014", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "201416S")]]
      fung <- FALSE
    }
    # for cwd
    if(grepl("ITS", curr.otu.name) == TRUE & grepl("2017", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "2017ITS_CWD")]]
      fung <- TRUE
    }
    if(grepl("16SB", curr.otu.name) == TRUE & grepl("2017", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "201716SB_CWD")]]
      fung <- FALSE
    }
    if(grepl("16SA", curr.otu.name) == TRUE & grepl("2017", curr.otu.name) == TRUE){
      curr.tab.taxa <- tab.taxa[[which(names(tab.taxa) == "201716SA_CWD")]]
      fung <- FALSE
    }
    
    otutaxa.list[[i]] <- list(curr.otu = curr.otu, curr.tab.taxa = curr.tab.taxa, fung = fung)
    
  }
  names(otutaxa.list) <- names(datalist)
  
  return(otutaxa.list)
  
}

sameSamples <- function(commY, commX){
    
    #make sure commX and commY have the same samples
    commy.otu <- commY$otu
    commx.otu <- commX$otu
    all.rownames <- c(row.names(commy.otu), row.names(commx.otu))
    common.rownames <- all.rownames[duplicated(all.rownames)]
    commy.otu.trim <- commy.otu[row.names(commy.otu) %in% common.rownames,]
    commx.otu.trim <- commx.otu[row.names(commx.otu) %in% common.rownames,]
    dim(commy.otu.trim); dim(commx.otu.trim)
    
    #check for empty samples
    commy.otu.trim1 <- remove_zeroSamprows(commy.otu.trim)
    commx.otu.trim1 <- remove_zeroSamprows(commx.otu.trim)
    dim(commy.otu.trim1); dim(commx.otu.trim1)
    
    #check that they still have the same samples
    if(dim(commy.otu.trim1)[1] != dim(commx.otu.trim1)[1]){
        
        all.rownames <- c(row.names(commy.otu.trim1), row.names(commx.otu.trim1))
        common.rownames <- all.rownames[duplicated(all.rownames)]
        commy.otu.trim1 <- commy.otu.trim1[row.names(commy.otu.trim1) %in% common.rownames,]
        commx.otu.trim1 <- commx.otu.trim1[row.names(commx.otu.trim1) %in% common.rownames,]
        dim(commy.otu.trim1); dim(commx.otu.trim1)
        
    }
    result <- list(commy.otu = commy.otu.trim1, commx.otu = commx.otu.trim1)
    
    return(result)
    
}

otu_to_taxon.tab <- function(mat.otu.pres, taxa.df){
  
  # turn the OTU table into a taxon table
  otu.t <- t(mat.otu.pres) # tranform
  
  # remove OTU taxa rows where then taxon id is na or "possible" (keep "Probable" and "Highly Probable")
  taxa.df %>%
    filter(!is.na(Taxon)) %>%
    filter(Confidence.Ranking %in% c("Probable","Highly Probable")) -> new.taxa.df
  
  # for each unique Taxon, select relevant OTUs
  TAX <- unique(new.taxa.df$Taxon)[!is.na(unique(new.taxa.df$Taxon))]
  tax.list <- list()
  for(i in 1:length(TAX)){
    new.taxa.df %>%
      filter(Taxon == TAX[i]) -> tmp
    
    if(sum(row.names(otu.t) %in% tmp$OTUid) != 0){
      sub.otu.t <- otu.t[row.names(otu.t) %in% tmp$OTUid,]
      
      if(sum(row.names(otu.t) %in% tmp$OTUid) != 1 & dim(tmp)[1] != 1){
        vec <- colSums(sub.otu.t) # sum the number of OTUs per Taxon per sample
      }else{
        vec <- sub.otu.t
      }
      df <- data.frame(tidySeq_sampleName = names(vec), 
                       nOTUs = vec, 
                       Taxon = rep(as.character(TAX[i]), length(vec)), 
                       row.names = NULL)
      tax.list[[i]] <- df
    }
  }
  tax.df <- list_to_df(tax.list)
  tax.df %>%
    spread(key = Taxon, value = nOTUs) -> tax.df.w
  row.names(tax.df.w) <- tax.df.w$tidySeq_sampleName
  tax.df.w <- tax.df.w[,-1]
  tax.df.w.pa <- 1 * (tax.df.w < 1) # transform into Taxon presence/absence
  
  return(tax.df.w.pa)
  
}


