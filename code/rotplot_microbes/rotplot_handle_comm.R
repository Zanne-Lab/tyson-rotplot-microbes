# rotplot-specific: requires prep_otutabs_all.R, summarize_OTUmat_all.R

#------------------------------------------------------------------#
# Find an appropriate OTU table normalization

plotSampleEffort_rotplotNobags <- function(output.path, tab.t100){
  
  jpeg(paste0(output.path, "sampleEffortCurves.jpeg"), width = 600, height = 600)
  par(mfrow=c(2,2))
  tmp <- tab.t100
  rarecurve(tmp[['2012ITS']], step=100,
            xlab="Number of reads per sample", ylab="Cumulative number of OTUs", label=FALSE,
            main = "2012 ITS")
  rarecurve(tmp[['201216S']], step=100,
            xlab="Number of reads per sample", ylab="Cumulative number of OTUs", label=FALSE,
            main = "2012 16S")
  rarecurve(tmp[['2014ITS']], step=100,
            xlab="Number of reads per sample", ylab="Cumulative number of OTUs", label=FALSE,
            main = "2014 ITS")
  rarecurve(tmp[['201416S']], step=100,
            xlab="Number of reads per sample", ylab="Cumulative number of OTUs", label=FALSE,
            main = "2014 16S")
  dev.off()
  
}

compareOTUnorms_rotplotNobags <- function(output.path, tab.t100){
  
  ## Rarify OTU table
  require("phyloseq")
  tab.r500 <- MakeRare.OTUtab(curr.otuTab.list = tab.t100, sampleSize = 500)
  
  ## CSS normalize OTU table
  require("metagenomeSeq")
  tab.css100 <- MakeCSSNorm.OTUtab(curr.otuTab.list = tab.t100)
  
  ## Transform into presence/absence
  tab.t100.pa <- lapply(tab.t100, function(x){1 * (x > 0)})
  tab.r500.pa <- lapply(tab.r500, function(x){1 * (x > 0)})
  
  ## Trim OTUs that show up in less than 20% of the samples
  tab.t100.trim <- lapply(tab.t100, trimLowSampleOTUs)
  
  # Make quick ordinations
  p.t100 <- collect_plots(otu.list = tab.t100, presAbs = F, title = "t100")
  p.t100.trim <- collect_plots(otu.list = tab.t100.trim, presAbs = F, title = "t100.trim")
  p.r500 <- collect_plots(otu.list = tab.r500, presAbs = F, title = "r500")
  p.t100.pa <- collect_plots(otu.list = tab.t100.pa, presAbs = T, title = "t100.pa")
  p.r500.pa <- collect_plots(otu.list = tab.r500.pa, presAbs = T, title = "r500.pa")
  p.css100 <- collect_plots(otu.list = tab.css100, presAbs = F, title = "css100")
  
  jpeg(paste0(output.path, "otuNorm_comparisons.jpeg"), width = 1000, height = 1000)
  grid.arrange(p.t100[['2012ITS']],
               p.t100.trim[['2012ITS']],
               p.r500[['2012ITS']],
               p.t100.pa[['2012ITS']],
               p.r500.pa[['2012ITS']],
               p.css100[['2012ITS']],
               
               p.t100[['2014ITS']],
               p.t100.trim[['2014ITS']],
               p.r500[['2014ITS']],
               p.t100.pa[['2014ITS']],
               p.r500.pa[['2014ITS']],
               p.css100[['2014ITS']],
               
               p.t100[['201216S']],
               p.t100.trim[['201216S']],
               p.r500[['201216S']],
               p.t100.pa[['201216S']],
               p.r500.pa[['201216S']],
               p.css100[['201216S']],
               
               p.t100[['201416S']],
               p.t100.trim[['201416S']],
               p.r500[['201416S']],
               p.t100.pa[['201416S']],
               p.r500.pa[['201416S']],
               p.css100[['201416S']],
               
               nrow = 4)
  dev.off()
  
  
}

#------------------------------------------------------------------#
# Summarize OTU tables

summarizeMicrobes_rotplot <- function(data.overlap, data.deploy, tab.taxa, otuPrep){
  
  # pair each OTU table with its taxa table
  datalist.overlap <- pair_otu_and_taxtab(datalist = data.overlap, tab.taxa = tab.taxa)
  datalist.deploy <- pair_otu_and_taxtab(datalist = data.deploy, tab.taxa = tab.taxa)
  
  # summarize each element in the list for OTU richness and taxon coverage
  summ.overlap <- lapply(datalist.overlap, function(x){
    summaryStats <- calc_otu_richnessCoverage(curr.otu = x$curr.otu, 
                                              curr.tab.taxa = x$curr.tab.taxa,
                                              fung = x$fung)
    return(summaryStats)
  })
  summ.overlap.df <- list_to_df(summ.overlap)
  
  summ.deploy <- lapply(datalist.deploy, function(x){
    summaryStats <- calc_otu_richnessCoverage(curr.otu = x$curr.otu, 
                                              curr.tab.taxa = x$curr.tab.taxa,
                                              fung = x$fung)
    return(summaryStats)
  })
  summ.deploy.df <- list_to_df(summ.deploy)
  
  otusumm <- rbind(summ.overlap.df, summ.deploy.df)
  tmp <- data.frame(otusumm, row.names = NULL) 
  tmp %>%
    spread(key = label, value = value) %>%
    select(source, meanRichness, seRichness, nRichness, totalOTUs, num.genus.unclass, num.guild.unclass) -> otusumm.df
  otusumm.df$otunorm <- otuPrep
  
  return(otusumm.df)
}

make_otusumm <- function(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa){
  
  # OTU richness, taxonomic and functional information coverage for each data subset
  otusumm.list <- list()
  for(i in 1:length(data.deploy.otuPreplist)){
    otusumm <- summarizeMicrobes_rotplot(data.overlap = data.overlap.otuPreplist[[i]], 
                                         data.deploy = data.deploy.otuPreplist[[i]], 
                                         tab.taxa = tab.taxa, 
                                         otuPrep = names(data.deploy.otuPreplist)[i])
    otusumm %>%
      separate(source, into = c("yearGene","subset")) %>%
      separate(yearGene, into = c("year","gene"), sep = 4) -> otusumm.list[[i]]
  }
  otusumm <- list_to_df(otusumm.list)
  
  return(otusumm)
  
}

compareNorms_totalOTUs <- function(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa){
  
  otusumm <- make_otusumm(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa)
  otusumm.totalOTUs <- unique(otusumm[,c("otunorm","year","gene","totalOTUs","num.genus.unclass","num.guild.unclass")])
  
  # for fungi, how many genus-ID OTUs also have a guild?
  otusumm.totalOTUs %>%
    filter(gene == "ITS") %>%
    mutate(num.genus.class = totalOTUs - num.genus.unclass) %>%
    mutate(num.guild.class = num.genus.class - num.guild.unclass) %>% # num.guild.unclass = if it has a genus-level id, how many do not have a guild id?
    mutate(per.guild.class = round((num.guild.class/num.genus.class) * 100, digits = 1))
  # if OTU has a genus-ID, then >95% also have a guild id
  
  # plot the total number of OTUs detected per gene, year, and otu normalization
  otusumm.totalOTUs$otunorm <- factor(otusumm.totalOTUs$otunorm, levels = c("t100pa","r100pa","r500pa","r1000pa"))
  otusumm.totalOTUs %>%
    mutate(num.genus.class = totalOTUs - num.genus.unclass) %>%
    select(otunorm, year, gene, num.genus.unclass, num.genus.class) %>%
    gather(key = "measure", value = "value", -c(otunorm, year, gene)) -> tmp
  p <- ggplot(data = tmp, 
              aes(x = otunorm, y = value, 
                  alpha = measure)) +
    geom_bar(stat = "identity") +
    facet_wrap(gene~year, scales = "free") +
    ylab("Total number of OTUs") + xlab("Normalization type") +
    scale_alpha_manual(name = "Genus-level identity", values = c(1,.5), labels = c("yes","no")) + 
    theme_bw()
  
  return(p)
  
}

compareNorms_numSamps <- function(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa){
  
  otusumm <- make_otusumm(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa)
  
  # make a table of the number of samples per gene, year, subset, and otu normalization
  otusumm.n <- otusumm[,c("year","gene","subset","otunorm","nRichness")]
  otusumm.n %>%
    spread(key = otunorm, value = nRichness) %>%
    select(year, gene, subset, t100pa, r100pa, r500pa, r1000pa) %>%
    arrange(gene, subset) -> otusumm.n
  
  return(otusumm.n)
}

compareNorms_sampRich <- function(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa){
  
  otusumm <- make_otusumm(data.deploy.otuPreplist, data.overlap.otuPreplist, tab.taxa)
  
  # plot mean OTU richness by gene, year, subset, and otu normalization
  otusumm %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> otusumm
  params <- plottingParams_rotplotNobags()
  otusumm$yearSubset <- factor(otusumm$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  otusumm$otunorm <- factor(otusumm$otunorm, levels = c("t100pa","r100pa","r500pa","r1000pa"))
  otusumm %>%
    filter(subset != "o") -> tmp
  p <- ggplot(data = tmp, aes(x = otunorm, y = meanRichness, 
                              color = yearSubset)) +
    geom_point() +
    facet_wrap(~gene, scales = "free") +
    geom_errorbar(aes(ymin = meanRichness - seRichness, 
                      ymax = meanRichness + seRichness), width =.2) +
    ylab("OTU richness per sample") + xlab("Normalization type") +
    theme_bw()
  
  return(p)
  
}

richnessDistrib_rotplot <- function(data.deploy){
  
  # calculate sample richness
  rich.list <- list()
  for(i in 1:length(data.deploy)){
    x <- data.deploy[[i]]
    richness <- apply(x$otu, MARGIN = 1,function(x) sum(x>0))
    rich.list[[i]] <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
  }
  names(rich.list) <- names(data.deploy)
  rich.df <- list_to_df(rich.list)
  rich.df %>%
    separate(source, into = c("yearGene","subset")) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) -> rich.df
  
  rich.df %>%
    group_by(year, gene) %>%
    summarize(mean = mean(rich),
              n = sum(!is.na(rich)),
              se = sd(rich)/sqrt(n))
  
  # plot
  rich.df %>%
    filter(gene =="ITS") -> tmp
  p.fung <- ggplot(tmp, aes(x = rich)) +
    geom_histogram() +
    facet_grid(~year) +
    ylab("Frequency") + xlab("Sample OTU richness") +
    theme_bw() + ggtitle("Fungi")
  
  rich.df %>%
    filter(gene =="16S") -> tmp
  p.bact <- ggplot(tmp, aes(x = rich)) +
    geom_histogram() +
    facet_grid(~year) +
    ylab("Frequency") + xlab("Sample OTU richness") +
    theme_bw() + ggtitle("Bacteria")
  
  p.list <- list(p.fung = p.fung, p.bact = p.bact)
  
  return(p.list)
}

freqPhyla <- function(tab.otu, tab.taxa){
  
  # commdata and tab.taxa need to be lists in the same order
  sum(names(tab.otu) != names(tab.taxa)) # this needs to be 0
  
  # calculate the number of OTU occurances
  tab.occur <- lapply(tab.otu, function(x){
    occur <- apply(x, 2, function(k){sum(k > 0)})
    occur.df <- data.frame(OTUid = names(occur), occur)
  })
  df.occur <- list_to_df(tab.occur)
  df.occur %>%
    rename('table' = 'source') -> df.occur
  
  # index OTUids for all tables
  tab.taxa.tmp <- lapply(tab.taxa, function(x){
    x %>%
      select(OTUid, genusSpecies, kingdom, phylum, class, order, family, genus) -> df
  })
  df.taxa <- list_to_df(tab.taxa.tmp)
  df.taxa %>%
    rename('table' = 'source') -> df.taxa
  
  # join OTU occurances with taxon info
  df.occur %>%
    left_join(df.taxa) -> df
  
  # identify highly cosmopolitan OTUs
  # determine the 95% quantile for OTU occurances in each table
  df %>%
    group_by(table) %>%
    summarize(q95 = quantile(occur, .95)) -> q95.indx
  q95.list <- list()
  for(i in 1:dim(q95.indx)[1]){
    df %>%
      filter(table == q95.indx[i,"table"]$table) %>%
      filter(occur > q95.indx[i, "q95"]$q95) %>%
      mutate(cat = "Highly cosmopolitan") -> q95.list[[i]]
  }
  names(q95.list) <- q95.indx$table
  q95.df <- list_to_df(q95.list)
  
  # identify rare OTUs (ie OTUs that only show up in 1 sample)
  rare.list <- list()
  for(i in 1:dim(q95.indx)[1]){
    df %>%
      filter(table == q95.indx[i,"table"]$table) %>%
      filter(occur == 1) %>%
      mutate(cat = "Rare") -> rare.list[[i]]
  }
  names(rare.list) <- q95.indx$table
  rare.df <- list_to_df(rare.list)
  
  result.df <- rbind(q95.df, rare.df)
  result.df %>%
    select(-source) -> result.df
  
  # simplify the phyla categories
  
  # how many OTUs from each phylum? (treat each table separately)
  df %>%
    group_by(table, phylum) %>%
    summarize(n = length(unique(OTUid))) -> tmp
  
  #bacteria
  tmp %>%
    filter(table == "201216S") %>%
    arrange(desc(n)) -> tmp.o
  tmp %>%
    filter(table == "201416S") %>%
    arrange(desc(n))
  keepPhyl <- tmp.o$phylum[1:8]
  keepPhyl1 <- c(as.character(keepPhyl), 'other')
  #fungi
  tmp %>%
    filter(table == "2012ITS") %>%
    arrange(desc(n))
  tmp %>%
    filter(table == "2014ITS") %>%
    arrange(desc(n))
  keepPhyl.fung <- c("Ascomycota","Basidiomycota","unidentified")
  keepPhyl1.fung <- c(as.character(keepPhyl.fung), 'other')
  
  #update df
  result.df %>%
    separate(table, into = c("yrseq","gene"), sep = 4, remove = F) -> result.df
  result.df %>%
    filter(gene == "16S") %>%
    mutate(phylum.select = ifelse(gene == "16S" & phylum %in% as.character(keepPhyl), 
                                  as.character(phylum), 'other')) -> result.df.bact
  result.df %>%
    filter(gene == "ITS") %>%
    mutate(phylum.select = ifelse(gene == "ITS" & phylum %in% as.character(keepPhyl.fung), 
                                  as.character(phylum), 'other')) -> result.df.fung
  result <- rbind(result.df.bact, result.df.fung)
  
  #summarize by source and phylum
  result %>%
    group_by(gene, yrseq, phylum.select, cat) %>%
    summarize(notus = length(unique(OTUid))) -> summ.df
  params <- plottingParams_rotplotNobags()
  
  #bacteria
  summ.df %>%
    filter(gene == "16S") -> tmp
  tmp$phylum.select <- factor(tmp$phylum.select, levels = keepPhyl1)
  p.bact <- ggplot(data = tmp, aes(x = cat, y = notus, fill = phylum.select)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab("Number of OTUs") +
    theme_bw() +
    scale_fill_discrete(name = "Phylum") +
    facet_grid(~yrseq)
  p.bact
  
  #fungi
  summ.df %>%
    filter(gene == "ITS") -> tmp
  tmp$phylum.select <- factor(tmp$phylum.select, levels = keepPhyl1.fung)
  p.fung <- ggplot(data = tmp, aes(x = cat, y = notus, fill = phylum.select)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab("Number of OTUs") +
    theme_bw() +
    scale_fill_discrete(name = "Phylum") +
    facet_grid(~yrseq)
  p.fung
  
  return.result <- list(p.fung = p.fung, p.bact = p.bact, df = result)
  
  return(return.result)
  }

# calc_sampleRichness_summByrotplot <- function(curr.otu, curr.meta){
#   
#   richness <- apply(curr.otu, MARGIN = 1,function(x) sum(x>0))
#   richness.df <- data.frame(tidySeq_sampleName = names(richness), richness)
#   richness.df %>%
#     left_join(curr.meta) %>%
#     select(tidySeq_sampleName, topo, logLoc, logtrt, richness) -> richness.df
#   
#   return(richness.df)
# }

plotRichness_rotplot_deploy <- function(data.deploy){
  
  # calculate sample richness and put it into a df with study metadata
  rich.deploy <- lapply(data.deploy, function(x){
    sampleRich <- calc_sampleRichness_summByrotplot(curr.otu = x$otu, curr.meta = x$meta)
    return(sampleRich)
  })
  rich.deploy.df <- list_to_df(rich.deploy)
  # identify source, add binomials
  logTrt <- load_logTrt()
  logTrt.indx <- logTrt[,c("logtrt","species","species4","binomial")]
  rich.deploy.df %>%
    separate(source, into = c("table","subset"), remove = F) %>%
    separate(table, into = c("year","gene"), sep = 4, remove = F) %>%
    left_join(logTrt.indx) -> rich.deploy.df
  
  #plot
  params <- plottingParams_rotplotNobags()

  # aggregate by topo, logLoc, and binomial
  rich.deploy.df %>%
    group_by(topo, logLoc, binomial, gene, year, subset) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(topoLogLoc = paste(topo, logLoc, sep = "_")) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) %>%
    filter(logLoc != "mush") -> rich.sum.df
  rich.sum.df$gene <- factor(rich.sum.df$gene, levels = c("ITS","16S"), labels = c("Fungi","Bacteria"))
  rich.sum.df$yearSubset <- factor(rich.sum.df$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  rich.sum.df$binomial <- factor(rich.sum.df$binomial, levels = params$binom.levels)
  p <- ggplot(data = rich.sum.df, aes(y = mean, x = binomial, 
                                            shape = logLoc,
                                            color = topo)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    #geom_line(aes(group = topoLogLoc), alpha = .5) +
    facet_wrap(gene~yearSubset, scales = "free_x", ncol = 4) + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    scale_shape_manual(name = "Log side",
                       values = params$logLoc$shapes) +
    scale_color_manual(name = "Plot location",
                       values = params$topo$colors) +
    theme_bw()
  p
  
  # aggregate by topo and logLoc
  rich.deploy.df %>%
    group_by(topo, logLoc, gene, year, subset) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(topoLogLoc = paste(topo, logLoc, sep = "_")) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) %>%
    filter(logLoc != "mush")-> rich.sum.df.pos
  
  rich.sum.df.pos$gene <- factor(rich.sum.df.pos$gene, levels = c("ITS","16S"), labels = c("Fungi","Bacteria"))
  rich.sum.df.pos$yearSubset <- factor(rich.sum.df.pos$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  rich.sum.df.pos$topoLogLoc <- factor(rich.sum.df.pos$topoLogLoc, 
                                       levels = c("H_t","H_b","L_t","L_b"), 
                                       labels = c("Ridge plot, log top",
                                                  "Ridge plot, log bottom",
                                                  "Valley plot, log top",
                                                  "Valley plot, log bottom"))
  p.pos <- ggplot(data = rich.sum.df.pos, aes(y = mean, x = topoLogLoc, 
                                      shape = logLoc,
                                      color = topo)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    facet_wrap(gene~yearSubset, scales = "free_x", ncol = 4) + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    scale_shape_manual(name = "Log side",
                       values = params$logLoc$shapes) +
    scale_color_manual(name = "Plot location",
                       values = params$topo$colors) +
    theme_bw() + guides(color = F, shape = F)
  p.pos
  
  # aggregate by binomial
  rich.deploy.df %>%
    group_by(binomial, gene, year, subset) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> rich.sum.df.binom
  rich.sum.df.binom$gene <- factor(rich.sum.df.binom$gene, levels = c("ITS","16S"), labels = c("Fungi","Bacteria"))
  rich.sum.df.binom$yearSubset <- factor(rich.sum.df.binom$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  rich.sum.df.binom$binomial <- factor(rich.sum.df.binom$binomial, levels = params$binom.levels)
  
  p.binom <- ggplot(data = rich.sum.df.binom, aes(y = mean, x = binomial)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    facet_wrap(gene~yearSubset, ncol = 4, scales = "free_x") + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    theme_bw()
  p.binom
  
  p.list <- list(p = p, p.pos = p.pos, p.binom = p.binom)
  return(p.list)
  
}

plotRichness_rotplot_deploy_otuPrep <- function(data.deploy.otuPreplist){
  
  rich.deploy.df.list <- list()
  logTrt <- load_logTrt()
  logTrt.indx <- logTrt[,c("logtrt","species","species4","binomial")]
  for(i in 1:length(data.deploy.otuPreplist)){
    
    data.deploy <- data.deploy.otuPreplist[[i]]
    # calculate sample richness and put it into a df with study metadata
    rich.deploy <- lapply(data.deploy, function(x){
      sampleRich <- calc_sampleRichness_summByrotplot(curr.otu = x$otu, curr.meta = x$meta)
      return(sampleRich)
    })
    rich.deploy.df <- list_to_df(rich.deploy)
    # identify source, add binomials
    rich.deploy.df %>%
      separate(source, into = c("table","subset")) %>%
      separate(table, into = c("year","gene"), sep = 4, remove = F) %>%
      left_join(logTrt.indx) -> rich.deploy.df
    rich.deploy.df.list[[i]] <- rich.deploy.df
    
  }
  names(rich.deploy.df.list) <- names(data.deploy.otuPreplist)
  rich.deploy.df <- list_to_df(rich.deploy.df.list)
  rich.deploy.df %>%
    rename('otuprep'='source') -> rich.deploy.df
  
  #plot
  params <- plottingParams_rotplotNobags()
  
  # aggregate by topo and logLoc
  rich.deploy.df %>%
    group_by(topo, logLoc, gene, year, subset, otuprep) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(topoLogLoc = paste(topo, logLoc, sep = "_")) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) %>%
    filter(logLoc != "mush")-> rich.sum.df.pos
  rich.sum.df.pos$gene <- factor(rich.sum.df.pos$gene, levels = c("ITS","16S"), labels = c("Fungi","Bacteria"))
  rich.sum.df.pos$yearSubset <- factor(rich.sum.df.pos$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  rich.sum.df.pos$topoLogLoc <- factor(rich.sum.df.pos$topoLogLoc, 
                                       levels = c("H_t","H_b","L_t","L_b"), 
                                       labels = c("Ridge plot, log top",
                                                  "Ridge plot, log bottom",
                                                  "Valley plot, log top",
                                                  "Valley plot, log bottom"))
  p.pos <- ggplot(data = rich.sum.df.pos, aes(y = mean, x = topoLogLoc, color = yearSubset, group = yearSubset)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    facet_wrap(gene~otuprep, scales = "free_y", ncol = 4)  +
    xlab("") + ylab("OTU richness per sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p.pos
  
  # aggregate by binomial
  rich.deploy.df %>%
    group_by(binomial, gene, year, subset, otuprep) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> rich.sum.df.binom
  rich.sum.df.binom$gene <- factor(rich.sum.df.binom$gene, levels = c("ITS","16S"), labels = c("Fungi","Bacteria"))
  rich.sum.df.binom$yearSubset <- factor(rich.sum.df.binom$yearSubset, levels = params$yearSubset$levels, labels = params$yearSubset$labels)
  rich.sum.df.binom$binomial <- factor(rich.sum.df.binom$binomial, levels = params$binom.levels)
  
  p.binom <- ggplot(data = rich.sum.df.binom, aes(y = mean, x = binomial, color = yearSubset, group = yearSubset)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    facet_wrap(gene~otuprep, ncol = 4, scales = "free_x") + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    theme_bw()
  p.binom
  
  result <- list(p.pos = p.pos, p.binom = p.binom)
  return(result)
  
}

plotRichness_rotplot_overlap <- function(data.overlap){
  
  # calculate sample richness and put it into a df with study metadata
  rich.overlap <- lapply(data.overlap, function(x){
    sampleRich <- calc_sampleRichness_summByrotplot(curr.otu = x$otu, curr.meta = x$meta)
    return(sampleRich)
  })
  rich.overlap.df <- list_to_df(rich.overlap)
  
  # identify source, add binomials
  logTrt <- load_logTrt()
  logTrt.indx <- logTrt[,c("logtrt","species","species4","binomial","yrdeploy")]
  rich.overlap.df %>%
    separate(source, into = c("table","subset"), remove = F) %>%
    separate(table, into = c("year","gene"), sep = 4, remove = F) %>%
    left_join(logTrt.indx) %>%
    mutate(deployInt = paste(yrdeploy, year, sep ="_")) -> rich.overlap.df
  
  #plot
  params <- plottingParams_rotplotNobags()
  #params$binom.levels
  
  # aggregate by topo, logLoc, and binomial
  rich.overlap.df %>%
    group_by(topo, logLoc, binomial, gene, deployInt) %>%
    summarize(mean = mean(richness),
              n = sum(!is.na(richness)),
              se = sd(richness)/sqrt(n)) %>%
    mutate(topoLogLoc = paste(topo, logLoc, sep = "_")) %>%
    filter(logLoc != "mush") %>%
    separate(deployInt, into = c("yrdeploy","year"), remove = F)-> rich.sum.df
  
  rich.sum.df %>%
    filter(gene == "ITS") -> tmp
  p.fung <- ggplot(data = tmp, aes(y = mean, x = yrdeploy, 
                                      shape = logLoc,
                                      color = topo)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    geom_line(aes(group = topoLogLoc), alpha = .5) +
    facet_grid(binomial~year) + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    scale_shape_manual(name = "Log side",
                       values = params$logLoc$shapes) +
    scale_color_manual(name = "Plot location",
                       values = params$topo$colors) +
    theme_bw()
  p.fung
  
  rich.sum.df %>%
    filter(gene == "16S") -> tmp
  p.bact <- ggplot(data = tmp, aes(y = mean, x = yrdeploy, 
                                   shape = logLoc,
                                   color = topo)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2, alpha =.5) +
    geom_line(aes(group = topoLogLoc), alpha = .5) +
    facet_grid(binomial~year) + coord_flip() +
    xlab("") + ylab("OTU richness per sample") +
    scale_shape_manual(name = "Log side",
                       values = params$logLoc$shapes) +
    scale_color_manual(name = "Plot location",
                       values = params$topo$colors) +
    theme_bw()
  p.bact
  
  p.list <- list(p.fung = p.fung, p.bact = p.bact)
  return(p.list)
}


