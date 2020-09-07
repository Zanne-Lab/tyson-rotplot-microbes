plottingParams_rotplotNobags<-function(){
  
  require(colorspace)
  
  # ----------------------------------------- #
  #logLoc = "Within stem"
  logLocLevels <- c("t","mush","b")
  logLocShapes <- c("t"=17, "mush"=16, "b"=6)
  logLocLabels <- c("top","mush","bottom")
  logLoc <- list(levels=logLocLevels, labels=logLocLabels, shapes=logLocShapes)
  
  # ----------------------------------------- #
  #topo = "Within watershed")
  topoLevels <- c("H","L")
  topoColors <- c("#808080", "#000000") # darkgray, black
  topoLabels <- c("ridge","valley")
  names(topoColors) <- topoLevels
  topo <- list(levels=topoLevels, labels=topoLabels, colors=topoColors)
  
  # ----------------------------------------- #
  #species4 = "Wood species"
  species4Levels <- c("ACRU","AEGL","AMAR","ASTR","CATO","CEDA","CEOC","CHRY","COFL","DIVI","FRAM","GLTR","HDMP","HICK","JUNI","JUVI","LOMA","PIEC","PINE","PIST","PLOC","PRSE/PRVI","QUAL","QUVE","ROAK","ROSE","ULRU","VIVU","ZEBRA")
  species4Colors <- rainbow_hcl(29)
  names(species4Colors) <- species4Levels
  species4 <- list(levels=species4Levels, labels=species4Levels, colors=species4Colors)
  
  # ----------------------------------------- #
  # overlap data subsets
  overlapSpLevels<-c("CEOC","JUVI","QUVE")
  overlapSpColors<-rainbow_hcl(3)
  names(overlapSpColors)<-overlapSpLevels
  overlapSp<-list(levels=overlapSpLevels, labels=overlapSpLevels, colors=overlapSpColors)
  # add colors for species4topo levels
  overlapSpColors6 <- c(rainbow_hcl(3), "#a05564","#385130","#355470")
  names(overlapSpColors6) <- c("CEOC_L", "JUVI_L","QUVE_L", "CEOC_H", "JUVI_H","QUVE_H")
    
  #identify yrdeploy2009 species
  species09Levels<-c("ACRU","AEGL","AMAR","ASTR" ,"CATO","CEOC","COFL","DIVI","GLTR","JUVI","PIST","PLOC" ,"PRSE/PRVI","QUVE","ULRU","VIVU")
  species09Colors<-rainbow_hcl(16)
  names(species09Colors)<-species09Levels
  species09<-list(levels=species09Levels, labels=species09Levels, colors=species09Colors)

  #identify yrdeploy2011 species
  species11Levels<-c("CEOC","FRAM", "JUNI", "JUVI", "LOMA", "PIEC", "QUAL", "QUVE")
  species11Colors<-rainbow_hcl(8)
  names(species11Colors)<-species11Levels
  species11<-list(levels=species11Levels, labels=species11Levels, colors=species11Colors)
    
  #species order for sample design plot
  uniq09Sp <- species09Levels[!species09Levels %in% species11Levels]
  uniq11Sp <- species11Levels[!species11Levels %in% species09Levels]
  orderSp <- c(uniq09Sp, overlapSp$levels, uniq11Sp)
  orderSp.df <- data.frame(species4 = orderSp, order = seq.int(1:length(orderSp)))
  logtrt.df <- load_logTrt()
  orderSp.df %>%
    left_join(logtrt.df) %>%
    select(species4, binomial, angio.gymno, order) %>%
    mutate(binomial = sub("_"," ", binomial)) %>%
    arrange(order) -> orderSp.df
  orderSp.df <- unique(orderSp.df)

    # # order the Binomial names by the phylo tree
    # zanneTree <- load_zanne_tree()
    # require(ggtree)
    # ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
    # ggt$data %>%
    # filter(isTip == T) %>%
    # select(y, label) %>%
    # arrange(y) -> binom.levels
    
    # order the deployment subsets
    yearSubsetLevels2 <- c("2011_2012","2009_2012","2011_2014","2009_2014")
    yearSubsetLevels <- c("2012_d11","2012_d09","2014_d11","2014_d09")
    yearSubsetLabels <- c("2011-12 (1 yr)","2009-12 (3 yrs)","2011-14 (3 yrs)","2009-14 (5 yrs)")
    yearSubsetLevels2_new <- c("2011_2012","2011_2014","2009_2012","2009_2014")
    yearSubsetLevels_new <- c("2012_d11","2014_d11","2012_d09","2014_d09")
    yearSubsetLabels_new <- c("2011-12 (1 yr)","2011-14 (3 yrs)","2009-12 (3 yrs)","2009-14 (5 yrs)")
    yearSubset <- list(levels = yearSubsetLevels, labels = yearSubsetLabels, levels2 = yearSubsetLevels2,
                       levels.new = yearSubsetLevels_new, labels.new = yearSubsetLabels_new, levels2.new = yearSubsetLevels2_new)
    
    # wood traits
    trait.levels <- c("density","crudeProt","lignin","cellulose","hemicellulose","P","Ca","Mn","parenchymaFrac","conduitDiam","conduitLeng")
    trait.labels <- c("Density","Crude protein", "Lignin","Cellulose","Hemicellulose","P","Ca","Mn","Parenchyma","Conduit diameter", "Conduit length")
    traits <- list(levels = trait.levels, labels = trait.labels)
    
    # otuprep types
    #otuprep.levels2 <- c("t100pa","r100pa","r500pa","r1000pa")
    #otuprep.levels <- c("t100","r100","r500","r1000")
    #otuprep.labels <- c("raw","r100","r500","r1000")
    otuprep.levels <- c("c100","r100","r500","r1000")
    otuprep.labels <- c("CLR","r100","r500","r1000")
    
    otuprep <- list(levels = otuprep.levels, labels = otuprep.labels)
    
    plotpars<-list(logLoc=logLoc, species4=species4, topo=topo,
                   overlapSp=overlapSp, overlapSp6=overlapSpColors6,
                   species09=species09, species11=species11,
                   yearSubset = yearSubset,
                   traits = traits,
                   otuprep = otuprep,
                   orderSp = orderSp.df)
    #binom.levels=binom.levels$label,
    
    return(plotpars)
}

make_ggplot_theme <- function(){
  
  require(ggplot2)
  
  #my ggplot template
  mytheme <- theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(panel.border = element_rect(colour = "black"),      #put a black box around the plotting area
          axis.line = element_line(colour = "black"),                 #axis lines are in black
          panel.grid.major = element_blank(),                         #turn off the gridlines
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(face='bold.italic'),         #facet labels hjust=0.05
          strip.text.y = element_text(face='bold.italic'),
          strip.background = element_rect(fill = 'white', colour='black'),    
          legend.key = element_blank(),                               #turn off box around legend
          plot.title=element_text(hjust=0, vjust=0.5, face='bold'), #style and position of the panel label
          plot.margin = unit(c(0.05,0.05,0.05,0.05),"in")
    )
  
  return(mytheme)
}


# study design #-----------------------------------------------#

make_sample_tileplot <- function(samples){
  
  params <- plottingParams_rotplotNobags()
  params$yearSubset
  # tile plot of samples (no mush)
  samples %>%
    filter(logLoc != "mush") %>%
    filter(gene == "16S") %>% # look at bacteria -- none of these samples were dropped due to sequencing fails
    mutate(dep.harv = factor(dep.harv, 
                             levels = params$yearSubset$levels2.new, 
                             labels = params$yearSubset$labels.new)) %>%
    mutate(topo = factor(topo,
                         levels = params$topo$levels,
                         labels = params$topo$labels)) %>%
    mutate(logLoc = factor(logLoc,
                           levels = params$logLoc$levels,
                           labels = params$logLoc$labels)) %>%
    mutate(logloc.topo = paste(logLoc, topo, sep = " & ")) -> samples
  
  #summarize into number of rep samples per treatment
  samples %>%
    group_by(gene, logloc.topo, dep.harv, species4) %>%
    summarize(n = length(table)) -> samples.summ
  
  #add the bionomial and gymno/angio that matches each species4
  samples.summ %>%
    left_join(params$orderSp) -> samples.summ 
    
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(data = samples.summ, aes(x = logloc.topo, y = reorder(binomial, -order), 
                                       label = n)) +
    geom_tile(fill = "white", color = "black") +
    geom_text() +
    facet_grid(~dep.harv) +
    theme_bw() +
    xlab("Within stem & watershed position") + ylab("Wood species") +
    scale_fill_distiller(name = "Samples", direction = 2) +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  return(p)
  # each sample was processed to detect ITS and 16S amplicons
  
}

woodspecies_table <- function(samples){
  
  params <- plottingParams_rotplotNobags()
  samples %>%
    select(species4, binomial, family, angio.gymno) %>%
    mutate(binomial = sub("_"," ", binomial)) -> sp.tab
  sp.tab <- unique(sp.tab)
  sp.tab$species4 <- factor(sp.tab$species4, levels = params$orderSp)
  sp.tab %>%
    arrange(species4) %>%
    rename('Code'='species4',
           'Genus species' = 'binomial',
           'Family' = 'family',
           'Angiosperm / Gymnosperm' = 'angio.gymno') -> sp.tab
  
  return(sp.tab)
  
}

make_mush_tileplot <- function(samples.missing){
  
  # missing samples -- too decayed and lost
  samples.missing %>%
    mutate(failType = ifelse(mushfail == TRUE & seqfail == FALSE, "Too decayed", "Unknown")) %>%
    mutate(failType = ifelse(mushfail == FALSE & seqfail == TRUE, "Sequencing fail", failType)) %>%
    mutate(dep.harv = paste(yrdeploy, yrharv, sep = "_")) %>%
    filter(failType != "Sequencing fail") %>%
    filter(gene == "16S") %>%
    filter(logLoc == "t") %>%
    group_by(topo, dep.harv, species4, failType) %>%
    summarize(n = length(topo)) -> summ.missing.beforeseq
  #order factors
  params <- plottingParams_rotplotNobags()
  uniq09Sp <- params$species09$levels[!params$species09$levels %in% params$species11$levels]
  uniq11Sp <- params$species11$levels[!params$species11$levels %in% params$species09$levels]
  orderSp <- c(uniq09Sp, params$overlapSp$levels, uniq11Sp)
  summ.missing.beforeseq$dep.harv <- factor(summ.missing.beforeseq$dep.harv, 
                                            levels = c("2011_2012","2011_2014","2009_2012","2009_2014"), 
                                            labels = params$yearSubset$labels)
  summ.missing.beforeseq$species4 <- factor(summ.missing.beforeseq$species4, levels = orderSp)
  #plot
  p <- ggplot(data = summ.missing.beforeseq, aes(x = topo, y = species4, label = n, fill = failType)) +
    geom_tile() +
    geom_text() +
    facet_grid(~dep.harv) +
    theme_bw() +
    xlab("Position within watershed and log") + ylab("Wood species") +
    scale_fill_brewer(name = "Reason missing", type = "qual")
  
  return(p)
  
}

make_seqfail_tileplot <- function(samples.missing){
  
  # missing samples -- sequencing failures
  samples.missing %>%
    filter(seqfail == TRUE) %>%
    mutate(dep.harv = paste(yrdeploy, yrharv, sep = "_")) %>%
    mutate(topo.logloc = paste(topo, logLoc, sep = "_")) %>%
    group_by(gene, topo.logloc, dep.harv, species4) %>%
    summarize(n = length(gene)) -> summ.missing.afterseq
  #order factors
  params <- plottingParams_rotplotNobags()
  uniq09Sp <- params$species09$levels[!params$species09$levels %in% params$species11$levels]
  uniq11Sp <- params$species11$levels[!params$species11$levels %in% params$species09$levels]
  orderSp <- c(uniq09Sp, params$overlapSp$levels, uniq11Sp)
  summ.missing.afterseq$dep.harv <- factor(summ.missing.afterseq$dep.harv, 
                                           levels = c("2011_2012","2011_2014","2009_2012","2009_2014"), 
                                           labels = params$yearSubset$labels)
  summ.missing.afterseq$species4 <- factor(summ.missing.afterseq$species4, levels = orderSp)
  samples.summ$topo.logloc <- factor(samples.summ$topo.logloc, levels = c("H_t","H_b","L_t","L_b"))
  #plot
  p <- ggplot(data = summ.missing.afterseq, aes(x = topo.logloc, y = species4, label = n, fill = 1)) +
    geom_tile() +
    geom_text() +
    facet_grid(~dep.harv) +
    theme_bw() +
    xlab("Position within watershed and log") + ylab("Wood species") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill = FALSE)
  
  return(p)
  
}


# community composition ordinations #-----------------------------------------------#

makeDF_ordplot_overlap <- function(data, presAbs){
  
  # do the ordination
  if(presAbs == TRUE){
    mod.obj <- capscale(data$otu ~ 1, data = NULL, distance = "jaccard", binary = TRUE)
    }
  if(presAbs == FALSE){
    mod.obj <- capscale(data$otu ~ 1, data = NULL, distance = "euclidean")
  }
  
  # do the environmental fit
  env <- data$meta
  env %>%
    mutate(grp = paste(topo, species4, yrdeploy, sep = "_")) -> env
  envfit.obj <- envfit(mod.obj ~ grp, env)
  
  # make a point dataframe
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  scrs <- cbind(scrs, tidySeq_sampleName = row.names(scrs))
  scrs %>%
    left_join(env) -> scrs
  
  # arrows that connect centroids
  centr.df <- data.frame(grp1 = row.names(envfit.obj$factors$centroids), envfit.obj$factors$centroids)
  centr.df %>%
    separate(grp1, into = c("drop","keep"), sep ="grp") %>%
    separate(keep, into = c("topo","species4","yrdeploy")) %>%
    select(-drop) -> centr.df.tmp
  centr.df.tmp %>%
    select(-MDS2) %>%
    spread(key = yrdeploy, value = MDS1) %>%
    rename('MDS1_2009' = `2009`,
           'MDS1_2011' = `2011`) -> mds1.df
  centr.df.tmp %>%
    select(-MDS1) %>%
    spread(key = yrdeploy, value = MDS2) %>%
    rename('MDS2_2009' = `2009`,
           'MDS2_2011' = `2011`) %>%
    left_join(mds1.df) -> centroids
  
  # new cols
  scrs %>%
    mutate(species4topo = paste(species4, topo, sep = "_")) -> scrs
  centroids %>%
    mutate(species4topo = paste(species4, topo, sep = "_")) %>%
    arrange(desc(topo), species4) -> centroids
  
  plot.data <- list(scrs = scrs, centroids = centroids,
                    mod.obj = mod.obj, envfit.obj = envfit.obj)
  return(plot.data)
  
}

ordplot_overlap <- function(plot.data.overlap, output.path, no_groups,
                            xlims, ylims){
  
  #plot.data.overlap
  #output.path
  #no_groups = T
  
  if(no_groups == TRUE){
    condition <- !grepl("sapros", names(plot.data.overlap)) & !grepl("basidios", names(plot.data.overlap))
    plot.data.overlap <- plot.data.overlap[condition]
  }
  
  # colors and linetypes
  if(no_groups == TRUE){
    panelLabels <- letters[1:length(plot.data.overlap)]
  }else{
    panelLabels <- names(plot.data.overlap)
  }
  pars <- plottingParams_rotplotNobags()
  species4topo_colors <- pars[['overlapSp6']]
  COLORS <- c(species4topo_colors[['CEOC_L']], species4topo_colors[['CEOC_H']], species4topo_colors[['CEOC_L']], species4topo_colors[['CEOC_H']],
              species4topo_colors[['JUVI_L']], species4topo_colors[['JUVI_H']], species4topo_colors[['JUVI_L']], species4topo_colors[['JUVI_H']],
              species4topo_colors[['QUVE_L']], species4topo_colors[['QUVE_H']], species4topo_colors[['QUVE_L']], species4topo_colors[['QUVE_H']])
  LINETYPES <- c(1, 1, 2, 2, 
                 1, 1, 2, 2,
                 1, 1, 2, 2)
  
  if(no_groups == TRUE){
    pdf(file = paste0(output.path, "ordplots_overlap.pdf"), width = 6, height = 6)
    par(mfrow=c(2,2), mar = c(0, 0, 0, 0), 
        oma = c(4, 4, 2, 2), tcl = -0.25, mgp = c(2, 0.6, 0))
  }else{
    pdf(file = paste0(output.path, "ordplots_overlap_allgroups.pdf"), width = 6, height = 6)
    par(mfrow=c(4,2), mar = c(0, 0, 0, 0), 
        oma = c(4, 4, 2, 2), tcl = -0.25, mgp = c(2, 0.6, 0))
  }
  #p<-1
  for(p in 1:length(panelLabels)){
    
    scrs <- plot.data.overlap[[p]]$scrs
    centroids <- plot.data.overlap[[p]]$centroids
    mod.obj <- plot.data.overlap[[p]]$mod.obj
    envfit.obj <- plot.data.overlap[[p]]$envfit.obj
    
    # plot
    with(scrs, plot(x = MDS1, y = MDS2, 
                    type = "n", axes = FALSE, asp = 1, 
                    xlab = " ", ylab= " ",
                    xlim = xlims, ylim = ylims))
    box()
    
    if(no_groups == TRUE){
      if(p %in% c(3,4)){
        axis(1)
      }
      if(p %in% c(1,3)){
        axis(2)
      }
    }
    
    # add ellipses
    for(i in 1:length(unique(scrs$grp))){
      with(envfit.obj, ordiellipse(mod.obj, scrs$grp, 
                                   kind="se", conf=0.95, lwd = 2,
                                   lty=LINETYPES[i], 
                                   col=COLORS[i], 
                                   show.groups = unique(scrs$grp)[i], 
                                   label=FALSE))
    }
    # add centroid arrows to connect groups through decay time
    with(centroids, 
         arrows(x0 = MDS1_2011, y0 = MDS2_2011, 
                x1 = MDS1_2009, y1 = MDS2_2009,
                col = species4topo_colors, lwd = 2))
    
    #add panel labels
    mtext(panelLabels[p], side = 3, line = -1, adj = 0.05, cex = 0.8, font = 2)
    
  }
  mtext("PCoA 1", side = 1, outer = TRUE, cex = 1, line = 2.2)
  mtext("PCoA 2", side = 2, outer = TRUE, cex = 1, line = 2.2)
  
  if(no_groups == TRUE){
    mtext("Bacteria                                         Fungi", 
          side=4, outer=TRUE, cex=1, line=1, font=2)
    mtext("2012                                              2014", 
          side=3, outer=TRUE, cex=1, line=1, font=2)
  }
  
  dev.off()
  
  #make legend separately
  st.df <- data.frame(species4topo = names(species4topo_colors),
                      colors = species4topo_colors, stringsAsFactors = FALSE)
  st.df %>%
    separate(species4topo, into=c("species4","topo"), remove = FALSE) %>%
    arrange(species4) -> st.df
  
  pdf(paste0(output.path, "ordplots_overlap_legend.pdf"), width = 5, height = 5)
  par(mfrow = c(1,1))
  plot(x = 1:10, y = 1:10, type = "n", axes = FALSE, xlab = "", ylab ="")
  legend(1, 10, legend = st.df$species4topo, col = st.df$colors, 
         lty = 1, lwd = 2, title = "Wood species & Log position")
  legend(1, 4, legend = c("2011", "2009"), col = 1, 
         lty = c(1, 2), lwd = 2, title = "Year deployed")
  dev.off()

  
}

makeDF_ordplot_deploy <- function(data, presAbs, GRPname){
  
  params <- plottingParams_rotplotNobags()
  data$meta %>%
    mutate(topo = factor(topo,
                         levels = params$topo$levels,
                         labels = params$topo$labels)) %>%
    mutate(logLoc = factor(logLoc,
                           levels = params$logLoc$levels,
                           labels = params$logLoc$labels)) %>%
    mutate(angio.gymno = tolower(angio.gymno)) -> data$meta
  
  # do the ordination
  if(presAbs == TRUE){
    mod.obj <- capscale(data$otu ~ 1, data = NULL, distance = "jaccard", binary = TRUE)
  }
  if(presAbs == FALSE){
    mod.obj <- capscale(data$otu ~ 1, data = NULL, distance = "euclidean")
  }
  
  # do the environmental fit
  env <- data$meta
  env %>%
    rename_("GRP_"= GRPname) -> env
  envfit.obj <- envfit(mod.obj ~ GRP_, env)
  
  # make a point dataframe
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  scrs <- cbind(scrs, tidySeq_sampleName = row.names(scrs))
  scrs %>%
    left_join(env) -> scrs
  
  plot.data <- list(scrs = scrs, mod.obj = mod.obj, envfit.obj = envfit.obj)
  return(plot.data)
  
}

ordplot_deploy<-function(plot.data.deploy.x, output.path, GRPname, no_groups, 
                         xlims, ylims){
  
  if(no_groups == TRUE){
    condition <- !grepl("sapros", names(plot.data.deploy.x)) & !grepl("basidios", names(plot.data.deploy.x))
    plot.data.deploy.x <- plot.data.deploy.x[condition]
    
    # re-order the data subsets
    names(plot.data.deploy.x)
    first <- grepl("d11", names(plot.data.deploy.x)) & grepl("ITS", names(plot.data.deploy.x))
    second <- grepl("d09", names(plot.data.deploy.x)) & grepl("ITS", names(plot.data.deploy.x))
    third<- grepl("d11", names(plot.data.deploy.x)) & grepl("16S", names(plot.data.deploy.x))
    fourth <- grepl("d09", names(plot.data.deploy.x)) & grepl("16S", names(plot.data.deploy.x))
    plot.data.deploy.x <- c(plot.data.deploy.x[first], plot.data.deploy.x[second],
                     plot.data.deploy.x[third], plot.data.deploy.x[fourth])
  }
  
  
  # colors and linetypes
  panelLabels <- letters[1:length(plot.data.deploy.x)]
  
  # gray, black
  params <- plottingParams_rotplotNobags()
  COLORS <- params$topo$colors
  
  pdf(file = paste0(output.path, "ordplots_deploy_", GRPname, ".pdf"), width = 8, height = 5)
  par(mfrow=c(2,4), mar = c(0, 0, 0, 0), 
      oma = c(4, 4, 2, 2), tcl = -0.25, mgp = c(2, 0.6, 0))
  for(p in 1:length(panelLabels)){
    
    scrs <- plot.data.deploy.x[[p]]$scrs
    mod.obj <- plot.data.deploy.x[[p]]$mod.obj
    envfit.obj <- plot.data.deploy.x[[p]]$envfit.obj
    curr.plotdata <- names(plot.data.deploy.x)[[p]]
    
    #if(GRPname == "logLoc"){
    #  scrs %>%
    #    filter(GRP_ != "mush") -> scrs
    #}
    
    # plot
    with(scrs, plot(x = MDS1, y = MDS2, 
                    type = "n", axes = FALSE, asp = 1, 
                    xlab = " ", ylab= " ",
                    xlim = xlims, ylim = ylims))
    box()
    if(p %in% c(5:8)){
      axis(1)
    }
    if(p %in% c(1,5)){
      axis(2)
    }
    # add ellipses
    for(i in 1:length(unique(scrs$GRP_))){
      with(envfit.obj, ordiellipse(mod.obj, scrs$GRP_, 
                                   kind="se", conf=0.95, lwd = 2,
                                   col=COLORS[i], 
                                   show.groups = unique(scrs$GRP_)[i], 
                                   label=TRUE))
    }
    #mtext(curr.plotdata, side = 3, line = -1, adj = 0.05, cex = 0.8, font = 2)
  }
  mtext("PCoA 1", side = 1, outer = TRUE, cex = 1, line = 2.2)
  mtext("PCoA 2", side = 2, outer = TRUE, cex = 1, line = 2.2)
  mtext("Bacteria                                         Fungi", 
        side=4, outer=TRUE, cex=1, line=0.5, font=2)
  
  ta <- "1 yr (2011-12)"
  tb <- "3 yrs (2011-14)"
  tc <- "3 yrs (2009-12)"
  td <- "5 yrs (2009-14)"
  mtext(paste(ta, tb, tc, td, sep="             "),
        side=3, outer=TRUE, cex=1, line=0.5, font=2)
  
  dev.off()
  
}


# coef plot of richness ~ study factors #-----------------------------------------------#

plot_richness_compareOTUprep <- function(rich.df, datasetType){
  
  mytheme <- make_ggplot_theme()
  if(datasetType == "overlap"){
    
    p <- ggplot(data = rich.df, aes(y = term, x = Estimate, alpha =  signif, color = otuprep)) +
      geom_point() +
      geom_errorbarh(aes(xmin = Estimate - Std..Error, xmax = Estimate + Std..Error)) +
      facet_wrap(gene~year, scales = "free_x") +
      geom_vline(xintercept = 0, linetype = 2) + 
      mytheme
  }
  
  if(datasetType == "deploy"){
    # deployed in 2009
    rich.df %>%
      filter(subset == "d09") -> tmp
    p.d09 <- ggplot(data = tmp, aes(y = term, x = Estimate, alpha =  signif, color = otuprep)) +
      geom_point() +
      geom_errorbarh(aes(xmin = Estimate - Std..Error, xmax = Estimate + Std..Error)) +
      facet_wrap(gene~year, scales = "free_x") +
      geom_vline(xintercept = 0, linetype = 2) + 
      mytheme
    # deployed in 2011
    rich.df %>%
      filter(subset == "d11") -> tmp
    p.d11 <- ggplot(data = tmp, aes(y = term, x = Estimate, alpha =  signif, color = otuprep)) +
      geom_point() +
      geom_errorbarh(aes(xmin = Estimate - Std..Error, xmax = Estimate + Std..Error)) +
      facet_wrap(gene~year, scales = "free_x") +
      geom_vline(xintercept = 0, linetype = 2) + 
      mytheme
    p <- list(p.d09 = p.d09, p.d11 = p.d11)
  }
  
  return(p)
}


# barplot of r2 and centroid distances #-----------------------------------------------#

preprocess_varpart_rd.traits <- function(rd.traits){
  df <- rd.traits$varpart.df
  # missing the pval column. why?
  df %>%
    separate(subset, into = c("yrharv","genething"), sep =4, remove = F) %>%
    separate(genething, into = c("gene","yrdeploy")) -> df
  rd.traits$varpart.df <- df
  
  return(rd.traits)
}

dbrda_order <- function(df, subsettype, trait_tf){
  
  if(subsettype == "o"){
    
    df %>%
      mutate(terms = recode_factor(terms,  yrdeploy = "Year deployed",
                                   logLoc = "Within stem",
                                   topo = "Within watershed",
                                   species4 = "Wood species", .ordered = T)) %>%
      mutate(otuprep = recode_factor(otuprep, 
                                     c100 = "CLR",
                                     r100 = "r100",
                                     r500 = "r500",
                                     r1000 = "r1000", .ordered = T)) %>%
      mutate(gene = recode_factor(gene, ITS = "Fungi", 
                                  `16S` = "Bacteria", .ordered = T)) %>%
      mutate(community = paste(gene, group, sep = "_")) %>%
      mutate(community = recode_factor(community, 
                                       Fungi_NA = "Fungi",
                                       Bacteria_NA = "Bacteria",
                                       Fungi_sapros = "Saprotrophs",
                                       Fungi_basidios = "Basidiomycetes")) -> df
    
  }
  
  if(subsettype == "d"){
    df %>%
      mutate(terms = recode_factor(terms,
                                   logLoc = "Within stem",
                                   topo = "Within watershed",
                                   species4 = "Wood species",  .ordered = T)) %>%
      mutate(gene = recode_factor(gene, ITS = "Fungi", 
                                  `16S` = "Bacteria", .ordered = T)) %>%
      mutate(otuprep = recode_factor(otuprep, 
                                     "c100" = "CLR",
                                     "r100" = "r100",
                                     "r500" = "r500",
                                     "r1000" = "r1000", .ordered = T)) -> df
    if(trait_tf == F){
      df %>%
        select(-group) %>%
        separate(subset, into = c("thing1","thing2","group"), fill = "right") %>%
        mutate(community = paste(gene, group, sep = "_")) %>%
        mutate(community = recode_factor(community, 
                                         Fungi_NA = "Fungi",
                                         Bacteria_NA = "Bacteria",
                                         Fungi_sapros = "Saprotrophs",
                                         Fungi_basidios = "Basidiomycetes")) %>%
        mutate(orig.subset = paste(thing1, thing2, sep = "_")) %>%
        mutate(subset = recode_factor(orig.subset, 
                                      
                                      `201416S_d09` = "5 yrs (2009-14)",
                                      `2014ITS_d09` = "5 yrs (2009-14)",
                                      
                                      `201216S_d09` = "3 yrs (2009-12)",
                                      `2012ITS_d09` = "3 yrs (2009-12)",
                                      
                                      `201416S_d11` = "3 yrs (2011-14)",
                                      `2014ITS_d11` = "3 yrs (2011-14)",
                                      
                                      `201216S_d11` = "1 yr (2011-12)",
                                      `2012ITS_d11` = "1 yr (2011-12)", .ordered = T)) -> df
    }
    if(trait_tf == T){
      df %>%
        mutate(subset = recode_factor(subset, 
                                      
                                      `201416S_d09` = "5 yrs (2009-14)",
                                      `2014ITS_d09` = "5 yrs (2009-14)",
                                      
                                      `201216S_d09` = "3 yrs (2009-12)",
                                      `2012ITS_d09` = "3 yrs (2009-12)",
                                      
                                      `201416S_d11` = "3 yrs (2011-14)",
                                      `2014ITS_d11` = "3 yrs (2011-14)",
                                      
                                      `201216S_d11` = "1 yr (2011-12)",
                                      `2012ITS_d11` = "1 yr (2011-12)", .ordered = T)) -> df
    }
    
    }
  
  return(df)

}

plot_varpartbars_compareOTUprep <- function(dbrda.intermed, subsettype, trait_tf){

  # dbrda.intermed = rd.traits 
  # subsettype = "d"
  # trait_tf = T
  
  mytheme <- make_ggplot_theme()
  params <- plottingParams_rotplotNobags()
  color.vec <- params$topo$colors # c("darkgray", "black") # for subsettype == "o"
  names(color.vec)<- NULL
  
  # colors
  require(RColorBrewer)
  light.red <- brewer.pal(n=8, name = "RdBu")[2] # 2011-12 (1yr)
  dark.red <- brewer.pal(n=8, name = "RdBu")[1] # 2011-2014 (3yrs)
  
  light.blue <- brewer.pal(n=8, name = "RdBu")[7] # 2009-2012 (3yrs)
  dark.blue <- brewer.pal(n=8, name = "RdBu")[8] # 2009-2014 (5yrs)
  deploy.colors <- c(dark.blue, light.blue, dark.red, light.red)
  
  light.purple <- brewer.pal(n=4, name = "Purples")[3] # overlap harvested in 2012
  dark.purple <- brewer.pal(n=4, name = "Purples")[4] # overlap harvested in 2014
  overlap.colors <- c(light.purple, dark.purple)
  
  ### r2
  # (1) clean and order factors
  varpart.df <- dbrda.intermed$varpart.df
  varpart.df <- dbrda_order(df = varpart.df, subsettype = subsettype, trait_tf = trait_tf)
  
  if(trait_tf == T){
    varpart.df %>%
      filter(otuprep == "CLR") -> varpart.df.raw
  }else{
    varpart.df %>%
      filter(otuprep == "CLR") %>%
      filter(is.na(group)) -> varpart.df.raw
  }
  
  # (2) plot w/ rarefaction sensitivity
  if(subsettype == "o"){
    p.r2 <- ggplot(data = varpart.df, 
                   aes(y = Adj.R.square, x = terms, fill = yrharv)) +
      geom_bar(stat = "identity", position = "dodge") + 
      facet_grid(otuprep~community) + 
      coord_flip() + 
      mytheme +
      xlab("") + ylab("Adjusted R squared") +
      scale_fill_manual(name = "Year harvested", 
                        values = overlap.colors)
    p.r2
  }
  if(subsettype == "d"){
    if(trait_tf == T){
      p.r2 <- ggplot(data = varpart.df, 
                     aes(y = Adj.R.square, 
                         x = terms, fill = subset)) +
        geom_bar(stat = "identity", position = "dodge") + 
        facet_grid(otuprep~gene) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Adjusted R squared") +
        scale_fill_manual(name = "Deployment subset",
                          values = deploy.colors)
      p.r2
      
    }else{
      p.r2 <- ggplot(data = varpart.df, 
                     aes(y = Adj.R.square, 
                         x = terms, fill = subset)) +
        geom_bar(stat = "identity", position = "dodge") + 
        facet_grid(otuprep~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Adjusted R squared") +
        scale_fill_manual(name = "Deployment subset",
                          values = deploy.colors)
      p.r2
    }
    
  }
  
  # (3) plot raw OTU only
  if(subsettype == "o"){
    
    p.r2.raw <- ggplot(data = varpart.df.raw, 
                 aes(y = Adj.R.square, x = terms)) +
      geom_bar(stat = "identity", position = "dodge") + 
      facet_grid(yrharv~community) + 
      coord_flip() + 
      mytheme +
      xlab("") + ylab("Adjusted R squared")
    p.r2.raw
  }
  if(subsettype == "d"){
    
    varpart.df.raw %>%
      mutate(subset = factor(subset, 
                             levels = c("1 yr (2011-12)","3 yrs (2009-12)",
                                        "3 yrs (2011-14)","5 yrs (2009-14)"))) -> varpart.df.raw
    if(trait_tf == T){
      p.r2.raw <- ggplot(data = varpart.df.raw, 
                         aes(y = Adj.R.square, x = terms)) +
        geom_bar(stat = "identity", position = "dodge") + 
        facet_grid(subset~gene) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Adjusted R squared")
      p.r2.raw
      
    }else{
      p.r2.raw <- ggplot(data = varpart.df.raw, 
                         aes(y = Adj.R.square, x = terms)) +
        geom_bar(stat = "identity", position = "dodge") + 
        facet_grid(subset~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Adjusted R squared")
      p.r2.raw
      
    }
    
    
  }
  
  if(trait_tf == T){
    p.dist = NULL
    p.dist.raw = NULL
    
  }else{
    ### centroid distance
    # (1) clean and order factors
    dist.df <- dbrda.intermed$dist.df
    dist.df <- dbrda_order(df = dist.df, subsettype = subsettype, trait_tf = trait_tf)
    dist.df %>%
      filter(otuprep == "CLR") %>%
      filter(is.na(group)) -> dist.df.raw
    
    # (2) plot w/ rarefaction sensitivity
    if(subsettype == "o"){
      p.dist <- ggplot(data = dist.df, 
                       aes(y = mean.distall, 
                           x = terms, 
                           fill = yrharv)) +
        geom_bar(stat = "identity", 
                 position = "dodge") + 
        geom_errorbar(aes(ymin = mean.distall - se.distall,
                          ymax = mean.distall + se.distall), 
                      position = "dodge") +
        facet_grid(otuprep~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Centroid distance") +
        scale_fill_manual(name = "Year harvested", 
                          values = overlap.colors)
      p.dist
      
    }
    if(subsettype == "d"){
      p.dist <- ggplot(data = dist.df, 
                       aes(y = mean.distall, 
                           x = terms, 
                           fill = subset)) +
        geom_bar(stat = "identity", 
                 position = "dodge") + 
        geom_errorbar(aes(ymin = mean.distall - se.distall,
                          ymax = mean.distall + se.distall), 
                      position = "dodge") +
        facet_grid(otuprep~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Centroid distance") +
        scale_fill_manual(name = "Deployment subset",
                          values = deploy.colors)
      p.dist
    }
    
    # (3) plot raw OTU only
    if(subsettype == "o"){
      p.dist.raw <- ggplot(data = dist.df.raw, 
                           aes(y = mean.distall, 
                               x = terms)) +
        geom_bar(stat = "identity", 
                 position = "dodge") + 
        geom_errorbar(aes(ymin = mean.distall - se.distall,
                          ymax = mean.distall + se.distall), 
                      position = "dodge") +
        facet_grid(yrharv~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Centroid distance") 
      p.dist.raw
    }
    if(subsettype == "d"){
      dist.df.raw %>%
        mutate(subset = factor(subset, 
                               levels = c("1 yr (2011-12)","3 yrs (2009-12)",
                                          "3 yrs (2011-14)","5 yrs (2009-14)"))) -> dist.df.raw
      p.dist.raw <- ggplot(data = dist.df.raw, 
                           aes(y = mean.distall, x = terms)) +
        geom_bar(stat = "identity", 
                 position = "dodge") + 
        geom_errorbar(aes(ymin = mean.distall - se.distall,
                          ymax = mean.distall + se.distall), 
                      position = "dodge") +
        facet_grid(subset~community) + 
        coord_flip() + 
        mytheme +
        xlab("") + ylab("Centroid distance")
      p.dist.raw
    }
    
  }
  
  
  result <- list(p.r2 = p.r2,
                 p.dist = p.dist,
                 p.r2.raw = p.r2.raw,
                 p.dist.raw = p.dist.raw)
  

  return(result)
}


# wood species trait space ordination #-----------------------------------------------#

make_trait_pca <- function(tmp, indx.logtrt){
  
  # do PCA and take the first 2 axes
  tmp <- tmp[complete.cases(tmp),]
  tmp1 <- apply(tmp[,-1], 2, scale)
  row.names(tmp1) <- tmp$species4
  pc <- prcomp(tmp1)
  df <- data.frame(species4 = row.names(pc$x), pc$x[,1:2])
  df %>%
    left_join(indx.logtrt) -> df
  
  #put together the loadings dataframe
  datapc <- data.frame(varnames=rownames(pc$rotation), pc$rotation[,1:2])
  mult <- min(
    (max(df[,"PC1"]) - min(df[,"PC2"])/(max(datapc[,"PC2"])-min(datapc[,"PC2"]))),
    (max(df[,"PC2"]) - min(df[,"PC1"])/(max(datapc[,"PC1"])-min(datapc[,"PC1"])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * datapc$PC1,
                      v2 = .7 * mult * datapc$PC2)
  
  #plot
  summary(pc)
  mytheme <- make_ggplot_theme()
  params <- plottingParams_rotplotNobags()
  color.vec <- c( "blue", "black")
  p <- ggplot(data = df, aes(x = PC1, y = PC2, label = species4, color = angio.gymno)) +
    geom_text(size = 3) +
    geom_segment(data = datapc, 
                 aes(x = 0, y = 0, xend = v1, yend = v2), 
                 arrow = arrow(length=unit(0.2,"cm")), color = "darkgray",
                 inherit.aes = F) + 
    geom_text(data = datapc,
              aes(x = v1, y = v2, label = varnames), 
              size = 3, color = "darkgray", hjust = "outward", vjust = "outward") + 
    coord_equal() +
    mytheme +
    scale_color_manual(values = color.vec)
  p.list <- list(p=p, pc=pc)
  
  return(p.list)
}

plot_woodTraits <- function(woodTraits){
  
  # create species and binomial indx
  logTrt <- load_logTrt()
  logTrt %>%
    filter(angio.gymno %in% c("Angiosperm","Gymnosperm")) -> logTrt
  indx.logtrt <- unique(logTrt[,c("species4","family","binomial","angio.gymno")])
  
  # subset by deployment
  params <- plottingParams_rotplotNobags()
  indx.logtrt %>%
    mutate(deploy09 = ifelse(species4 %in% params$species09$levels, TRUE, FALSE)) %>%
    mutate(deploy11 = ifelse(species4 %in% params$species11$levels, TRUE, FALSE)) -> indx.logtrt
  # 09
  indx.logtrt %>%
    filter(deploy09 == TRUE) -> indx.09
  woodTraits %>%
    filter(species4 %in% indx.09$species4) -> tmp.09
  # 11
  indx.logtrt %>%
    filter(deploy11 == TRUE) -> indx.11
  woodTraits %>%
    filter(species4 %in% indx.11$species4) -> tmp.11
  
  # trait correlations
  require(corrplot)
  tmp.09 %>%
    select(-species4) -> tmp1.09
  cor.09 <- corrplot(cor(tmp1.09))
  tmp.11 %>%
    select(-species4) -> tmp1.11
  cor.11 <- corrplot(cor(tmp1.11))
  
  # pca plot
  pca <- make_trait_pca(tmp = tmp.09, indx.logtrt = indx.logtrt)
  summary(pca$pc)
  p.09 <- pca$p + xlab("PC1 (37%)") + ylab("PC2 (20%)")
  
  pca <- make_trait_pca(tmp = tmp.11, indx.logtrt = indx.logtrt)
  summary(pca$pc)
  p.11 <- pca$p + xlab("PC1 (52%)") + ylab("PC2 (16%)")
  
  result <- list(p.09 = p.09, p.11 = p.11, cor.09 = cor.09, cor.11 = cor.11)
  return(result)
}

# community composition constrained by wood trait environment #-----------------------------------------------#

dbrdaTest_woodTraits_finalmod_plots <- function(rd.all, woodTraits, data.deploy){
  
  table.names <- names(data.deploy)
  trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
  for(i in 1:length(data.deploy)){
    
    terms <- extract_rhs.terms(rd.all = rd.all, table.name = table.names[i])
    mat.otu.pres <- data.deploy[[i]][["otu"]]
    envVars <- data.deploy[[i]][["meta"]]
    
    # trait data for species4 == PLOC is missing (only in d09), so exclude it
    envVars %>%
      filter(species4 != "PLOC") -> envVars
    mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
    
    # set up the model
    envVars %>%
      left_join(woodTraits) %>%
      select(terms) -> envVars.to
    cap.env <-  capscale(mat.otu.pres ~ ., data = envVars.to, distance = 'jaccard', binary = TRUE)
    
    # save plot
    params <- plottingParams_rotplotNobags()
    if(grepl("d11", table.names[i])){
      envVars$species4 <- factor(envVars$species4, levels = params$species11$labels, labels = params$species11$colors)
      sp.legend <- params$species11$labels
      sp.legend.col <- params$species11$colors
    }else{
      envVars$species4 <- factor(envVars$species4, levels = params$species09$labels, labels = params$species09$colors)
      sp.legend <- params$species11$labels
      sp.legend.col <- params$species11$colors
    }
    envVars %>%
      mutate(shape = NA) %>%
      mutate(shape = ifelse(logLoc == "t" & topo == "H", 24, shape)) %>%
      mutate(shape = ifelse(logLoc == "t" & topo == "L", 2, shape)) %>%
      mutate(shape = ifelse(logLoc == "b" & topo == "H", 25, shape)) %>%
      mutate(shape = ifelse(logLoc == "b" & topo == "L", 6, shape)) -> envVars
    plot.legend <- c("tH","tL","bH","bL")
    plot.legend.pch <- c(24,2,25,6)
    
    #plot
    file.name <- paste0(output.path, "dbrda_woodtraits_",table.names[i],".pdf")
    pdf(file = file.name, width = 6, height = 6)
    plot(cap.env, type = "n")
    points(cap.env, display = "sites", 
           col = as.character(envVars$species4), pch = envVars$shape, bg = as.character(envVars$species4))
    points(cap.env, display = "bp")
    text(cap.env, display = "bp")
    #legend("topright", legend = sp.legend, col = sp.legend.col, pch = 16)
    title(main = table.names[i])
    dev.off()
    
  }
  
}




