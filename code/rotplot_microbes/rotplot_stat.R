
# taq tests #-----------------------------------------------#

test_taqdiffs <- function(df.r, yLab){
  
  # taqType differences (yrseq == 2017)
  df.r %>%
    filter(yrseq == 2017) -> df.r.17
  
  # for yrharv == 2012
  df.r.17 %>%
    filter(yrharv == 2012) -> df.r.17.12
  mod.12 <- lmer(resp ~ taqType + (1|sampid), data = df.r.17.12)
  
  # for yrharv == 2014
  df.r.17 %>%
    filter(yrharv == 2014) -> df.r.17.14
  mod.14 <- lmer(resp ~ taqType + (1|sampid), data = df.r.17.14)
  
  #plot
  mytheme <- make_ggplot_theme()
  p.taq.17 <- ggplot(df.r.17, aes(x = taqType, y = resp, group = sampid)) +
    geom_point() +geom_line() +
    facet_grid(~yrharv) +
    xlab("Taq type") + ylab(yLab) + mytheme +
    ggtitle("Samples re-run in 2017")
  #p.taq.17
  #ggsave(p.taq.17, file = paste0(output.path, "sampreads_byTaq.pdf"), width = 6, height = 3)
  
  result <- list(mod.12 = mod.12, 
                 mod.14 = mod.14,
                 p = p.taq.17)
  
}

test_yrseq_gotaq <- function(df.r, yLab){
  
  df.r %>%
    filter(yrharv == 2012) %>%
    filter(taqType == "GoTaq") -> df.r.go
  
  df.r.go %>%
    group_by(sampid) %>%
    summarize(n = length(yrharv)) %>%
    filter(n != 2) -> exclude.sampids
  
  df.r.go %>%
    filter(!sampid %in% c(exclude.sampids$sampid)) -> df.r.go
  
  mod <- lmer(resp ~ yrseq + (1|sampid), data = df.r.go)
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(df.r.go, aes(x = factor(yrseq), y = resp, group = sampid)) +
    geom_point() + geom_line() + 
    xlab("Year sequenced") + ylab(yLab) + mytheme
  
  result <- list(mod = mod, p = p)
  return(result)
  
}

test_yrseq_qtaq <- function(df.r, yLab){
  
  df.r %>%
    filter(yrharv == 2014) %>%
    filter(taqType != "GoTaq") -> df.r.q
  
  mod <- lmer(resp ~ yrseq + (1|sampid), data = df.r.q)
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(df.r.q, aes(x = factor(yrseq), y = resp, group = sampid)) +
    geom_point() + geom_line() + 
    xlab("Year sequenced") + ylab(yLab) + mytheme
  
  result <- list(mod = mod, p = p)
  return(result)
  
}

test_taqdiffs_composition <- function(mat.otu.pres, envVars, logTrt){
  
  #mat.otu.pres = tab.17
  #envVars = samp.17
  #logTrt
  
  # test
  cap.env <- capscale(mat.otu.pres ~ taqType + yrharv,
                      data = envVars, distance = 'jaccard', binary = TRUE)
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  anova.df
  
  ######
  # plot
  
  # do the unconstrained ordination
  mod.obj <- capscale(mat.otu.pres ~ 1, data = NULL, distance = "jaccard", binary = TRUE)
  
  # do the environmental fit
  envVars %>%
    mutate(grp = paste(taqType, yrharv, sep = "_")) -> envVars
  envfit.obj <- envfit(mod.obj ~ grp, envVars)
  
  # make a point dataframe
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  scrs <- cbind(scrs, tidySeq_sampleName = row.names(scrs))
  scrs %>%
    left_join(envVars) %>%
    left_join(logTrt) %>%
    mutate(dep.harv = paste(yrdeploy, yrharv, sep = "_")) -> scrs
  
  # make segment dataframe
  scrs %>%
    select(MDS1, sampid, taqType) %>%
    spread(key = "taqType", value = MDS1) %>%
    rename('g.MDS1'='GoTaq',
           'q.MDS1'='Q5') %>%
    arrange(sampid) -> seg.1
  scrs %>%
    select(MDS2, sampid, taqType) %>%
    spread(key = "taqType", value = MDS2) %>%
    rename('g.MDS2'='GoTaq',
           'q.MDS2'='Q5') %>%
    arrange(sampid) %>%
    left_join(seg.1) %>%
    separate(sampid, into = c("plot","logtrt","logLoc","yrharv"), sep = "_", remove = F) %>%
    left_join(logTrt) %>%
    mutate(dep.harv = paste(yrdeploy, yrharv, sep = "_")) -> seg
  
  params <- plottingParams_rotplotNobags()
  scrs$dep.harv1 <- factor(scrs$dep.harv, 
                          levels = as.character(params$yearSubset$levels2.new),
                          labels = params$yearSubset$labels.new)
  seg$dep.harv1 <- factor(seg$dep.harv,
                         levels = as.character(params$yearSubset$levels2.new),
                         labels = params$yearSubset$labels.new)
  
  # colors
  light.red <- brewer.pal(n=8, name = "RdBu")[2] # 2011-12 (1yr)
  dark.red <- brewer.pal(n=8, name = "RdBu")[1] # 2011-2014 (3yrs)
  
  light.blue <- brewer.pal(n=8, name = "RdBu")[7] # 2009-2012 (3yrs)
  dark.blue <- brewer.pal(n=8, name = "RdBu")[8] # 2009-2014 (5yrs)
  
  light.purple <- brewer.pal(n=4, name = "Purples")[3] # overlap harvested in 2012
  dark.purple <- brewer.pal(n=4, name = "Purples")[4] # overlap harvested in 2014
  
  
  # plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(scrs, aes(x = MDS1, y = MDS2, color = dep.harv1, shape = taqType)) +
    geom_point() +
    geom_segment(data = seg, 
                 aes(x = g.MDS1, y = g.MDS2, 
                     xend = q.MDS1, yend = q.MDS2, 
                     color = dep.harv1), inherit.aes = F) +
    coord_fixed(ratio = 1) +
    xlab("PCoA 1") + ylab("PCoA 2") + mytheme +
    scale_color_manual(name = "Deployment subset", 
                       values = c(light.red, dark.red, light.blue, dark.blue)) +
    scale_shape_discrete(name = "Taq type") 
  p
  
  result <- list(anova.df = anova.df, p = p)
  return(result)
  
}

test_yrseq_composition <- function(mat.otu.pres, envVars, logTrt, yrharv){
  
  # test
  cap.env <- capscale(mat.otu.pres ~ yrseq, data = envVars, distance = 'jaccard', binary = TRUE)
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  anova.df
  
  ######
  # plot
  
  # do the unconstrained ordination
  mod.obj <- capscale(mat.otu.pres ~ 1, data = NULL, distance = "jaccard", binary = TRUE)
  
  # do the environmental fit
  envfit.obj <- envfit(mod.obj ~ yrseq, envVars)
  
  # make a point dataframe
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  scrs <- cbind(scrs, tidySeq_sampleName = row.names(scrs))
  scrs %>%
    left_join(envVars) %>%
    left_join(logTrt) -> scrs
  
  # make segment dataframe
  if(yrharv == 2012){
    scrs %>%
      select(MDS1, sampid, yrseq) %>%
      spread(key = "yrseq", value = MDS1) %>%
      rename('yf.MDS1'='2012',
             'yl.MDS1'='2017') %>%
      arrange(sampid) -> seg.1
    scrs %>%
      select(MDS2, sampid, yrseq) %>%
      spread(key = "yrseq", value = MDS2) %>%
      rename('yf.MDS2'='2012',
             'yl.MDS2'='2017') %>%
      arrange(sampid) %>%
      left_join(seg.1) %>%
      separate(sampid, into = c("plot","logtrt","logLoc","yrharv"), sep = "_", remove = F) %>%
      left_join(logTrt) -> seg
  }
  if(yrharv == 2014){
    scrs %>%
      select(MDS1, sampid, yrseq) %>%
      spread(key = "yrseq", value = MDS1) %>%
      rename('yf.MDS1'='2014',
             'yl.MDS1'='2017') %>%
      arrange(sampid) -> seg.1
    scrs %>%
      select(MDS2, sampid, yrseq) %>%
      spread(key = "yrseq", value = MDS2) %>%
      rename('yf.MDS2'='2014',
             'yl.MDS2'='2017') %>%
      arrange(sampid) %>%
      left_join(seg.1) %>%
      separate(sampid, into = c("plot","logtrt","logLoc","yrharv"), sep = "_", remove = F) %>%
      left_join(logTrt) -> seg
  }
  
  
  # plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(scrs) +
    geom_point(mapping = aes(x = MDS1, y = MDS2, shape = factor(yrseq), color = species4, alpha = yrdeploy)) +
    geom_segment(data = seg, mapping = aes(x = yf.MDS1, y = yf.MDS2, 
                                           xend = yl.MDS1, yend = yl.MDS2, color = species4, alpha = yrdeploy)) +
    coord_fixed(ratio = 1) +  ## need aspect ratio of 1!
    xlab("PCoA 1") + ylab("PCoA 2") + mytheme +
    scale_alpha_manual(name = "Year deployed", values = c(1, .5)) +
    scale_shape_discrete(name = "Year sequenced") +
    scale_color_discrete(name = "Wood species")
  p
  
  result <- list(anova.df = anova.df, p = p)
  return(result)
  
  
}


# richness ~ study factors #-----------------------------------------------#

richnessTest <- function(data, datasetType){

  # calculate sample richness and run lm
  summ.coef.list <- list()
  summ.fit.list <- list()
  for(i in 1:length(data)){
    
    # calc richness
    x <- data[[i]]
    richness <- apply(x$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    x$meta %>%
      left_join(richness.df) %>%
      filter(logLoc != "mush") -> rich.df
    
    # run lm
    if(datasetType == "overlap"){
      mod <- lm(rich ~ logLoc + topo + species4 + yrdeploy, data = rich.df)
    }
    if(datasetType == "deploy"){
      mod <- lm(rich ~ logLoc + topo + species4, data = rich.df)
    }
    
    #save coefs and fit stats
    summ.coef.list[[i]] <- data.frame(term =  row.names(summary(mod)$coefficients), summary(mod)$coefficients)
    summ.fit.list[[i]] <- data.frame(fstat = summary(mod)$fstatistic[['value']],
                                     numdf = summary(mod)$fstatistic[['numdf']],
                                     dendf = summary(mod)$fstatistic[['dendf']],
                                     r.squared = summary(mod)$r.squared)
  }
  names(summ.coef.list) <- names(data)
  names(summ.fit.list) <- names(data)
  
  coef.df <- list_to_df(summ.coef.list)
  coef.df %>%
    separate(source, into = c("yearGene","subset")) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    rename('p.value' = 'Pr...t..') %>%
    mutate(signif = ifelse(p.value < .05, TRUE, FALSE)) -> coef.df
  fit.df <- list_to_df(summ.fit.list)
  fit.df %>%
    rename('subset'='source') -> fit.df
  result <- list(coef.df = coef.df, fit.df = fit.df)
  
  return(result)
}

otuprep_richnessTest <- function(data.otuPreplist, datasetType){
  
  #calc and extract lm richness coefs
  rich.list <- lapply(data.otuPreplist, function(x){
    tmp <- richnessTest(data = x, datasetType = datasetType)
    return(tmp$coef.df)
    })
  rich.df <- list_to_df(rich.list)
  rich.df %>%
    rename('otuprep'='source') -> rich.df
  
  # calc and extract lm fit statistics
  rich.fit.list <- lapply(data.otuPreplist, function(x){
    tmp <- richnessTest(data = x, datasetType = datasetType)
    return(tmp$fit.df)
  })
  rich.fit.df <- list_to_df(rich.fit.list)
  rich.fit.df %>%
    rename('otuprep'='source') -> rich.fit.df
  
  result <- list(rich.df = rich.df, rich.fit.df = rich.fit.df)
  
  return(result)
}


# composition ~ study factors #-----------------------------------------------#

calcCentroidDist <- function(cap.env, FACTORS){
  
  centroids <- data.frame(level = row.names(summary(cap.env)$centroids), summary(cap.env)$centroids)
  centroids
  all.caps <- colnames(centroids)[grepl("CAP", colnames(centroids))]
  #FACTORS <- c("logLoc","topo","species4","yrdeploy")
  dist.list <- list()
  for(i in 1:length(FACTORS)){
    
    # all constrained axes
    distmat.allcaps <- as.matrix(dist(centroids[grepl(FACTORS[i], centroids$level), all.caps]))
    dist.allcaps <- distmat.allcaps[upper.tri(distmat.allcaps)]
    
    # first 2 constrained axes
    distmat <- as.matrix(dist(centroids[grepl(FACTORS[i], centroids$level), c("CAP1","CAP2")]))
    dist <- distmat[upper.tri(distmat)]
    
    # if more than 2 levels, calc mean distance
    nlevels <- sum(grepl(FACTORS[i], centroids$level))
    dist.list[[i]] <- data.frame(mean.distall = mean(dist.allcaps),
                                 se.distall = sd(dist.allcaps)/sqrt(nlevels),
                                 mean.dist = mean(dist),
                                 se.dist = sd(dist)/sqrt(nlevels),
                                 n = nlevels)
  }
  names(dist.list) <- FACTORS
  dist.df <- list_to_df(dist.list)
  dist.df %>%
    rename('terms'='source') -> dist.df
  
  return(dist.df)
  
}

dbrdaTest_presAbs <- function(data, datasetType){
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  
  # remove logLoc == "mush"
  condition <- sum(grepl("mush", envVars$logLoc))
  if(condition !=0){
    envVars %>%
      filter(logLoc == "mush") -> remove
    mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
    envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
  }
  
  #full model anova
  if(datasetType == "overlap"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + species4 + yrdeploy,
                        data = envVars, distance = 'jaccard', binary = TRUE)
  }
  if(datasetType == "deploy"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + species4,
                        data = envVars, distance = 'jaccard', binary = TRUE)
  }
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  
  #calculate centroid distances - if factor has more than 2 levels, take the mean distance
  if(datasetType == "overlap"){
    select.factors <- c("logLoc","topo","species4","yrdeploy")
  }
  if(datasetType == "deploy"){
    select.factors <- c("logLoc","topo","species4")
  }
  dist.df <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
  
  #do variance partitioning
  if(datasetType == "overlap"){
    endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                        ~ species4,
                        ~ yrdeploy,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:4,]
    terms <- c("species4","yrdeploy","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ logLoc 
                          + Condition(envVars$yrdeploy) + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ yrdeploy 
                          + Condition(envVars$logLoc) + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[2,"pval"] <- an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ topo 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[3,"pval"] <- an3$`Pr(>F)`[1]
    
    # significance of partition X4
    an4 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ species4 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$topo), data = envVars))
    var.fract[4,"pval"] <- an4$`Pr(>F)`[1]
  }
  if(datasetType == "deploy"){
    endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                        ~ species4,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:3,]
    terms <- c("species4","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ logLoc
                          + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ topo
                          + Condition(envVars$logLoc)
                          + Condition(envVars$species4), data = envVars))
    var.fract[2,"pval"]<-an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ species4
                          + Condition(envVars$logLoc)
                          + Condition(envVars$topo), data = envVars))
    var.fract[3,"pval"]<-an3$`Pr(>F)`[1]
  }
  
  result <- list(anova.df = anova.df, var.fract = var.fract, dist.df = dist.df)
  return(result)
}

dbrdaTest_clr <- function(data, datasetType){
  
  #data = x
  #datasetType = datasetType
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  
  # remove logLoc == "mush"
  condition <- sum(grepl("mush", envVars$logLoc))
  if(condition !=0){
    envVars %>%
      filter(logLoc == "mush") -> remove
    mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
    envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
  }
  
  #full model anova
  datasetType
  if(datasetType == "overlap"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + species4 + yrdeploy,
                        data = envVars, distance = 'euclidean')
  }
  if(datasetType == "deploy"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + species4,
                        data = envVars, distance = 'euclidean')
  }
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  
  #calculate centroid distances - if factor has more than 2 levels, take the mean distance
  if(datasetType == "overlap"){
    select.factors <- c("logLoc","topo","species4","yrdeploy")
  }
  if(datasetType == "deploy"){
    select.factors <- c("logLoc","topo","species4")
  }
  dist.df <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
  
  #do variance partitioning
  if(datasetType == "overlap"){
    endo.var <- varpart(vegdist(mat.otu.pres, method = 'euclidean'),
                        ~ species4,
                        ~ yrdeploy,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:4,]
    terms <- c("species4","yrdeploy","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ logLoc 
                          + Condition(envVars$yrdeploy) + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ yrdeploy 
                          + Condition(envVars$logLoc) + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[2,"pval"] <- an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ topo 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[3,"pval"] <- an3$`Pr(>F)`[1]
    
    # significance of partition X4
    an4 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ species4 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$topo), data = envVars))
    var.fract[4,"pval"] <- an4$`Pr(>F)`[1]
  }
  if(datasetType == "deploy"){
    endo.var <- varpart(vegdist(mat.otu.pres, method = 'euclidean'),
                        ~ species4,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:3,]
    terms <- c("species4","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ logLoc
                          + Condition(envVars$topo) 
                          + Condition(envVars$species4), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ topo
                          + Condition(envVars$logLoc)
                          + Condition(envVars$species4), data = envVars))
    var.fract[2,"pval"]<-an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ species4
                          + Condition(envVars$logLoc)
                          + Condition(envVars$topo), data = envVars))
    var.fract[3,"pval"]<-an3$`Pr(>F)`[1]
  }
  
  result <- list(anova.df = anova.df, var.fract = var.fract, dist.df = dist.df)
  return(result)
}

otuprep_dbrdaTest <- function(data.otuPreplist, datasetType){
  
  an.list <- list()
  varpart.list <- list()
  dist.list <- list()
  
  for(i in 1:length(data.otuPreplist)){
    
    data <- data.otuPreplist[[i]]
    print(names(data.otuPreplist)[i])
    if(i == 1){ # just for CLR data
      tmp.list <- lapply(data, function(x){
        # removes "mush" level first
        tmp <- dbrdaTest_clr(data = x, datasetType = datasetType)
        # rename the column name so that we can merge everything into a dataframe downstream
        tmp$anova.df %>%
          dplyr::rename('SumOfSqs'='Variance') -> tmp$anova.df 
        return(tmp)
      })
    }else{
      # calculate full model anova, centroid distances, and varpart
      tmp.list <- lapply(data, function(x){
        # removes "mush" level first
        tmp <- dbrdaTest_presAbs(data = x, datasetType = datasetType)
        return(tmp)
      })
    }
    
    tmp.df <- list_to_df(lapply(tmp.list, function(x){x$anova.df}))
    tmp.df %>%
      rename('subset' = 'source') -> an.list[[i]]
    tmp.df <- list_to_df(lapply(tmp.list, function(x){x$var.fract}))
    tmp.df %>%
      rename('subset' = 'source') -> varpart.list[[i]]
    tmp.df <- list_to_df(lapply(tmp.list, function(x){x$dist.df}))
    tmp.df %>%
      rename('subset' = 'source') -> dist.list[[i]]
    
  }
  names(an.list) <- names(data.otuPreplist)
  names(varpart.list) <- names(data.otuPreplist)
  names(dist.list) <- names(data.otuPreplist)
  
  #anova
  an.df <- list_to_df(an.list)
  an.df %>%
    rename('otuprep'='source') -> an.df
  
  #varpart
  varpart.df <- list_to_df(varpart.list)
  varpart.df %>%
    separate(subset, into = c("yrharv","geneDrop"), sep = 4, remove = F) %>%
    separate(geneDrop, into = c("gene","group","drop"), fill = "right") %>%
    mutate(group = replace(group, group == "o", NA)) %>%
    select(-drop) %>%
    rename('otuprep'='source') -> varpart.df
  
  #centroid distances
  dist.df <- list_to_df(dist.list)
  dist.df %>%
    separate(subset, into = c("yrharv","geneDrop"), sep = 4, remove = F) %>%
    separate(geneDrop, into = c("gene","group","drop"), fill = "right") %>%
    mutate(group = replace(group, group == "o", NA)) %>%
    select(-drop) %>%
    rename('otuprep'='source') -> dist.df
  
  result <- list(an.df = an.df, varpart.df = varpart.df, dist.df = dist.df)
  return(result)
}

dbrdaTest_ag_presAbs <- function(data, datasetType){
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  
  # remove logLoc == "mush"
  envVars %>%
    filter(logLoc == "mush") -> remove
  mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
  envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
  
  #full model anova
  if(datasetType == "overlap"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + angio.gymno + yrdeploy,
                        data = envVars, distance = 'jaccard', binary = TRUE)
  }
  if(datasetType == "deploy"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + angio.gymno,
                        data = envVars, distance = 'jaccard', binary = TRUE)
  }
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  
  #calculate centroid distances - if factor has more than 2 levels, take the mean distance
  if(datasetType == "overlap"){
    select.factors <- c("logLoc","topo","angio.gymno","yrdeploy")
  }
  if(datasetType == "deploy"){
    select.factors <- c("logLoc","topo","angio.gymno")
  }
  dist.df <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
  
  #do variance partitioning
  if(datasetType == "overlap"){
    endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                        ~ angio.gymno,
                        ~ yrdeploy,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:4,]
    terms <- c("angio.gymno","yrdeploy","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ logLoc 
                          + Condition(envVars$yrdeploy) + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ yrdeploy 
                          + Condition(envVars$logLoc) + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[2,"pval"] <- an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ topo 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[3,"pval"] <- an3$`Pr(>F)`[1]
    
    # significance of partition X4
    an4 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ angio.gymno 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$topo), data = envVars))
    var.fract[4,"pval"] <- an4$`Pr(>F)`[1]
  }
  if(datasetType == "deploy"){
    endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                        ~ angio.gymno,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:3,]
    terms <- c("angio.gymno","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ logLoc
                          + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ topo
                          + Condition(envVars$logLoc)
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[2,"pval"]<-an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE)	~ angio.gymno
                          + Condition(envVars$logLoc)
                          + Condition(envVars$topo), data = envVars))
    var.fract[3,"pval"]<-an3$`Pr(>F)`[1]
  }
  result <- list(anova.df = anova.df, var.fract = var.fract, dist.df = dist.df)
  return(result)
}

dbrdaTest_ag_clr <- function(data, datasetType){
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  
  # remove logLoc == "mush"
  envVars %>%
    filter(logLoc == "mush") -> remove
  mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
  envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
  
  #full model anova
  if(datasetType == "overlap"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + angio.gymno + yrdeploy,
                        data = envVars, distance = 'euclidean')
  }
  if(datasetType == "deploy"){
    cap.env <- capscale(mat.otu.pres ~ logLoc + topo + angio.gymno,
                        data = envVars, distance = 'euclidean')
  }
  anova.df <- anova(cap.env, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(cap.env)$constr.chi
  anova.df$unconst.chi <- summary(cap.env)$unconst.chi
  anova.df$tot.chi <- summary(cap.env)$tot.chi
  
  #calculate centroid distances - if factor has more than 2 levels, take the mean distance
  if(datasetType == "overlap"){
    select.factors <- c("logLoc","topo","angio.gymno","yrdeploy")
  }
  if(datasetType == "deploy"){
    select.factors <- c("logLoc","topo","angio.gymno")
  }
  dist.df <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
  
  #do variance partitioning
  if(datasetType == "overlap"){
    endo.var <- varpart(vegdist(mat.otu.pres, method = 'euclidean'),
                        ~ angio.gymno,
                        ~ yrdeploy,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:4,]
    terms <- c("angio.gymno","yrdeploy","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ logLoc 
                          + Condition(envVars$yrdeploy) + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ yrdeploy 
                          + Condition(envVars$logLoc) + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[2,"pval"] <- an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ topo 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[3,"pval"] <- an3$`Pr(>F)`[1]
    
    # significance of partition X4
    an4 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ angio.gymno 
                          + Condition(envVars$logLoc) + Condition(envVars$yrdeploy) 
                          + Condition(envVars$topo), data = envVars))
    var.fract[4,"pval"] <- an4$`Pr(>F)`[1]
  }
  if(datasetType == "deploy"){
    endo.var <- varpart(vegdist(mat.otu.pres, method = 'euclidean'),
                        ~ angio.gymno,
                        ~ topo,
                        ~ logLoc, data = envVars)
    var.fract <- endo.var$part$indfract[1:3,]
    terms <- c("angio.gymno","topo","logLoc")
    var.fract$terms <- terms
    
    # significance of partition X1
    an1 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ logLoc
                          + Condition(envVars$topo) 
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[1,"pval"] <- an1$`Pr(>F)`[1]
    
    # significance of partition X2
    an2 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ topo
                          + Condition(envVars$logLoc)
                          + Condition(envVars$angio.gymno), data = envVars))
    var.fract[2,"pval"]<-an2$`Pr(>F)`[1]
    
    # significance of partition X3
    an3 <- anova(capscale(vegdist(mat.otu.pres, method = 'euclidean')	~ angio.gymno
                          + Condition(envVars$logLoc)
                          + Condition(envVars$topo), data = envVars))
    var.fract[3,"pval"]<-an3$`Pr(>F)`[1]
  }
  result <- list(anova.df = anova.df, var.fract = var.fract, dist.df = dist.df)
  return(result)
}

# richness ~ wood traits #-----------------------------------------------#

richnessTest_woodtraits_rf <- function(data, woodTraits){
  
  require(randomForest)
  
  # calculate sample richness and run rf will all wood traits
  imp.list <- list()
  r2.list <- list()
  for(i in 1:length(data)){
    
    # calc richness
    x <- data[[i]]
    richness <- apply(x$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    x$meta %>%
      left_join(richness.df) %>%
      select(species4, rich) %>%
      left_join(woodTraits) %>%
      select(-c(species4)) -> rich.df
   
    # do random forest and extract variable importance and pseudo-r2
    rf <- randomForest(rich ~., data = rich.df, importance = TRUE, na.action = na.exclude)
    # variable importance
    tmp <- data.frame(term = row.names(rf$importance), rf$importance, rf$importanceSD)
    tmp$X.IncMSE.s <- scale(tmp$X.IncMSE, center = F)
    tmp$rf.importanceSD.s <- tmp$rf.importanceSD / attr(tmp$X.IncMSE.s, 'scaled:scale')
    tmp$X.IncMSE.s <- as.numeric(tmp$X.IncMSE.s)
    imp.list[[i]] <- tmp
    # pseudo-r2
    r2.list[[i]] <- data.frame(r2.mean = mean(rf$rsq), r2.sd = sd(rf$rsq))
  }
  names(imp.list) <- names(data)
  names(r2.list) <- names(data)
  
  imp.df <- list_to_df(imp.list)
  imp.df %>%
    separate(source, into = c("yearGene","subset")) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> imp.df
  
  r2.df <- list_to_df(r2.list)
  r2.df
  
  result <- list(imp.df = imp.df, r2.df = r2.df)
  
  return(result)
}

otuprep_richnessTest_woodtraits_rf <- function(data.otuPreplist, woodTraits){
  
  #calc rf richness variable importance and fit stats
  rich.list <- lapply(data.otuPreplist, function(x){
    tmp <- richnessTest_woodtraits_rf(data = x, woodTraits)
    return(tmp)
  })
  
  #extract coefs
  rich.df <- lapply(rich.list, function(x){x$imp.df})
  rich.df <- list_to_df(rich.df)
  rich.df %>%
    rename('otuprep'='source') -> rich.df
  
  #extract fit stats
  rich.fit.df <- lapply(rich.list, function(x){x$r2.df})
  rich.fit.df <- list_to_df(rich.fit.df)
  rich.fit.df %>%
    rename('otuprep'='source') -> rich.fit.df
  
  result <- list(rich.df = rich.df, rich.fit.df = rich.fit.df)
  
  return(result)
}

richnessTest_woodtraits_step <- function(data, woodTraits){
  
  # calculate sample richness, set up regression, do model selection 
  an.list <- list()
  for(i in 1:length(data)){
    
    # calc richness
    x <- data[[i]]
    richness <- apply(x$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    x$meta %>%
      left_join(richness.df) %>%
      select(species4, rich) %>%
      left_join(woodTraits) %>%
      select(-c(species4)) -> rich.df
    
    # set up model
    mod <- lm(rich ~ ., data = rich.df)
    # do stepwise model selection
    mod.s <- step(mod)
    # save final model
    an <- anova(mod.s)
    summ <- summary(mod.s)
    an.list[[i]] <- data.frame(term = row.names(an), an, adj.r2 = summ$adj.r.squared)
    
  }
  names(an.list) <- names(data)
  an.df <- list_to_df(an.list)
  an.df %>%
    separate(source, into = c("yearGene","subset")) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> an.df
  
  return(an.df)
  
}


# composition ~ wood traits #-----------------------------------------------#

ordiR2step_woodtraits <- function(data, woodTraits){
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
  
  # trait data for species4 == PLOC is missing (only in d09), so exclude it
  envVars %>%
    filter(species4 != "PLOC") -> envVars
  mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
  
  # make the fully-constrained and unconstrained ordinations
  envVars %>%
    left_join(woodTraits) %>%
    select(trait.names) -> envVars.to
  cap <- capscale(mat.otu.pres ~ 1, data = envVars.to, distance = 'jaccard', binary = TRUE)
  cap.env <-  capscale(mat.otu.pres ~ ., data = envVars.to, distance = 'jaccard', binary = TRUE)
  
  # do model selection
  os <- ordiR2step(cap, cap.env)
  anova.df <- anova(os, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(os)$constr.chi
  anova.df$unconst.chi <- summary(os)$unconst.chi
  anova.df$tot.chi <- summary(os)$tot.chi
  
  result <- list(ordistep.sum = os$anova, anova = anova.df)
  return(result)

}

ordiR2step_woodtraits_clr <- function(data, woodTraits){
  
  mat.otu.pres <- data[["otu"]]
  envVars <- data[["meta"]]
  trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
  
  # trait data for species4 == PLOC is missing (only in d09), so exclude it
  envVars %>%
    filter(species4 != "PLOC") -> envVars
  mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
  
  # make the fully-constrained and unconstrained ordinations
  envVars %>%
    left_join(woodTraits) %>%
    select(trait.names) -> envVars.to
  cap <- capscale(mat.otu.pres ~ 1, data = envVars.to, distance = 'euclidean')
  cap.env <-  capscale(mat.otu.pres ~ ., data = envVars.to, distance = 'euclidean')
  
  # do model selection
  os <- ordiR2step(cap, cap.env)
  anova.df <- anova(os, by = "terms")
  anova.df <- data.frame(term = row.names(anova.df), anova.df)
  anova.df$constr.chi <- summary(os)$constr.chi
  anova.df$unconst.chi <- summary(os)$unconst.chi
  anova.df$tot.chi <- summary(os)$tot.chi
  
  result <- list(ordistep.sum = os$anova, anova = anova.df)
  return(result)
  
}

otuprep_ordiR2step_woodtraits <- function(data.otuPreplist, woodTraits){
  
  anova.list <- list()
  for(i in 1:length(data.otuPreplist)){
    
    data <- data.otuPreplist[[i]]
    print(names(data.otuPreplist)[i])
    
    # do ordistep
    tmp.list <- lapply(data, function(x){
      tmp <- ordiR2step_woodtraits(data = x, woodTraits = woodTraits)
      return(tmp)
    })
    tmp.list
    tmp.df <- list_to_df(lapply(tmp.list, function(x){x$anova}))
    tmp.df %>%
      rename('subset' = 'source') -> anova.list[[i]]
    
  }
  names(anova.list) <- names(data.otuPreplist)
  
  #os
  anova.list
  an.df <- list_to_df(anova.list)
  an.df %>%
    rename('otuprep'='source') -> an.df
  
  return(an.df)
}

extract_rhs.terms <- function(rd, table.name, otuprep.name){

  rd %>%
    mutate(tmp = paste0(yrharv, gene)) %>%
    mutate(table = paste(tmp, yrdeploy, sep = "_")) %>%
    filter(table == table.name) %>%
    filter(otuprep == otuprep.name) %>%
    filter(term != "Residual") -> curr
  terms <- c(as.character(curr$term))
  string <- paste(terms, collapse = " + ")
  
  return(terms)
}

dbrdaTest_woodTraits_finalmod_v <- function(rd.all, woodTraits, data.deploy){
  
  # exclude datasets for sapros and basidios
  condition <- !grepl("sapro", names(data.deploy)) & !grepl("basidio", names(data.deploy))
  data.deploy <- data.deploy[condition]
  
  table.names <- names(data.deploy)
  trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
  cca.biplot <- list()
  v.list <- list()
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
    
    #vector info
    cca.biplot[[i]] <- data.frame(terms = row.names(cap.env$CCA$biplot), cap.env$CCA$biplot)
    
    # top OTUs
    v.list[[i]] <- data.frame(OTUId = row.names(cap.env$CCA$v), cap.env$CCA$v)
    
  }
  names(cca.biplot) <- table.names
  names(v.list) <- table.names
  
  return(v.list)
  
}

extract_keyOTUs <- function(v.list, cap.choice, table.choice, tab.taxa){
  
  tax.name <- strsplit(table.choice, "_")[[1]][1]
  tab.taxa.df <- tab.taxa[[tax.name]]
  v.df <- v.list[[table.choice]]
  v.df %>%
    select(OTUId, cap.choice) %>%
    rename('OTUid'='OTUId') %>%
    left_join(tab.taxa.df) %>%
    arrange_(cap.choice) -> tmp
  
  return(tmp)
}


# composition ~ study factors and key wood traits #-----------------------------------------------#

dbrdaTest_woodTraits_presAbs_deploy <- function(rd, woodTraits, data.deploy){
  
  table.names <- names(data.deploy)
  trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
  
  an.list <- list()
  dist.list <- list()
  var.fract.list <- list()
  for(i in 1:length(data.deploy)){
    
    terms <- extract_rhs.terms(rd.all = rd.all, table.name = table.names[i])
    mat.otu.pres <- data.deploy[[i]][["otu"]]
    envVars <- data.deploy[[i]][["meta"]]
    
    # trait data for species4 == PLOC is missing (only in d09), so exclude it
    envVars %>%
      filter(species4 != "PLOC") -> envVars
    mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
    
    # remove logLoc == "mush"
    envVars %>%
      filter(logLoc == "mush") -> remove
    mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
    envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
    
    # set up the model
    envVars %>%
      left_join(woodTraits) %>%
      select(terms) -> envVars.to
    envVars %>%
      select(logLoc, topo) -> add.envVars
    envVars.tmp <- cbind(envVars.to, add.envVars)
    cap.env <-  capscale(mat.otu.pres ~ ., 
                         data = envVars.tmp, 
                         distance = 'jaccard', binary = TRUE)
    
    # anova
    anova.df <- anova(cap.env, by = "terms")
    anova.df <- data.frame(term = row.names(anova.df), anova.df)
    anova.df$constr.chi <- summary(cap.env)$constr.chi
    anova.df$unconst.chi <- summary(cap.env)$unconst.chi
    anova.df$tot.chi <- summary(cap.env)$tot.chi
    an.list[[i]] <- anova.df
    
    #calculate centroid distances - if factor has more than 2 levels, take the mean distance
    select.factors <- c("logLoc","topo")
    dist.list[[i]] <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
    
    #do variance partitioning
    endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                        ~ as.matrix(envVars.to),
                        ~ envVars.tmp$topo,
                        ~ envVars.tmp$logLoc)
    endo.var$part$indfract[1:3,]
    
    var.fract <- endo.var$part$indfract[1:3,]
    terms <- c("woodtraits","topo","logLoc")
    var.fract$terms <- terms
    var.fract.list[[i]] <- var.fract
    
  }
  names(an.list) <- names(dist.list) <- names(var.fract.list) <- names(data.deploy)
  
  an.df <- list_to_df(an.list)
  dist.df <- list_to_df(dist.list)
  var.fract.df <- list_to_df(var.fract.list)
  
  result <- list(anova.df = anova.df, var.fract = var.fract, dist.df = dist.df)
  return(result)
}

otuprep_dbrdaTest_woodTraits <- function(rd, woodTraits, data.otuPreplist){
  
  an.list2 <- list()
  varpart.list2 <- list()
  dist.list2 <- list()
  for(i in 1:length(data.otuPreplist)){
    
    data <- data.otuPreplist[[i]]
    # exclude datasets for sapros and basidios
    condition <- !grepl("sapro", names(data)) & !grepl("basidio", names(data))
    data <- data[condition]
    
    # calculate full model anova, centroid distances, and varpart
    # for each data element...
    table.names <- names(data)
    trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
    print(names(data.otuPreplist)[i])
    otuprep.name <- paste0(names(data.otuPreplist)[i], "pa")
    
    an.list <- list()
    dist.list <- list()
    var.fract.list <- list()
    for(k in 1:length(data)){
      terms <- extract_rhs.terms(rd = rd, 
                                 table.name = table.names[k], 
                                 otuprep.name = otuprep.name)
      mat.otu.pres <- data[[k]][["otu"]]
      envVars <- data[[k]][["meta"]]
      
      # trait data for species4 == PLOC is missing (only in d09), so exclude it
      envVars %>%
        filter(species4 != "PLOC") -> envVars
      mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
      
      # remove logLoc == "mush"
      envVars %>%
        filter(logLoc == "mush") -> remove
      mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
      envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
      
      # set up the model
      envVars %>%
        left_join(woodTraits) %>%
        select(terms) -> envVars.to
      envVars %>%
        select(logLoc, topo) -> add.envVars
      envVars.tmp <- cbind(envVars.to, add.envVars)
      cap.env <-  capscale(mat.otu.pres ~ ., 
                           data = envVars.tmp, 
                           distance = 'jaccard', binary = TRUE)
      
      # anova
      anova.df <- anova(cap.env, by = "terms")
      anova.df <- data.frame(term = row.names(anova.df), anova.df)
      anova.df$constr.chi <- summary(cap.env)$constr.chi
      anova.df$unconst.chi <- summary(cap.env)$unconst.chi
      anova.df$tot.chi <- summary(cap.env)$tot.chi
      an.list[[k]] <- anova.df
      
      #calculate centroid distances - if factor has more than 2 levels, take the mean distance
      select.factors <- c("logLoc","topo")
      dist.list[[k]] <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
      
      #do variance partitioning
      endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                          ~ as.matrix(envVars.to),
                          ~ envVars.tmp$topo,
                          ~ envVars.tmp$logLoc)
      endo.var$part$indfract[1:3,]
      
      var.fract <- endo.var$part$indfract[1:3,]
      terms <- c("woodtraits","topo","logLoc")
      var.fract$terms <- terms
      var.fract.list[[k]] <- var.fract
      
    }
    names(an.list) <- names(dist.list) <- names(var.fract.list) <- table.names
    
    tmp <- list_to_df(an.list)
    tmp %>%
      rename('subset'='source') -> an.list2[[i]]
    tmp <- list_to_df(dist.list)
    tmp %>%
      rename('subset'='source') -> dist.list2[[i]]
    tmp <- list_to_df(var.fract.list)
    tmp %>%
      rename('subset'='source') -> varpart.list2[[i]]
  }
  names(an.list2) <- names(dist.list2) <- names(varpart.list2) <- names(data.otuPreplist)
  
  tmp <- list_to_df(an.list2)
  tmp %>%
    rename('otuprep'='source') -> an.df
  tmp <- list_to_df(varpart.list2)
  tmp %>%
    rename('otuprep'='source') -> varpart.df
  tmp <- list_to_df(dist.list2)
  tmp %>%
    rename('otuprep'='source') -> dist.df
  
  result <- list(an.df = an.df, varpart.df = varpart.df, dist.df = dist.df)
  return(result)
}

otuprep_dbrdaTest_woodTraits_clr <- function(rd, woodTraits, data.otuPreplist){
  
  rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d_clr.RData")
  woodTraits <- readRDS(file = "data_Rsynth/rotplot_intermediates/woodTraits.RData")
  deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))
  data.otuPreplist <- deploy.otuPreplist
  
  
  an.list2 <- list()
  varpart.list2 <- list()
  dist.list2 <- list()

  for(i in 1:length(data.otuPreplist)){
    
    data <- data.otuPreplist[[i]]
    # exclude datasets for sapros and basidios
    condition <- !grepl("sapro", names(data)) & !grepl("basidio", names(data))
    data <- data[condition]
    
    # calculate full model anova, centroid distances, and varpart
    # for each data element...
    table.names <- names(data)
    trait.names <- colnames(woodTraits)[colnames(woodTraits) != "species4"]
    print(names(data.otuPreplist)[i])
    
    an.list <- list()
    dist.list <- list()
    var.fract.list <- list()
    
    if(i == 1){
      otuprep.name <- paste0(names(data.otuPreplist)[i])
      for(k in 1:length(data)){
        terms <- extract_rhs.terms(rd = rd, 
                                   table.name = table.names[k], 
                                   otuprep.name = otuprep.name)
        mat.otu.pres <- data[[k]][["otu"]]
        envVars <- data[[k]][["meta"]]
        
        # trait data for species4 == PLOC is missing (only in d09), so exclude it
        envVars %>%
          filter(species4 != "PLOC") -> envVars
        mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
        
        # remove logLoc == "mush"
        envVars %>%
          filter(logLoc == "mush") -> remove
        mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
        envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
        
        # set up the model
        envVars %>%
          left_join(woodTraits) %>%
          select(terms) -> envVars.to
        envVars.to.s <- scale(envVars.to)
        envVars %>%
          select(logLoc, topo) -> add.envVars
        envVars.tmp <- cbind(envVars.to.s, add.envVars)
        cap.env <-  capscale(mat.otu.pres ~ ., 
                             data = envVars.tmp, 
                             method = 'euclidean')
        # anova
        anova.df <- anova(cap.env, by = "terms")
        anova.df <- data.frame(term = row.names(anova.df), anova.df)
        anova.df$constr.chi <- summary(cap.env)$constr.chi
        anova.df$unconst.chi <- summary(cap.env)$unconst.chi
        anova.df$tot.chi <- summary(cap.env)$tot.chi
        an.list[[k]] <- anova.df
        
        #calculate centroid distances - if factor has more than 2 levels, take the mean distance
        select.factors <- c("logLoc","topo")
        dist.list[[k]] <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
        
        #do variance partitioning
        endo.var <- varpart(vegdist(mat.otu.pres, method = 'euclidean'),
                            ~ as.matrix(envVars.to),
                            ~ envVars.tmp$topo,
                            ~ envVars.tmp$logLoc)
        endo.var$part$indfract[1:3,]
        
        var.fract <- endo.var$part$indfract[1:3,]
        terms <- c("woodtraits","topo","logLoc")
        var.fract$terms <- terms
        var.fract.list[[k]] <- var.fract
        
        
      }
      
    }else{
      otuprep.name <- paste0(names(data.otuPreplist)[i], "pa")
      for(k in 1:length(data)){
        unique(rd$otuprep)
        terms <- extract_rhs.terms(rd = rd, 
                                   table.name = table.names[k], 
                                   otuprep.name = otuprep.name)
        mat.otu.pres <- data[[k]][["otu"]]
        envVars <- data[[k]][["meta"]]
        
        # trait data for species4 == PLOC is missing (only in d09), so exclude it
        envVars %>%
          filter(species4 != "PLOC") -> envVars
        mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
        
        # remove logLoc == "mush"
        envVars %>%
          filter(logLoc == "mush") -> remove
        mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
        envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
        
        # set up the model
        envVars %>%
          left_join(woodTraits) %>%
          select(terms) -> envVars.to
        envVars.to.s <- scale(envVars.to)
        envVars %>%
          select(logLoc, topo) -> add.envVars
        envVars.tmp <- cbind(envVars.to.s, add.envVars)
        cap.env <-  capscale(mat.otu.pres ~ ., 
                             data = envVars.tmp, 
                             distance = 'jaccard', binary = TRUE)
        
        # anova
        anova.df <- anova(cap.env, by = "terms")
        anova.df <- data.frame(term = row.names(anova.df), anova.df)
        anova.df$constr.chi <- summary(cap.env)$constr.chi
        anova.df$unconst.chi <- summary(cap.env)$unconst.chi
        anova.df$tot.chi <- summary(cap.env)$tot.chi
        an.list[[k]] <- anova.df
        
        #calculate centroid distances - if factor has more than 2 levels, take the mean distance
        select.factors <- c("logLoc","topo")
        dist.list[[k]] <- calcCentroidDist(cap.env = cap.env, FACTORS = select.factors)
        
        #do variance partitioning
        endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
                            ~ as.matrix(envVars.to),
                            ~ envVars.tmp$topo,
                            ~ envVars.tmp$logLoc)
        endo.var$part$indfract[1:3,]
        
        var.fract <- endo.var$part$indfract[1:3,]
        terms <- c("woodtraits","topo","logLoc")
        var.fract$terms <- terms
        var.fract.list[[k]] <- var.fract
        
      
    }
    }
    names(an.list) <- names(dist.list) <- names(var.fract.list) <- table.names
    
    tmp <- list_to_df(an.list)
    tmp %>%
      rename('subset'='source') -> an.list2[[i]]
    tmp <- list_to_df(dist.list)
    tmp %>%
      rename('subset'='source') -> dist.list2[[i]]
    tmp <- list_to_df(var.fract.list)
    tmp %>%
      rename('subset'='source') -> varpart.list2[[i]]
  }
  names(an.list2) <- names(dist.list2) <- names(varpart.list2) <- names(data.otuPreplist)
  an.list2$c100 %>%
    dplyr::rename('SumOfSqs'='Variance') -> an.list2$c100
  tmp <- list_to_df(an.list2)
  tmp %>%
    rename('otuprep'='source') -> an.df
  tmp <- list_to_df(varpart.list2)
  tmp %>%
    rename('otuprep'='source') -> varpart.df
  tmp <- list_to_df(dist.list2)
  tmp %>%
    rename('otuprep'='source') -> dist.df
  
  result <- list(an.df = an.df, varpart.df = varpart.df, dist.df = dist.df)
  return(result)
}


# fungal richness ~ bacterial richness #-----------------------------------------------#

richnessFB_df <- function(data){
  
  # identify the unique yrharv and subset tables
  df.tmp <- data.frame(table.name = names(data))
  df.tmp %>%
    separate(table.name, into = c("yearGene","subset"), sep = "_", remove = F) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> df.tmp
  head(df.tmp)
  YRSUB <- unique(df.tmp$yearSubset)
  rich.df.w <- list()
  for(i in 1:length(YRSUB)){
    
    #fungi
    table <- as.character(df.tmp[df.tmp$yearSubset == YRSUB[i] & df.tmp$gene == "ITS","table.name"])
    richness <- apply(data[[table]]$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    data[[table]]$meta %>%
      left_join(richness.df) -> fung.rich.df
    
    #bacteria
    table <- as.character(df.tmp[df.tmp$yearSubset == YRSUB[i] & df.tmp$gene == "16S","table.name"])
    richness <- apply(data[[table]]$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    data[[table]]$meta %>%
      left_join(richness.df) -> bact.rich.df
    
    rich.df <- rbind(fung.rich.df, bact.rich.df)
    rich.df %>%
      select(gene, yrharv, yrdeploy, tidySeq_sampleName, rich, logLoc, topo, species4, angio.gymno) %>%
      spread(key = gene, value = rich) %>%
      rename('rich.bact' = '16S',
             'rich.fung' = 'ITS') %>%
      mutate(topoLogloc = paste(topo, logLoc, sep = "_")) -> rich.df.w[[i]]
  }
  names(rich.df.w) <- YRSUB
  rich.df <- list_to_df(rich.df.w)
  rich.df %>%
    rename('yearSubset' = 'source') -> rich.df
  
  return(rich.df)
}

richnessFB_tests <- function(data){
  
  # identify the unique yrharv and subset tables
  df.tmp <- data.frame(table.name = names(data))
  df.tmp %>%
    separate(table.name, into = c("yearGene","subset"), sep = "_", remove = F) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) -> df.tmp
  YRSUB <- unique(df.tmp$yearSubset)
  an.list <- list()
  for(i in 1:length(YRSUB)){
    
    #fungi
    table <- as.character(df.tmp[df.tmp$yearSubset == YRSUB[i] & df.tmp$gene == "ITS","table.name"])
    richness <- apply(data[[table]]$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    data[[table]]$meta %>%
      left_join(richness.df) -> fung.rich.df
    
    #bacteria
    table <- as.character(df.tmp[df.tmp$yearSubset == YRSUB[i] & df.tmp$gene == "16S","table.name"])
    richness <- apply(data[[table]]$otu, MARGIN = 1,function(x) sum(x>0))
    richness.df <- data.frame(tidySeq_sampleName = names(richness), rich = richness, stringsAsFactors = F)
    data[[table]]$meta %>%
      left_join(richness.df) -> bact.rich.df
    
    rich.df <- rbind(fung.rich.df, bact.rich.df)
    rich.df %>%
      select(gene, yrharv, yrdeploy, tidySeq_sampleName, rich, logLoc, topo, species4, angio.gymno) %>%
      spread(key = gene, value = rich) %>%
      rename('rich.bact' = '16S',
             'rich.fung' = 'ITS') %>%
      mutate(topoLogloc = paste(topo, logLoc, sep = "_")) -> rich.df.w
    
    # tests
    # fungi
    mod.fung <- lm(rich.fung ~ rich.bact + logLoc + topo + species4, data = rich.df.w)
    an.df.fung <- data.frame(resp = "fung.rich",
                             term = row.names(anova(mod.fung)), 
                             anova(mod.fung), 
                             adj.r2 = summary(mod.fung)$adj.r.squared)
    # bacteria
    mod.bact <- lm(rich.bact ~ rich.fung + logLoc + topo + species4, data = rich.df.w)
    an.df.bact <- data.frame(resp = "bact.rich",
                             term = row.names(anova(mod.bact)), 
                             anova(mod.bact), 
                             adj.r2 = summary(mod.bact)$adj.r.squared)
    an.list[[i]] <- rbind(an.df.fung, an.df.bact)
  }
  names(an.list) <- YRSUB
  an.df <- list_to_df(an.list)
  an.df %>%
    rename('yearSubset' = 'source') -> an.df
  
  return(an.df)
}


# mantel test for fungal composition ~ bacterial composition #-----------------------------------------------#

mantel_test_jaccard <- function(commY, commX){
    
    #make sure commX and commY have the same samples
    samelist <- sameSamples(commY = commY, commX = commX)
    commy.otu <- samelist$commy.otu
    commx.otu <- samelist$commx.otu
    
    #does predictor matrix effect community composition?
    dist.y <- vegdist(commy.otu, method = 'jaccard', binary = T)
    dist.x <- vegdist(commx.otu, method = 'jaccard', binary = T)
    mantel.stat <- mantel(dist.y, dist.x)
    manteldf <- data.frame(r = mantel.stat$statistic, p = mantel.stat$signif , n = dim(commy.otu)[1])
    
    return(manteldf)
}

mantel_test_euclid <- function(commY, commX){
  
  #make sure commX and commY have the same samples
  samelist <- sameSamples(commY = commY, commX = commX)
  commy.otu <- samelist$commy.otu
  commx.otu <- samelist$commx.otu
  
  #does predictor matrix effect community composition?
  dist.y <- vegdist(commy.otu, method = 'euclidean')
  dist.x <- vegdist(commx.otu, method = 'euclidean')
  mantel.stat <- mantel(dist.y, dist.x)
  manteldf <- data.frame(r = mantel.stat$statistic, p = mantel.stat$signif , n = dim(commy.otu)[1])
  
  return(manteldf)
}

batch_write_manteltests_rotplots <- function(data.overlap.choosePrep, data.deploy.choosePrep, presAbs){
  
  if(presAbs == TRUE){
    r.2012_o <- mantel_test_jaccard(commY = data.overlap.choosePrep[['2012ITS_o']], 
                                    commX = data.overlap.choosePrep[['201216S_o']])
    r.2014_o <- mantel_test_jaccard(commY = data.overlap.choosePrep[['2014ITS_o']], 
                                    commX = data.overlap.choosePrep[['201416S_o']])
    r.2012_d11 <- mantel_test_jaccard(commY = data.deploy.choosePrep[['2012ITS_d11']], 
                                      commX = data.deploy.choosePrep[['201216S_d11']])
    r.2012_d09 <- mantel_test_jaccard(commY = data.deploy.choosePrep[['2012ITS_d09']], 
                                      commX = data.deploy.choosePrep[['201216S_d09']])
    r.2014_d11 <- mantel_test_jaccard(commY = data.deploy.choosePrep[['2014ITS_d11']], 
                                      commX = data.deploy.choosePrep[['201416S_d11']])
    r.2014_d09 <- mantel_test_jaccard(commY = data.deploy.choosePrep[['2014ITS_d09']], 
                                      commX = data.deploy.choosePrep[['201416S_d09']])
  }else{
    r.2012_o <- mantel_test_euclid(commY = data.overlap.choosePrep[['2012ITS_o']], 
                                    commX = data.overlap.choosePrep[['201216S_o']])
    r.2014_o <- mantel_test_euclid(commY = data.overlap.choosePrep[['2014ITS_o']], 
                                    commX = data.overlap.choosePrep[['201416S_o']])
    r.2012_d11 <- mantel_test_euclid(commY = data.deploy.choosePrep[['2012ITS_d11']], 
                                      commX = data.deploy.choosePrep[['201216S_d11']])
    r.2012_d09 <- mantel_test_euclid(commY = data.deploy.choosePrep[['2012ITS_d09']], 
                                      commX = data.deploy.choosePrep[['201216S_d09']])
    r.2014_d11 <- mantel_test_euclid(commY = data.deploy.choosePrep[['2014ITS_d11']], 
                                      commX = data.deploy.choosePrep[['201416S_d11']])
    r.2014_d09 <- mantel_test_euclid(commY = data.deploy.choosePrep[['2014ITS_d09']], 
                                      commX = data.deploy.choosePrep[['201416S_d09']])
  }
  
  
  r.df <- rbind(r.2012_o, r.2014_o,
                r.2012_d11, r.2012_d09,
                r.2014_d11, r.2014_d09)
  r.df$yearSubset <- c("2012_o","2014_o",
                       "2012_d11","2012_d09",
                       "2014_d11","2014_d09")
  return(r.df)
  
}

# composition ~ co-occurring clade composition + study factors #-----------------------------------------------#

# this takes forever, it needs to happen on the cluster
modelSelectionMDS_jaccard <- function(commx.otu){
    
    # represent the community X matrix with MDS axes
    commx.pc <- capscale(commx.otu ~ 1, distance = 'jaccard', binary = T)
    lastMDS <- dim(commx.otu)[1] - 1
    commx.scores <- scores(commx.pc, choices = 1:lastMDS, display = "sites")
    
    # figure out which MDS axes are needed to represent community X using model selection
    cap.pc <- capscale(commx.otu ~ ., data=data.frame(commx.scores), distance = 'jaccard', binary = T) # full case
    
    mod0.pc <- capscale(commx.otu ~ 1, data=data.frame(commx.scores), distance = 'jaccard', binary = T) # set up the null case
    step.pc <- ordistep(mod0.pc, scope=formula(cap.pc)) # model selection - this takes forever
    #step.pc$anova # only use the significant axes
    
    return(step.pc)
}

#commY <- data.deploy.t100pa[['2012ITS_d11']]
#commX <- data.deploy.t100pa[['201216S_d11']]
#commX.signifaxes = c(1,2)
dbrda_varpart_studyfactorsVScomm <- function(commY, commY.name, commX, commX.name, commX.signifaxes){
  
  #make sure the communities have the same number of samples
  same.list <- sameSamples(commY, commX)
  commy.otu <- same.list$commy.otu
  commx.otu <- same.list$commx.otu
  
  #identify the environmental data
  meta <- commY$meta
  colnames(meta)
  meta %>%
    filter(tidySeq_sampleName %in% row.names(commy.otu)) %>%
    select(tidySeq_sampleName, logLoc, species4, topo) -> envVars
  
  #identify the significant predictor community axes
  commx.pc <- capscale(commx.otu ~ 1, distance = 'jaccard', binary = T)
  commx.sub <- scores(commx.pc, choices = commX.signifaxes, display = "sites")
  
  # partition variation among 2 predictor tables:
  #   1) study factors : 
  #   2) significant community x axes
  endo.var <- varpart(vegdist(commy.otu, distance = 'jaccard', binary = T), 
                      ~ logLoc + species4 + topo,
                      commx.sub, data=envVars)
  var.fract <- endo.var$part$fract[1:2,]
  var.fract.ind <- endo.var$part$indfract[c(1,3,2,4),]
  var.fract <- rbind(var.fract, var.fract.ind)
  terms <- c("env","commx","env","commx", "both","residuals")
  fractType <- c("trad","trad","ind","ind","common","residuals")
  var.fract$terms <- terms
  var.fract$fractType <- fractType
  var.fract$pval <- NA
  var.fract
  # significance of X1
  an10<-anova(capscale(vegdist(commy.otu, distance='bray')	~ logLoc + species4 + topo,
                       data=envVars))
  var.fract[var.fract$terms == "env" & var.fract$fractType == "trad","pval"] <- an10$`Pr(>F)`[1]
  
  # significance of X2
  an20<-anova(capscale(vegdist(commy.otu, distance='bray')	~ commx.sub,
                       data=envVars))
  var.fract[var.fract$terms == "commx" & var.fract$fractType == "trad","pval"] <- an20$`Pr(>F)`[1]
  
  
  # significance of X1|X2
  an1<-anova(capscale(vegdist(commy.otu, distance='bray')	~ logLoc + species4 + topo
                      + Condition(commx.sub),
                      data=envVars))
  var.fract[var.fract$terms == "env" & var.fract$fractType == "ind","pval"] <- an1$`Pr(>F)`[1]
  
  # significance of X2|X1
  an2<-anova(capscale(vegdist(commy.otu, distance='bray')	~ commx.sub
                      + Condition(logLoc) + Condition(species4) + Condition(topo),
                      data=envVars))
  var.fract[var.fract$terms == "commx" & var.fract$fractType == "ind","pval"] <- an2$`Pr(>F)`[1]
  var.fract
  # add community names
  var.fract$commY.name <- commY.name
  var.fract$commX.name <- commX.name
  
  return(var.fract)
}

dbrda_varpart_studyfactorsVScomm_clr <- function(commY, commY.name, commX, commX.name, commX.signifaxes,
                                                 presAbs){
  
  #make sure the communities have the same number of samples
  same.list <- sameSamples(commY, commX)
  commy.otu <- same.list$commy.otu
  commx.otu <- same.list$commx.otu
  
  #identify the environmental data
  meta <- commY$meta
  colnames(meta)
  meta %>%
    filter(tidySeq_sampleName %in% row.names(commy.otu)) %>%
    select(tidySeq_sampleName, logLoc, species4, topo) -> envVars
  
  #identify the significant predictor community axes
  if(presAbs == TRUE){
    commx.pc <- capscale(commx.otu ~ 1, distance = 'jaccard', binary = T)
    commx.sub <- scores(commx.pc, choices = commX.signifaxes, display = "sites")
    
    # partition variation among 2 predictor tables:
    #   1) study factors : 
    #   2) significant community x axes
    endo.var <- varpart(vegdist(commy.otu, distance = 'jaccard', binary = T), 
                        ~ logLoc + species4 + topo,
                        commx.sub, data=envVars)
    var.fract <- endo.var$part$fract[1:2,]
    var.fract.ind <- endo.var$part$indfract[c(1,3,2,4),]
    var.fract <- rbind(var.fract, var.fract.ind)
    terms <- c("env","commx","env","commx", "both","residuals")
    fractType <- c("trad","trad","ind","ind","common","residuals")
    var.fract$terms <- terms
    var.fract$fractType <- fractType
    var.fract$pval <- NA
    var.fract
    # significance of X1
    an10<-anova(capscale(vegdist(commy.otu, distance='bray')	~ logLoc + species4 + topo,
                         data=envVars))
    var.fract[var.fract$terms == "env" & var.fract$fractType == "trad","pval"] <- an10$`Pr(>F)`[1]
    
    # significance of X2
    an20<-anova(capscale(vegdist(commy.otu, distance='bray')	~ commx.sub,
                         data=envVars))
    var.fract[var.fract$terms == "commx" & var.fract$fractType == "trad","pval"] <- an20$`Pr(>F)`[1]
    
    
    # significance of X1|X2
    an1<-anova(capscale(vegdist(commy.otu, distance='bray')	~ logLoc + species4 + topo
                        + Condition(commx.sub),
                        data=envVars))
    var.fract[var.fract$terms == "env" & var.fract$fractType == "ind","pval"] <- an1$`Pr(>F)`[1]
    
    # significance of X2|X1
    an2<-anova(capscale(vegdist(commy.otu, distance='bray')	~ commx.sub
                        + Condition(logLoc) + Condition(species4) + Condition(topo),
                        data=envVars))
    var.fract[var.fract$terms == "commx" & var.fract$fractType == "ind","pval"] <- an2$`Pr(>F)`[1]
    var.fract
    # add community names
    var.fract$commY.name <- commY.name
    var.fract$commX.name <- commX.name
  }else{
    commx.pc <- capscale(commx.otu ~ 1, method = 'euclidean')
    commx.sub <- scores(commx.pc, choices = commX.signifaxes, display = "sites")
    
    # partition variation among 2 predictor tables:
    #   1) study factors : 
    #   2) significant community x axes
    endo.var <- varpart(vegdist(commy.otu, method = 'euclidean'), 
                        ~ logLoc + species4 + topo,
                        commx.sub, data=envVars)
    var.fract <- endo.var$part$fract[1:2,]
    var.fract.ind <- endo.var$part$indfract[c(1,3,2,4),]
    var.fract <- rbind(var.fract, var.fract.ind)
    terms <- c("env","commx","env","commx", "both","residuals")
    fractType <- c("trad","trad","ind","ind","common","residuals")
    var.fract$terms <- terms
    var.fract$fractType <- fractType
    var.fract$pval <- NA
    var.fract
    # significance of X1
    an10<-anova(capscale(vegdist(commy.otu, method = 'euclidean')	~ logLoc + species4 + topo,
                         data=envVars))
    var.fract[var.fract$terms == "env" & var.fract$fractType == "trad","pval"] <- an10$`Pr(>F)`[1]
    
    # significance of X2
    an20<-anova(capscale(vegdist(commy.otu, method = 'euclidean')	~ commx.sub,
                         data=envVars))
    var.fract[var.fract$terms == "commx" & var.fract$fractType == "trad","pval"] <- an20$`Pr(>F)`[1]
    
    
    # significance of X1|X2
    an1<-anova(capscale(vegdist(commy.otu, method = 'euclidean')	~ logLoc + species4 + topo
                        + Condition(commx.sub),
                        data=envVars))
    var.fract[var.fract$terms == "env" & var.fract$fractType == "ind","pval"] <- an1$`Pr(>F)`[1]
    
    # significance of X2|X1
    an2<-anova(capscale(vegdist(commy.otu, method = 'euclidean')	~ commx.sub
                        + Condition(logLoc) + Condition(species4) + Condition(topo),
                        data=envVars))
    var.fract[var.fract$terms == "commx" & var.fract$fractType == "ind","pval"] <- an2$`Pr(>F)`[1]
    var.fract
    # add community names
    var.fract$commY.name <- commY.name
    var.fract$commX.name <- commX.name
  }
  
  
  return(var.fract)
}

#tab.names <- names(data.deploy.t100pa)
create_fung_bact_pairs_rotplot <- function(tab.names){
  
  indx <- data.frame(tab = tab.names)
  indx %>%
    separate(tab, into = c("yearGene","subset"), remove = F) %>%
    separate(yearGene, into = c("year","gene"), sep = 4) %>%
    mutate(yearSubset = paste(year, subset,  sep = "_")) %>%
    select(tab, gene) -> indx
  
  fung.comm <- c("ITS")
  bact.comm <- c("16S")
  tmp <- expand.grid(fung.comm, bact.comm)
  tmp %>%
    mutate(col1 = paste(Var1, Var2, sep = "_<-_")) %>%
    mutate(col2 = paste(Var2, Var1, sep = "_<-_")) -> tmp
  pair.vec <- c(tmp$col1, tmp$col2)
  
  pair.df <- data.frame(pair.vec = pair.vec)
  pair.df %>%
    separate(pair.vec, into = c("gene","commX.name"), sep = "_<-_", remove = FALSE) %>%
    left_join(indx) %>%
    rename('commY.tab'='tab',
           'commY.name'='gene',
           'gene'='commX.name') %>%
    separate(commY.tab, into = c("yearGene","subset"), remove = F) %>%
    separate(yearGene, into = c("year","geneY"), sep = 4) %>%
    mutate(commX.tab = paste0(year, gene, "_", subset)) %>%
    rename('commX.name'='gene') %>%
    mutate(yearSubset = paste(year, subset, sep = "_")) %>%
    select(pair.vec, yearSubset, commY.name, commX.name, commY.tab, commX.tab) -> pair.df
  
  return(pair.df)
}

batch_varpart_commx_rotplot <- function(pair.df, commX.signifaxes.list, datasets, presAbs){
  
  # pair.df
  # datasets <- curr.data
  # commX.signifaxes.list
  # presAbs <- FALSE

  var.fract.list <- list()
  for(i in 1:dim(pair.df)[1]){
    commY.name <- pair.df[i,"commY.name"]
    commX.name <- pair.df[i,"commX.name"]
    commY.tab <- as.character(pair.df[i,"commY.tab"])
    commX.tab <- as.character(pair.df[i,"commX.tab"])
    commX.signifaxes <- commX.signifaxes.list[[i]]
    var.fract.df <- dbrda_varpart_studyfactorsVScomm_clr(commY = datasets[[commY.tab]], commY.name, 
                                                            commX = datasets[[commX.tab]], commX.name,
                                                            commX.signifaxes, presAbs = presAbs)
    var.fract.df$yearSubset <- pair.df[i,'yearSubset']
    var.fract.list[[i]] <- var.fract.df
  }
  var.fract.df <- list_to_df(var.fract.list)
  return(var.fract.df)
  
}



