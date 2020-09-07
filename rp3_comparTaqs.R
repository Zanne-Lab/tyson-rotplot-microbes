# rp3_compareTaqs.R

# Different taqs were used for the 2012 and 2014 sequencing runs
# What was the effect of the different taqs on...
# sampling depth (# reads per sample), 
# OTU richness (# OTUs per sample), and 
# OTU composition (distance between pres/abs OTU tables)

#----------------------------------------------------------#
# Compare reads per sample

  r <- make_taqtest_objs(samp.taqredos, samp.meta, tab.t100, taqredo.t100)
  tab.taqcomp <- r$otulist
  samp.indx <- r$sampdf
  
  # calculate the total number of reads per sample
  tab.r <- lapply(tab.taqcomp, function(x){
    vec <- rowSums(x)
    df <- data.frame(tidySeq_sampleName = names(vec), resp = vec)
    return(df)
  })
  df.r <- list_to_df(tab.r)
  df.r %>%
    rename('table'='source') %>%
    left_join(samp.indx) -> df.r
  
  yLab <- "Reads per sample"
  
  # test differences 
  library(lme4)
  library(lmerTest)
  #install.packages("pbkrtest")
  
  # taq differences
  r.taq <- test_taqdiffs(df.r, yLab)
  anova(r.taq$mod.12, ddf="Kenward-Roger") # ns
  anova(r.taq$mod.14, ddf="Kenward-Roger") #ns
  r.taq$p
  ggsave(r.taq$p, file = paste0(output.path, "sampreads_byTaq.pdf"), width = 6, height = 3)
  
  # yrseq differences for Go taq
  r.yrseq_gotaq <- test_yrseq_gotaq(df.r, yLab)
  anova(r.yrseq_gotaq$mod, ddf="Kenward-Roger") #ns
  
  # yrseq differences for Q taq
  r.yrseq_qtaq <- test_yrseq_qtaq(df.r, yLab)
  anova(r.yrseq_qtaq$mod, ddf="Kenward-Roger") # signif!
  pdf(file = paste0(output.path, "sampreads_byYrSeq.pdf"), width = 6, height = 3)
  grid.arrange(
    r.yrseq_gotaq$p + ggtitle("GoTaq") + ylim(c(0,50000)),
    r.yrseq_qtaq$p + ggtitle("Q5") + ylim(c(0,50000)),
    ncol = 2)
  dev.off()
  
#----------------------------------------------------------#
# Compare sample OTU richness

  # calculate OTU richness per sample
  tab.r <- lapply(tab.taqcomp, function(x){
    x.pa <- 1 * (x > 0) # transform to presence/absence
    vec <- rowSums(x.pa) # calculate sample OTU richness
    df <- data.frame(tidySeq_sampleName = names(vec), resp = vec)
    return(df)
  })
  df.r <- list_to_df(tab.r)
  df.r %>%
    rename('table'='source') %>%
    left_join(samp.indx) -> df.r
  
  yLab <- "OTU richness per sample"
  
  # test differences 
  
  # taq differences
  r.taq <- test_taqdiffs(df.r, yLab)
  anova(r.taq$mod.12, ddf="Kenward-Roger") # signif more otus with Q5
  anova(r.taq$mod.14, ddf="Kenward-Roger") # signif more otus with Q5
  r.taq$p
  ggsave(r.taq$p, file = paste0(output.path, "sampotus_byTaq.pdf"), width = 6, height = 3)
  
  # yrseq differences for Go taq
  r.yrseq_gotaq <- test_yrseq_gotaq(df.r, yLab)
  anova(r.yrseq_gotaq$mod, ddf="Kenward-Roger") # signif fewer otus in 2017
  
  # yrseq differences for Q taq
  r.yrseq_qtaq <- test_yrseq_qtaq(df.r, yLab)
  anova(r.yrseq_qtaq$mod, ddf="Kenward-Roger") # signif more otus in 2017
  pdf(file = paste0(output.path, "sampotus_byYrSeq.pdf"), width = 6, height = 3)
  grid.arrange(
    r.yrseq_gotaq$p + ggtitle("GoTaq") + ylim(c(0,200)),
    r.yrseq_qtaq$p + ggtitle("Q5") + ylim(c(0,200)),
    ncol = 2)
  dev.off()
  
  
#----------------------------------------------------------#
# OTU composition - taq differences

  # taq differences
  tab.17 <- tab.taqcomp[['2017ITS']]
  tab.17 <- 1 * (tab.17 > 0) # transform to presence/absence
  samp.indx %>%
    filter(yrseq == 2017) -> samp.17
  
  r <- test_taqdiffs_composition(mat.otu.pres = tab.17, 
                                 envVars = samp.17, 
                                 logTrt)
  r$anova.df
  r$p
  ggsave(r$p, file = paste0(output.path, "otucomp_byTaq.pdf"), width = 6, height = 3)
