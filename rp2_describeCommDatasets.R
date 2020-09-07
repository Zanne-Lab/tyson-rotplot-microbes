# rp2_describeCommDatasets.R

#----------------------------------------------------------#
# Reads per sample in each sequencing run
tab.t100 <- otuPreplist[['t100']]

result.t100 <- lapply(tab.t100, function(x){
  x <- x$otu
  vec <- rowSums(x) # calc the total number of reads per sample
  df <- data.frame(min = range(vec)[1],
                   max = range(vec)[2],
                   mean = mean(vec),
                   n = length(vec),
                   se = sd(vec)/sqrt(length(vec)))
  return(df)
})
result.t100 <- list_to_df(result.t100)
write.csv(result.t100, file = paste0(output.path, "otuTable_numReads.csv"))



#----------------------------------------------------------#
# Samples and OTUs per normalization procedure

# For each sequencing run and OTU normalization procedure report... 
# (a) the total number of samples, 
# (b) the total number of OTUS, 
# (c) the min, max, mean, and se for sample OTU richness

result <- lapply(otuPreplist, function(x){
  lapply(x, function(k){
    vec <- rowSums(k$otu) # calc the number of OTUs per sample
    df <- data.frame(
      nSamps = dim(k$otu)[1], # calc the total number of samples
      nOTUs = dim(k$otu)[2], # calc the total number of OTUs
      minOTUs.persamp = range(vec)[1],
      maxOTUs.persamp = range(vec)[2],
      meanOTUs.persamp = mean(vec),
      seOTUs.persamp = sd(vec)/sqrt(length(vec))
      )
    return(df)
  })
})
df <- list_to_df(
  lapply(result, function(k){
    df <- list_to_df(k)
    df %>%
      rename('table'='source') -> df
    })
)
df %>%
  rename('otuprep'='source') %>%
  mutate(table = recode_factor(table, `2012ITS` ='2012ITS', `2014ITS` = '2014ITS', `201216S` = '201216S',`201416S` = '201416S')) %>%
  mutate(otuprep = recode_factor(otuprep, t100 = 't100', r100 = 'r100', r500 = 'r500', r1000 = 'r1000')) %>%
  arrange(table, otuprep) -> df
df
#write.csv(df, file = paste0(output.path, "otuTable_summary.csv"))



#----------------------------------------------------------#
# What proportion of OTUs in each sequencing run are identified to phyla and genus?

result <- lapply(tab.taxa, function(x){
  nOTUs <- dim(x)[1]
  nGenus <- sum(x$genus != "unidentified")
  nFamily <- sum(x$family != "unidentified")
  nOrder <- sum(x$order != "unidentified")
  nClass <- sum(x$class != "unidentified")
  nPhylum <- sum(x$phylum != "unidentified")
  nKingdom <- sum(x$kingdom != "unidentified")
  df <- data.frame(nOTUs = nOTUs,
             percGenus = nGenus / nOTUs,
             percFamily = nFamily / nOTUs,
             percOrder = nOrder / nOTUs,
             percClass = nClass / nOTUs,
             percPhylum = nPhylum / nOTUs,
             percKingdom = nKingdom / nOTUs)
  return(df)
})
df <- list_to_df(result)
df %>%
  rename('table'='source') %>%
  mutate(table = recode_factor(table, 
                               `2012ITS` ='2012ITS', 
                               `2014ITS` = '2014ITS', 
                               `201216S` = '201216S',
                               `201416S` = '201416S')) %>%
  arrange(table) -> df
df
write.csv(df, file = paste0(output.path, "otuTable_taxonids.csv"))


#----------------------------------------------------------#
# What proportion of fungal OTUs in each sequencing run are identified to functional guild?

tab.taxa.its <- tab.taxa[grepl("ITS", names(tab.taxa))]
result <- lapply(tab.taxa.its, function(x){
  nOTUs <- dim(x)[1]
  nTroph <- sum(!is.na(x$Trophic.Mode))
  nTroph.HP <- sum(!is.na(x$Trophic.Mode) & x$Confidence.Ranking %in% c("Probable", "Highly Probable"))
  df <- data.frame(nOTUs = nOTUs,
             percTroph = nTroph / nOTUs,
             percTroph.ProbHProb = nTroph.HP / nOTUs)
  return(df)
})
df <- list_to_df(result)
df %>%
  rename('table'='source') %>%
  mutate(table = recode_factor(table, `2012ITS` ='2012ITS', `2014ITS` = '2014ITS', `201216S` = '201216S',`201416S` = '201416S')) %>%
  arrange(table) -> df
df
write.csv(df, file = paste0(output.path, "otuTable_guildids.csv"))


#----------------------------------------------------------#
# Frequency of OTUs by phylum in each sequencing run

r <- freqPhyla(tab.otu = tab.t100, tab.taxa = tab.taxa)
mytheme <- make_ggplot_theme()

pdf(file = paste0(output.path, "freqPhyla.pdf"), width = 6, height = 6)
grid.arrange(r$p.fung + ggtitle("Fungi") + mytheme,
             r$p.bact + ggtitle("Bacteria") + mytheme)
dev.off()


#----------------------------------------------------------#
# 10 most cosmopolitan OTUs in each sequencing run

topten.list <- list()
for(i in 1:length(unique(r$df$table))){
  r$df %>%
    filter(cat == "Highly cosmopolitan") %>%
    filter(table == unique(r$df$table)[i]) %>%
    arrange(desc(occur)) -> tmp
  topten.list[[i]] <- tmp[1:10,]
}
names(topten.list) <- unique(r$df$table)
topten.df <- list_to_df(topten.list)
topten.df %>%
  select(gene, yrseq, OTUid, occur, genusSpecies, phylum, class, order, family) -> topten.df
write.csv(topten.df, file = paste0(output.path, "otuSummary_topten.csv"))


