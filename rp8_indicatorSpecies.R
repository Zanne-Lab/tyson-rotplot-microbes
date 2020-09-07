# rp8_indicatorSpecies.R

#----------------------------------------------------------#
#Set up

#libraries
library("tidyverse") # ggplot2, tidyr, readr, dplyr
library("vegan") #rarecurve(), ordistep(), etc
library("gridExtra") #grid.arrange()
library("labdsv")

#functions
source("code/helpers.R")
source("code/load_woodTraits.R")
source("code/load_woodPhylo.R")
sourceDir("code/microbes")
sourceDir("code/rotplot_microbes")

#paths
intermed.path <- "data_Rsynth/rotplot_intermediates/"
output.path <- "output/rotplot_microbes/"

#----------------------------------------------------------#
#Load OTU tables
overlap.otuPreplist <- readRDS(file = paste0(intermed.path, "overlap_otuPreplist.RData"))
deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))

#Load taxonomy tables
tab.taxa <- read_tidy_taxon(table.names = c("2012ITS","2014ITS","201216S","201416S"))
# create an informative uniqueID for each OTU
TABLE <- names(tab.taxa)
for(i in 1:length(TABLE)){
  tab.taxa[[i]] <- add_niceOTUname(curr.tab.taxa = tab.taxa[[TABLE[i]]])
}
lapply(tab.taxa, colnames)

#Load wood traits
woodTraits <- readRDS(file = "data_Rsynth/rotplot_intermediates/woodTraits.RData")
# pull out just the wood traits to include in indicator species tests
woodTraits %>%
  select(species4, Mn, lignin) %>%
  filter(!is.na(Mn)) -> select.traits
hist(select.traits$lignin)
hist(select.traits$Mn)
select.traits %>%
  mutate(lignin.cat2 = cut(as.numeric(lignin),2,labels=FALSE)) %>%
  mutate(Mn.cat2 = cut(as.numeric(Mn),2,labels=FALSE)) -> select.traits

#----------------------------------------------------------#

# create function to do the analysis
do_IndSp_w.r500pa <- function(otuPreplist, 
                              selection, select.traits, overlap){
  
  #load data
  curr.otu <- otuPreplist$r500[[selection]]$otu # presence-absence
  curr.meta <- otuPreplist$r500[[selection]]$meta
  
  # add wood traits to the meta table
  curr.meta %>%
    left_join(select.traits) -> curr.meta
  # remove OTU cols that do not appear in any samples
  curr.otu <- curr.otu[,colSums(curr.otu) != 0]
  dim(curr.otu)
  
  # make look up table for levels
  topo.levels <- levels(factor(curr.meta$topo))
  topo.df <- data.frame(level = topo.levels, 
             num = seq(1:length(topo.levels)))
  logLoc.levels <- levels(factor(curr.meta$logLoc))
  logLoc.df <- data.frame(level = logLoc.levels, 
                        num = seq(1:length(logLoc.levels)))
  species.levels <- levels(factor(curr.meta$species4))
  species.df <- data.frame(level = species.levels, 
                          num = seq(1:length(species.levels)))
  ag.levels <- levels(factor(curr.meta$angio.gymno))
  ag.df <- data.frame(level = ag.levels, 
                           num = seq(1:length(ag.levels)))
  lookup.list <- list(topo = topo.df,
       logLoc = logLoc.df,
       species = species.df,
       ag = ag.df)
  
  # do indicator species analysis
  r.topo <- indval(curr.otu, as.numeric(factor(curr.meta$topo)))
  r.logLoc <- indval(curr.otu, as.numeric(factor(curr.meta$logLoc)))
  r.species <- indval(curr.otu, as.numeric(factor(curr.meta$species4)))
  r.ag <- indval(curr.otu, as.numeric(factor(curr.meta$angio.gymno)))
  r.lignin <- indval(curr.otu, curr.meta$lignin.cat2)
  r.Mn <- indval(curr.otu, curr.meta$Mn.cat2)
  indsp <- list(r.topo = r.topo, 
                r.logLoc = r.logLoc, 
                r.species = r.species, 
                r.ag = r.ag, 
                r.lignin = r.lignin, 
                r.Mn = r.Mn)
  
  if(overlap == TRUE){
    r.yrdeploy <- indval(curr.otu, as.numeric(factor(curr.meta$yrdeploy)))
    indsp <- list(r.topo = r.topo, 
                  r.logLoc = r.logLoc, 
                  r.species = r.species, 
                  r.ag = r.ag, 
                  r.lignin = r.lignin, 
                  r.Mn = r.Mn,
                  r.yrdeploy = r.yrdeploy)
  }
  
  result <- list(indsp = indsp, 
       lookup = lookup.list)
  
  return(result)
}

# create function to put analysis results into a dataframe
makedf_IndSp <- function(data){
  gradients <- names(data)
  df.list <- list()
  for(i in 1:length(gradients)){
    df <- data.frame(maxcls = data[[gradients[i]]]$maxcls, 
                     indcls = data[[gradients[i]]]$indcls,
                     pval = data[[gradients[i]]]$pval)
    df.list[[i]]<- data.frame(OTUid = row.names(df), df, row.names = NULL)
  }
  names(df.list) <- gradients
  df <- list_to_df(df.list)  
  df %>%
    filter(pval < 0.05) -> df.f
  
  return(df.f)
}

# annotate df with OTU taxonomy
addTax <- function(df, tab.taxa, tax.selection){
  
  curr.tax <- tab.taxa[[tax.selection]]
  colnames(curr.tax)
  curr.tax %>%
    select(OTUid, Trophic.Mode, Confidence.Ranking, best.name) -> curr.tax.tmp
  df %>%
    left_join(curr.tax.tmp) -> tmp
  return(tmp)
}



#----------------------------------------------------------#
# Overlap datasets

# do analyses and make dataframes
# ITS 2012
r <- do_IndSp_w.r500pa(otuPreplist = overlap.otuPreplist,
                  selection = "2012ITS_o",
                  select.traits = select.traits,
                  overlap = TRUE)
#summary(r$indsp$r.Mn) #no indicators
df <- makedf_IndSp(data = r$indsp)
df <- addTax(df = df, tab.taxa = tab.taxa, tax.selection = "2012ITS")
df %>%
  mutate(selection = "2012ITS_o") %>%
  separate(source, into = c("drop","gradient"), sep = "r.") %>%
  mutate(gradient = ifelse(gradient == "y", "yrdeploy", gradient)) %>%
  mutate(maxcls.code = as.character(maxcls)) %>%
  select(-drop) -> tmp

#unique(tmp$gradient)
#r$lookup$topo
tmp[tmp$gradient == "topo","maxcls.code"] <- recode(tmp[tmp$gradient == "topo","maxcls.code"], "1"="H", "2"="L")
#r$lookup$logLoc
tmp[tmp$gradient == "logLoc","maxcls.code"] <- recode(tmp[tmp$gradient == "logLoc","maxcls.code"], "1"="t", "2"="b")
#r$lookup$species
tmp[tmp$gradient == "species","maxcls.code"] <- recode(tmp[tmp$gradient == "species","maxcls.code"], 
                                                       "1"="CEOC", "2"="JUVI", "3"="QUVE")
#r$lookup$ag
tmp[tmp$gradient == "ag","maxcls.code"] <- recode(tmp[tmp$gradient == "ag","maxcls.code"], 
                                                       "1"="Angiosperm", "2"="Gymnosperm")
df.its12 <- tmp


# ITS 2014
r <- do_IndSp_w.r500pa(otuPreplist = overlap.otuPreplist,
                       selection = "2014ITS_o",
                       select.traits = select.traits,
                       overlap = TRUE)
#summary(r$indsp$r.Mn) #no indicators
df <- makedf_IndSp(data = r$indsp)
df <- addTax(df = df, tab.taxa = tab.taxa, tax.selection = "2014ITS")
df %>%
  mutate(selection = "2014ITS_o") %>%
  separate(source, into = c("drop","gradient"), sep = "r.") %>%
  mutate(gradient = ifelse(gradient == "y", "yrdeploy", gradient)) %>%
  mutate(maxcls.code = as.character(maxcls)) %>%
  select(-drop) -> tmp

unique(tmp$gradient)
r$lookup$topo
tmp[tmp$gradient == "topo","maxcls.code"] <- recode(tmp[tmp$gradient == "topo","maxcls.code"], "1"="H", "2"="L")
r$lookup$logLoc
tmp[tmp$gradient == "logLoc","maxcls.code"] <- recode(tmp[tmp$gradient == "logLoc","maxcls.code"], "1"="t", "2"="b","3"="mush")
r$lookup$species
tmp[tmp$gradient == "species","maxcls.code"] <- recode(tmp[tmp$gradient == "species","maxcls.code"], 
                                                       "1"="CEOC", "2"="JUVI", "3"="QUVE")
r$lookup$ag
tmp[tmp$gradient == "ag","maxcls.code"] <- recode(tmp[tmp$gradient == "ag","maxcls.code"], 
                                                  "1"="Angiosperm", "2"="Gymnosperm")
df.its14 <- tmp

#combine
df.its.o <- rbind(df.its12, df.its14)

df.its.o %>%
  group_by(gradient, selection, maxcls.code) %>%
  summarize(n = length(OTUid)) -> tmp

#----------------------------------------------------------#
# plot number of significant OTU indicators per gradient+category

ggplot(tmp, aes(x = maxcls.code, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(gradient ~ selection, scales = "free_x") +
  theme_bw() +
  xlab("Gradient category") + 
  ylab("Number of indicator OTUs (p<0.05")
ggsave(filename = file.path(output.path, "num_IndSp.jpeg"))

# table of indicators

# topo
df.its.o %>%
  filter(gradient == "topo") %>%
  select(gradient, selection, maxcls.code, indcls,
         best.name, Trophic.Mode, Confidence.Ranking, pval) %>%
  arrange(selection, maxcls.code) -> df.its.o.topo
write.csv(df.its.o.topo, file = file.path(output.path, "indSp_topo.csv"))

# logLoc
df.its.o %>%
  filter(gradient == "logLoc") %>%
  select(gradient, selection, maxcls.code, indcls,
         best.name, Trophic.Mode, Confidence.Ranking, pval) %>%
  arrange(selection, maxcls.code) -> df.its.o.logLoc
write.csv(df.its.o.logLoc, file = file.path(output.path, "indSp_logLoc.csv"))


