# rp1_compileCommDatasets.R

# otuPreplist
# overlap.otuPreplist
# deploy.otuPreplist
# taqredos.otuPreplist

#----------------------------------------------------------#
#Set up

#libraries
library("tidyverse") # ggplot2, tidyr, readr, dplyr
library("vegan") #rarecurve(), ordistep(), etc
library("gridExtra") #grid.arrange()

#functions
source("code/helpers.R")
source("code/load_woodTraits.R")
source("code/load_woodPhylo.R")
sourceDir("code/microbes")
sourceDir("code/rotplot_microbes")

#paths
intermed.path <- "data_Rsynth/rotplot_intermediates/"
output.path <- "output/rotplot_microbes/"

# #----------------------------------------------------------#
# #Load OTU tables
# 
# # load and process the raw OTU table
# tab.t100 <- prep_rotplotNobags_otuTabs() # tidy the sample names, subtract negative controls, trim failed samples
# lapply(tab.t100, dim)
# 
# # load taq redos
# taqredo.t100 <- prep_taqRedo_2017ITS() # OTU table only for ITS
# lapply(taqredo.t100, dim)
# 
# #----------------------------------------------------------#
# # Do OTU table normalization
# 
# MakeCLMNorm.OTUtab <- function(curr.otuTab.list){
# 
#   require(compositions)
# 
#   curr.otuTab.list.norm <- list()
#   for(i in 1:length(curr.otuTab.list)){
# 
#     curr.otuTab <- curr.otuTab.list[[i]]
#     asv_clr <- clr(curr.otuTab) # transform
#     asv_clr <- data.frame(asv_clr)
# 
#     curr.otuTab.list.norm[[i]] <- asv_clr
#   }
#   names(curr.otuTab.list.norm) <- names(curr.otuTab.list)
# 
#   return(curr.otuTab.list.norm)
# }
# 
# make_pa <- function(tab){
#   tab.pa <- 1 * (tab > 0)
#   return(tab.pa)
# }
# 
# require("phyloseq") # for rarify
# 
# # do otu table normalizations
# tab.r100 <- MakeRare.OTUtab(curr.otuTab.list = tab.t100, sampleSize = 100)
# tab.r500 <- MakeRare.OTUtab(curr.otuTab.list = tab.t100, sampleSize = 500)
# tab.r1000 <- MakeRare.OTUtab(curr.otuTab.list = tab.t100, sampleSize = 1000)
# tab.clr100 <- MakeCLMNorm.OTUtab(curr.otuTab.list = tab.t100)
# tab.list <- list(c100 = tab.clr100,
#                  r100 = lapply(tab.r100, make_pa),
#                  r500 = lapply(tab.r500, make_pa),
#                  r1000 = lapply(tab.r1000, make_pa))
# lapply(tab.list, function(x) lapply(x, dim))
# 
# 
# # do otu table normalization -- for taq redos
# taqredo.r500 <- MakeRare.OTUtab(curr.otuTab.list = taqredo.t100, sampleSize = 500)
# taqredo.r100 <- MakeRare.OTUtab(curr.otuTab.list = taqredo.t100, sampleSize = 100)
# taqredo.r1000 <- MakeRare.OTUtab(curr.otuTab.list = taqredo.t100, sampleSize = 1000)
# taqredo.clr100 <- MakeCLMNorm.OTUtab(curr.otuTab.list = taqredo.t100)
# taqredo.list <- list(c100 = taqredo.clr100,
#                      r100 = lapply(taqredo.r500, make_pa),
#                      r500 = lapply(taqredo.r100, make_pa),
#                      r1000 = lapply(taqredo.r1000, make_pa))
# lapply(taqredo.list, function(x) lapply(x, dim))
# 
# 
# #----------------------------------------------------------#
# # Load OTU taxon tables
# 
# tab.taxa <- read_tidy_taxon(table.names = names(tab.t100))
# taqredo.taxa <- read_tidy_taxon(table.names = "2017ITS") # just the taq redos
# 
# # trim the tidy taxon files -- OTUs that are not longer included in the OTU table should not appear in the taxon table
# tab.taxa <- trim_tab.taxa(tab.otus = tab.t100, tab.taxa = tab.taxa)
# taqredo.taxa <- trim_tab.taxa(tab.otus = taqredo.t100, tab.taxa = taqredo.taxa)
# 
# #----------------------------------------------------------#
# # Compile community datasets by adding in OTU taxon info and study metadata
# 
# ## Sample meta data
# logTrt <- load_logTrt()
# plotLoc <- load_plotLoc()
# samps <- make_miSeq_sampleList()
# samps %>%
#   filter(sampleType == "rotplot" & is.na(bag)) %>%
#   left_join(logTrt) %>%
#   left_join(plotLoc) -> samp.meta
# samp.taqredos <- make_taqRedo_sampleList()
# 
# # ## Compile data subsets
# 
# # # by sequencing run
# otuPreplist <- lapply(tab.list, function(x){
#   compile_datasubsets_wgroups(tab.otu = x, samp.data = samp.meta, tab.taxa = tab.taxa)
# })
# names(otuPreplist[[1]])
# lapply(otuPreplist[[1]], function(x) dim(x$otu))
# #
# # # by overlap
# overlap.otuPreplist <- lapply(tab.list, function(x){
#   make_overlap_subsets(tab.otu = x, samp.meta = samp.meta, add_group_tf = T, tab.taxa = tab.taxa)
# })
# names(overlap.otuPreplist[[1]])
# lapply(overlap.otuPreplist[[1]], function(x) dim(x$otu))
# #
# # # by deploy
# deploy.otuPreplist <- lapply(tab.list, function(x){
#   make_deploy_subsets(tab.otu = x, samp.meta = samp.meta, add_group = T, tab.taxa = tab.taxa)
# })
# names(deploy.otuPreplist[[1]])
# lapply(deploy.otuPreplist[[1]], function(x) dim(x$otu))
# #
# # # taq redos
# taqredos.otuPreplist <- lapply(taqredo.list, function(x){
#   compile_datasubsets(tab.otu = x, samp.data = samp.taqredos)
# })
# names(taqredos.otuPreplist[[1]])
# 
# # # write community datasets
# saveRDS(otuPreplist, file = paste0(intermed.path, "otuPreplist.RData"))
# saveRDS(overlap.otuPreplist, file = paste0(intermed.path, "overlap_otuPreplist.RData"))
# saveRDS(deploy.otuPreplist, file = paste0(intermed.path, "deploy_otuPreplist.RData"))
# saveRDS(taqredos.otuPreplist, file = paste0(intermed.path, "taqredos_otuPreplist.RData"))

# read community datasets
otuPreplist <- readRDS(file = paste0(intermed.path, "otuPreplist.RData"))
otuPreplist$c100$`2012ITS`$meta

overlap.otuPreplist <- readRDS(file = paste0(intermed.path, "overlap_otuPreplist.RData"))
deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))
taqredos.otuPreplist <- readRDS(file = paste0(intermed.path, "taqredos_otuPreplist.RData"))



