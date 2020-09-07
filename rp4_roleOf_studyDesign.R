# rp4_roleOf_studyDesign.R
# generates ro, rd, rd.ag

# Rotplot study design as a predictor of microbial composition

# -------------------------------------------------------------------#
## Overlap subsets
### *Hyp (overlap)* Microclimate, wood species identity, and decay time predict microbial composition

# --------------------#
# NOTE: updated otuprep_dbrdaTest() to accomadate CLR matrix -- be careful order of c100 in list matters!
# db-RDA(composition ~ logLoc + topo + species4 + yrdeploy)
# ro <- otuprep_dbrdaTest(data.otuPreplist = overlap.otuPreplist,
#                         datasetType = "overlap") # this takes a while
#saveRDS(ro, file = "data_Rsynth/rotplot_intermediates/dbrda_o_clr.RData")
ro <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_o_clr.RData")
overlap.otuPreplist <- readRDS(file = paste0(intermed.path, "overlap_otuPreplist.RData"))

# --------------------#
# Make variance partitioning plots
result <- plot_varpartbars_compareOTUprep(dbrda.intermed = ro, subsettype = "o", trait_tf = F)
result$p.r2
result$p.dist
result$p.r2.raw
result$p.dist.raw

ggsave(result$p.r2,
       file = paste0(output.path, "varpart_r2_o_allrare.pdf"), 
       width = 7, height = 4)

ggsave(result$p.dist, 
       file = paste0(output.path, "varpart_dist_o_allrare.pdf"), 
       width = 7, height = 6)

ggsave(result$p.r2.raw,
       file = paste0(output.path, "varpart_r2_o.pdf"),
       width = 4, height = 3)

ggsave(result$p.dist.raw,
       file = paste0(output.path, "varpart_dist_o.pdf"),
       width = 4, height = 3)

# --------------------#
# Make ordinations
#overlap.t100pa <- overlap.otuPreplist[['t100']]
overlap.c100 <- overlap.otuPreplist[['c100']]

plot.data.overlap <- lapply(overlap.c100, 
                            function(x){makeDF_ordplot_overlap(data = x, presAbs = FALSE)})
ordplot_overlap(plot.data.overlap, output.path, no_groups = T, 
                xlims = c(-3.5, 7), ylims = c(-4, 3.5)) #produces 1 plot and 1 legend

# ordplot_overlap(plot.data.overlap, output.path, 
#                 xlims = c(-3.75, 7), ylims = c(-4, 3.5),
#                 no_groups = F) #produces 1 plot and 1 legend


# -------------------------------------------------------------------#
## Deployment subsets
### *Hyp (deploy)* Microclimate and wood species identity predict microbial composition

# --------------------#
# db-RDA(composition ~ logLoc + topo + species4)
# rd <- otuprep_dbrdaTest(data.otuPreplist = deploy.otuPreplist,
#                         datasetType = "deploy") # this takes a while
#saveRDS(rd, file = "data_Rsynth/rotplot_intermediates/dbrda_d.RData")
rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d.RData")

# --------------------#
# Make variance partitioning plots
result <- plot_varpartbars_compareOTUprep(dbrda.intermed = rd, 
                                          subsettype = "d", 
                                          trait_tf = F)

result$p.r2
result$p.dist
result$p.r2.raw
result$p.dist.raw

ggsave(result$p.r2,
       file = paste0(output.path, "varpart_r2_d_allrare.pdf"), 
       width = 7, height = 4)

ggsave(result$p.dist, 
       file = paste0(output.path, "varpart_dist_d_allrare.pdf"), 
       width = 7, height = 6)

ggsave(result$p.r2.raw,
       file = paste0(output.path, "varpart_r2_d.pdf"),
       width = 4, height = 5)

ggsave(result$p.dist.raw,
       file = paste0(output.path, "varpart_dist_d.pdf"),
       width = 4, height = 5)


# --------------------#
# # Just check out varpart for c100
# varpartTest_presAbs_deploy <- function(data){
#   
#   mat.otu.pres <- data[["otu"]]
#   envVars <- data[["meta"]]
#   
#   # remove logLoc == "mush"
#   condition <- sum(grepl("mush", envVars$logLoc))
#   if(condition !=0){
#     envVars %>%
#       filter(logLoc == "mush") -> remove
#     mat.otu.pres <- mat.otu.pres[!row.names(mat.otu.pres) %in% remove$tidySeq_sampleName,]
#     envVars <- envVars[!envVars$tidySeq_sampleName %in% remove$tidySeq_sampleName,]
#   }
#   
#   endo.var <- varpart(vegdist(mat.otu.pres, distance = 'jaccard', binary = TRUE),
#                       ~ species4,
#                       ~ topo,
#                       ~ logLoc, data = envVars)
#   
#   return(endo.var)
#   
# }
# 
# tmp <- varpartTest_presAbs_deploy(data = deploy.otuPreplist[['c100']][['2012ITS_d11']])
# tmp
# 
# tmp <- varpartTest_presAbs_deploy(data = deploy.otuPreplist[['c100']][['201216S_d11']])
# tmp

# --------------------#
# Make ordinations

#logLoc - composition plot
#data.deploy.t100pa <- deploy.otuPreplist[['t100']]
data.deploy.c100 <- deploy.otuPreplist[['c100']]

plot.data.deploy.logLoc <- lapply(data.deploy.c100, function(x){
  makeDF_ordplot_deploy(data = x, presAbs = FALSE, GRPname = "logLoc")
})

ordplot_deploy(plot.data.deploy.x = plot.data.deploy.logLoc, 
               output.path = output.path, 
               GRPname = "logLoc",
               no_groups = TRUE,
               xlims = c(-2,2), ylims = c(-2,2)) #produces 1 plot

#topo - composition plot
plot.data.deploy.topo <- lapply(data.deploy.c100, function(x){
  makeDF_ordplot_deploy(data = x, presAbs = FALSE, GRPname = "topo")
})
ordplot_deploy(plot.data.deploy.x = plot.data.deploy.topo, 
               output.path = output.path, 
               GRPname = "topo",
               no_groups = TRUE,
               xlims = c(-2,2), ylims = c(-2,2)) #produces 1 plot


# -------------------------------------------------------------------#
### *sub hyp (deploy)* The influence of wood species identity is largely underpinned by differences between angiosperm and gymnosperm taxa

# --------------------#
# db-RDA(composition ~ logLoc + topo + angio.gymno)
# rd.ag <- lapply(data.deploy.c100, function(x){
#  dbrdaTest_ag_clr(data = x, datasetType = "deploy")
#  })
#saveRDS(rd.ag, file = "data_Rsynth/rotplot_intermediates/dbrda_ag_d_clr.RData")
rd.ag <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_ag_d_clr.RData")

anova.df <- list_to_df(lapply(rd.ag, function(x){x$anova.df}))
anova.df
var.fract <- list_to_df(lapply(rd.ag, function(x){x$var.fract}))
var.fract
dist.df <- list_to_df(lapply(rd.ag, function(x){x$dist.df}))
dist.df

# --------------------#
# Make ordinations

#angiosperm/gymnosperm - composition plot
plot.data.deploy.ag <- lapply(data.deploy.c100, function(x){
  makeDF_ordplot_deploy(data = x, presAbs = FALSE, GRPname = "angio.gymno")
})
ordplot_deploy(plot.data.deploy.x = plot.data.deploy.ag, 
               output.path = output.path, 
               GRPname = "angiogymno",
               no_groups = TRUE,
               xlims = c(-2,2), ylims = c(-2,2)) #produces 1 plot


# --------------------#
# Summarize the % of OTUs shared between species, within species for within stem and within watershed position

# make a function to summarize the number of shared OTUs within watershed and within stem
shared.otus <- function(mat.otu, meta, SPECIES, YRDEPLOY){
  
  # divide the OTU tables into treatment levels
  # create composite treatment levels
  meta %>%
    mutate(species_loglocTopo_yrdeploy = paste(species, logLoc, topo, yrdeploy, sep = "_")) -> meta.tmp
  
  # make a sample list for each treatment level
  LEVEL <- unique(meta.tmp$species_loglocTopo_yrdeploy)
  test <- list()
  for(i in 1:length(LEVEL)){
    
    #pull samples
    meta.tmp %>%
      filter(species_loglocTopo_yrdeploy == LEVEL[i]) -> selection.df
    selection.df
    selection <- selection.df$tidySeq_sampleName
    #summarize pres/abs
    mat.otu.select <- mat.otu[row.names(mat.otu) %in% selection,]
    if(is.null(dim(mat.otu.select))){
      presabs.vec <- mat.otu.select
    }else{
      presabs.vec <- colSums(mat.otu.select)
    }
    test[[i]]<- presabs.vec
  }
  
  
  df <- data.frame(sapply(test, rbind))
  colnames(df) <- LEVEL
  row.names(df) <- colnames(mat.otu)
  df.pa <- (df != 0)*1 # turn into presence/absence
  mat.trt <- t(df.pa)
  
  # just look at 1 species
  sel <- grepl(SPECIES, row.names(mat.trt)) & grepl(YRDEPLOY, row.names(mat.trt))
  tmp.sel <- row.names(mat.trt)[sel]
  tmp.mat.trt<- mat.trt[row.names(mat.trt) %in% tmp.sel,]
  
  #filter out OTUs that are not in this species+yrdeploy
  tmp.mat.trt <- tmp.mat.trt[,colSums(tmp.mat.trt) != 0]
  total.shared.otus <- dim(tmp.mat.trt)[2]
  
  # of the 641 species shared among CEOC deployed in 2009 and sampled in 2012...
  #1. how many otus are shared in b vs t 
  # in L
  tmp <- tmp.mat.trt[grepl("L", row.names(tmp.mat.trt)),]
  b.vs.t_Lshared <- sum(colSums(tmp) == 2) #136 of 641
  # in H
  tmp <- tmp.mat.trt[grepl("H", row.names(tmp.mat.trt)),]
  b.vs.t_Hshared <- sum(colSums(tmp) == 2) #112 of 641
  
  # 2. how many otus are shared in L vs H
  # in b
  tmp <- tmp.mat.trt[grepl("b", row.names(tmp.mat.trt)),]
  L.vs.H_bshared <- sum(colSums(tmp) == 2) #110 of 641
  
  # in t
  tmp <- tmp.mat.trt[grepl("t", row.names(tmp.mat.trt)),]
  L.vs.H_tshared <- sum(colSums(tmp) == 2) #99 of 641
  
  # summarize in a df
  summ <- data.frame(total.shared.otus, b.vs.t_Lshared, b.vs.t_Hshared, 
                     L.vs.H_bshared, L.vs.H_tshared)
  summ
  return(summ)
}

names(overlap.t100pa)
names(overlap.t100pa$`201216S_o`)

shared.otus1 <- shared.otus(mat.otu = overlap.t100pa$`2012ITS_o`$otu, 
                            meta = overlap.t100pa$`2012ITS_o`$meta, 
                            SPECIES = "QUVE", 
                            YRDEPLOY = 2009)
shared.otus1

shared.otus2 <- shared.otus(mat.otu = overlap.t100pa$`201416S_o`$otu, 
                            meta = overlap.t100pa$`201416S_o`$meta, 
                            SPECIES = "CEOC", 
                            YRDEPLOY = 2009)
shared.otus2

df<- rbind(shared.otus1, shared.otus2)
df$comm <- c("fungi","bacteria")
df$b.vs.t_Lshared/df$total.shared.otus
df$b.vs.t_Hshared/df$total.shared.otus
df$L.vs.H_bshared/df$total.shared.otus
df$L.vs.H_tshared/df$total.shared.otus

