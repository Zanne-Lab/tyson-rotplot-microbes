# rp5_roleOf_woodTraits.R

# Wood traits as predictors of microbial composition

# -------------------------------------------------------------------#
# First, check out how wood traits are correlated
wood.traits <- load_fungalEnvirTraits()
result <- plot_woodTraits(woodTraits = wood.traits)

pdf(file = paste0(output.path, "woodtraits_cor_09.pdf"),
    width = 7, height = 7)
#corrplot(result$cor.09, type = "upper", diag = F)
corrplot(result$cor.09, type = "upper", diag = F, method = "number")
dev.off()

pdf(file = paste0(output.path, "woodtraits_cor_11.pdf"), 
    width = 7, height = 7)
# #corrplot(result$cor.11, type = "upper", diag = F)
corrplot(result$cor.11, type = "upper", diag = F, method = "number")
dev.off()

# conduit diameter and conduit length
# conduit diameter and Ca
# parenchymaFrac and conduit length
# parenchymaFrac and density
# parenchymaFrac and lignin

# get rid of conduit dimater and parenchymaFrac for models?



# -------------------------------------------------------------------#
#Next, check out how wood species fall out in multivariate trait space

ggsave(result$p.09,
       file = paste0(output.path, "woodtraits_ord_09.pdf"),
       width = 7, height = 4)

ggsave(result$p.11,
       file = paste0(output.path, "woodtraits_ord_11.pdf"),
       width = 7, height = 4)



# -------------------------------------------------------------------#
### *Hyp (deploy)* Some wood traits are better than others for explaining microbial composition
# Only investigating this with "deployment subsets" because it covers more wood species and comprises a wide range of trait space unlike the "overlap subsets" with only 3 wood species.

# --------------------#
# ordiR2step composition ~ wood traits (deploy)

#this takes forever, so I ran it on the cluster
#rd <- otuprep_ordiR2step_woodtraits(data.otuPreplist = deploy.otuPreplist, 
#                                    woodTraits = wood.traits)
#saveRDS(rd, file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d.RData")
#rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d.RData")
#woodTraits <- readRDS(file = "data_Rsynth/rotplot_intermediates/woodTraits.RData")

# logTrt <- load_logTrt()
# plotLoc <- load_plotLoc()
# samps <- make_miSeq_sampleList()
# samps %>%
#   filter(sampleType == "rotplot" & is.na(bag)) %>%
#   left_join(logTrt) %>%
#   left_join(plotLoc) -> samp.meta

# add c100 data
# deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))
# data <- deploy.otuPreplist[['c100']]
# # do ordistep
# tmp.list <- lapply(data, function(x){
#     tmp <- ordiR2step_woodtraits_clr(data = x, woodTraits = woodTraits)
#     return(tmp)
# })
# tmp.list
# tmp.df <- list_to_df(lapply(tmp.list, function(x){x$anova}))
# tmp.df %>%
#     rename('subset' = 'source') -> anova.df
# # generate index for annotating the data subset
# df.indx <- unique(rd[,c("yrharv","gene","yrdeploy","yearSubset")])
# df.indx %>%
#   mutate(subset = paste0(yrharv, gene, "_", yrdeploy)) -> df.indx
# anova.df %>%
#   left_join(df.indx) %>%
#   filter(!is.na(yrharv)) %>%
#   dplyr::rename('SumOfSqs'='Variance') %>%
#   select(-subset) %>%
#   mutate(otuprep = "c100") -> anova.df
# 
# # add c100 data to rd
# rd.c100 <- rbind(rd, anova.df)
# saveRDS(rd.c100, file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d_clr.RData")
rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d_clr.RData")

# --------------------#
#Do trait model selection for just the c100 OTU table -- but scale the trait matrix first
# deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))
# data <- deploy.otuPreplist$c100

# get rid of conduit dimater and parenchymaFrac for models
woodTraits <- readRDS(file = "data_Rsynth/rotplot_intermediates/woodTraits.RData")
woodTraits %>%
  select(-c(parenchymaFrac,conduitDiam)) -> woodTraits
woodTraits.s <- scale(woodTraits[,-1])
woodTraits.s <- cbind(woodTraits[,1], woodTraits.s)

# do ordistep
# tmp.list.s <- lapply(data, function(x){
#   tmp <- ordiR2step_woodtraits_clr(data = x, woodTraits = woodTraits.s)
#   return(tmp)
# })
# tmp.list.s
#saveRDS(tmp.list.s, file = "data_Rsynth/rotplot_intermediates/traitmodelsel_tmplist_s_clr.RData")
# tmp.list.s <- readRDS(file = "data_Rsynth/rotplot_intermediates/traitmodelsel_tmplist_s_clr.RData")
# 
# tmp.df <- list_to_df(lapply(tmp.list.s, function(x){x$anova}))
# tmp.df %>%
#   rename('subset' = 'source') -> anova.list

# --------------------#
# Make a grid for traits in the best model

# extract signif terms
# anova.list %>%
#   filter(term != "Residual") %>%
#   select(term, Pr..F., subset) -> rd.terms
# 
# # annotate with model fit
# rd.fit <- unique(anova.list[,c("constr.chi","tot.chi","subset")])
# rd.terms %>%
#   left_join(rd.fit) %>%
#   mutate(percConst = (constr.chi / tot.chi) * 100 ) -> rd.terms
# 
# # update subset info
# rd.terms %>%
#   separate(subset, into = c("yrharvGene","deploy","commsub"), sep = "_", remove = F) %>%
#   separate(yrharvGene, into = c("yrharv","gene"), sep = 4) %>%
#   mutate(yearSubset = paste(yrharv, deploy, sep = "_")) -> rd.terms
# 
# #factor order
# levels.tmp <- c("2012_d11", "2014_d11","2012_d09","2014_d09")
# labels.tmp <- c("2011-12 (1 yr)","2011-14 (3 yrs)", "2009-12 (3 yrs)", "2009-14 (5 yrs)")
# rd.terms$yearSubset <- factor(rd.terms$yearSubset, 
#                               level = levels.tmp, 
#                               labels = labels.tmp)
# rd.terms$gene <- factor(rd.terms$gene, level = c("ITS","16S"), labels = c("Fungi","Bacteria"))
# params <- plottingParams_rotplotNobags()
# rd.terms$term <- factor(rd.terms$term, level = rev(params$traits$levels), labels = rev(params$traits$labels))
# 
# # colors
# light.red <- brewer.pal(n=8, name = "RdBu")[2] # 2011-12 (1yr)
# dark.red <- brewer.pal(n=8, name = "RdBu")[1] # 2011-2014 (3yrs)
# light.blue <- brewer.pal(n=8, name = "RdBu")[7] # 2009-2012 (3yrs)
# dark.blue <- brewer.pal(n=8, name = "RdBu")[8] # 2009-2014 (5yrs)
# deploy.colors <- c(light.red, dark.red, light.blue, dark.blue)
# 
# #remove commsub
# rd.terms %>%
#   filter(is.na(commsub)) -> rd.terms

#plot
# p <- ggplot(rd.terms, aes(x = yearSubset, y = term, fill = yearSubset))+
#   geom_tile(color = "black") +
#   facet_grid(~gene) +
#   xlab(" ") + ylab("Wood trait") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_manual(name = "Deployment subset",
#                     values = deploy.colors)
# p
# ggsave(p, filename = paste0(output.path, "woodTraits_bestmodGrid.pdf"), 
#        width = 5, height = 5)


# --------------------#
# Summarize the variance explained by each trait in the final models

woodTraits.s %>%
  rename("Conduit length"="conduitLeng") %>%
  rename("Density"="density") %>%
  rename("Crude protein"="crudeProt") %>%
  rename("Lignin"="lignin") %>%
  rename("Cellulose"="cellulose") %>%
  rename("Hemicellulose"="hemicellulose") -> woodTrait.s.names

# tmp.list.s <- readRDS(file = "data_Rsynth/rotplot_intermediates/traitmodelsel_tmplist_s.RData")
# 
# rhs.list <- list()
# terms.list <- list()
# SUBSET <- unique(rd.terms$subset)
# for(i in 1:length(SUBSET)){
#   rd.terms %>%
#     filter(subset == SUBSET[i]) %>%
#     filter(term != "Residual") -> curr
#   curr
#   terms <- c(as.character(curr$term))
#   rhs.list[[i]] <- paste(terms, collapse = " + ")
#   terms.list[[i]]<- terms
# }
# 
# i<-1
# mat.otu.pres <- data[[i]][["otu"]]
# envVars <- data[[i]][["meta"]]
# terms <- terms.list[[i]]
# rhs <- rhs.list[[i]]
# 
# # trait data for species4 == PLOC is missing (only in d09), so exclude it
# envVars %>%
#   filter(species4 != "PLOC") -> envVars
# mat.otu.pres <- mat.otu.pres[row.names(mat.otu.pres) %in% envVars$tidySeq_sampleName,]
# 
# # set up the model
# envVars %>%
#   left_join(woodTrait.s.names) %>%
#   select(terms) -> envVars.to
# 
# cap.env <-  capscale(mat.otu.pres ~ ., data = envVars.to, distance = 'jaccard', binary = TRUE)
# 
# #vector info
# cap.env$CCA$envcentre
# 
# cca.biplot[[i]] <- data.frame(terms = row.names(cap.env$CCA$biplot), cap.env$CCA$biplot)
# 
# plot(cap.env)


# --------------------#
# Make heat map for best models across traits and otupreps

# extract signif terms
rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d_clr.RData")
rd %>%
  filter(term != "Residual") %>%
  select(term, Pr..F., yearSubset, gene, otuprep) -> rd.terms

# annotate with model fit
rd.fit <- unique(rd[,c("constr.chi","tot.chi","gene","yearSubset","otuprep")])
rd.terms %>%
  left_join(rd.fit) %>%
  mutate(percConst = (constr.chi / tot.chi) * 100 ) -> rd.terms

#factor order
params <- plottingParams_rotplotNobags()
mytheme <- make_ggplot_theme()
levels.tmp <- c("2012_d11", "2014_d11","2012_d09","2014_d09")
labels.tmp <- c("2011-12 (1 yr)","2011-14 (3 yrs)", "2009-12 (3 yrs)", "2009-14 (5 yrs)")
rd.terms$yearSubset <- factor(rd.terms$yearSubset, 
                              level = levels.tmp, 
                              labels = labels.tmp)
rd.terms$gene <- factor(rd.terms$gene, level = c("ITS","16S"), labels = c("Fungi","Bacteria"))
rd.terms$term <- factor(rd.terms$term, level = rev(params$traits$levels), 
                        labels = rev(params$traits$labels))
unique(rd.terms$otuprep)
rd.terms %>%
  filter(otuprep != "t100pa") -> rd.terms

new.levels <- c("c100", paste0(params$otuprep$levels[-1], "pa"))
new.levels
rd.terms$otuprep <- factor(rd.terms$otuprep, 
                           level = new.levels, 
                           labels = params$otuprep$labels)

p <- ggplot(rd.terms, aes(x = otuprep, y = term, fill = percConst))+
  geom_tile(color = "black") +
  facet_grid(gene~yearSubset) +
  xlab("Rarefaction level") + ylab("Wood trait") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_distiller(name = "% Constrained", direction = 1) +
  mytheme
p

ggsave(p, 
       filename = paste0(output.path, "woodTraits_composition_ordistep_clr.pdf"), 
       width = 6.5, height = 5)

# black tiles
p <- ggplot(rd.terms, aes(x = otuprep, y = term))+
  geom_tile(color = "black", fill = "gray") +
  facet_grid(gene~yearSubset) +
  xlab("Rarefaction level") + ylab("Wood trait") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_distiller(name = "% Constrained", direction = 1)
p
ggsave(p, 
       filename = paste0(output.path, "woodTraits_composition_ordistep_blacktiles_clr.pdf"), 
       width = 5, height = 5)


# remove rarefaction level
rd.terms %>%
  filter(otuprep == "CLR")-> tmp.rd.terms

p <- ggplot(tmp.rd.terms, aes(x = yearSubset, y = term, 
                              fill = percConst))+
  geom_tile(color = "black") +
  facet_grid(~gene) +
  xlab(" ") + ylab("Wood trait") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_distiller(name = "% Constrained", direction = 1)
p
ggsave(p, 
       filename = paste0(output.path, "woodTraits_composition_ordistep.pdf"), 
       width = 7, height = 5)


# -------------------------------------------------------------------#
### *Hyp (deploy)* Wood traits account for a majority of the microbial community variation that is explained by wood species identity

# --------------------#
# terms that are included across all otupreps
rd %>%
  filter(term != "Residual") %>%
  select(term, Pr..F., yrharv, gene, yrdeploy, otuprep) %>%
  group_by(yrharv, gene, yrdeploy, term) %>%
  summarize(n = length(Pr..F.)) %>%
  filter(n == 4) %>%
  mutate(yearGene = paste(yrharv, gene, sep ="")) %>%
  mutate(table = paste(yearGene, yrdeploy, sep = "_")) -> rd.all


# --------------------#
#Plot and print contrained ordinations using the best set of traits and non-rarefied data -- skip this

# # plot the constrained ordinations for otuprep = t100pa
# dbrdaTest_woodTraits_finalmod_plots(rd.all = rd.all, 
#                                     woodTraits = wood.traits, 
#                                     data.deploy = deploy.otuPreplist[['t100']])
# v.list <- dbrdaTest_woodTraits_finalmod_v(rd.all = rd.all, 
#                                           woodTraits = wood.traits, 
#                                           data.deploy = deploy.otuPreplist[['t100']])
# 
# # percent constrained
# rd.all %>%
#   select(term, yrharv, gene, yrdeploy) %>%
#   left_join(rd) %>%
#   filter(otuprep == "t100pa") -> tmp
# rd.all.mod <- unique(tmp[,c("yrharv","gene","yrdeploy","constr.chi","unconst.chi","tot.chi")])
# rd.all.mod %>%
#   mutate(perc.const = round(constr.chi / tot.chi, digits =2)) %>%
#   select(gene, yrharv, yrdeploy, perc.const) %>%
#   spread(key = gene, value = perc.const)
# 
# rd.all.mod
# write.csv(rd.all.mod, file = paste0(output.path, "dbrda_woodtraits_d.csv"))
# 
# 
# # extract OTU associations
# v.df <- extract_keyOTUs(v.list = v.list, cap.choice = "CAP1", table.choice = "2012ITS_d11", tab.taxa = tab.taxa)
# quant <- quantile(v.df$CAP1, c(.03, .97))
# v.df %>%
#   filter(CAP1 < quant[1]) -> v.df.low
# v.df %>%
#   filter(CAP1 > quant[2]) -> v.df.high
# #View(v.df.high)

# --------------------#
# Test db-RDAs using the best traits for each otuprep
rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/ordiR2step_woodtraits_d_clr.RData")
woodTraits <- readRDS(file = "data_Rsynth/rotplot_intermediates/woodTraits.RData")
deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))

# db-RDA(composition ~ logLoc + topo + traits)
# rd.traits <- otuprep_dbrdaTest_woodTraits_clr(rd = rd,
#                                           woodTraits = woodTraits,
#                                           data.otuPreplist = deploy.otuPreplist) # this takes a while
# saveRDS(rd.traits, file = "data_Rsynth/rotplot_intermediates/dbrda_d_traits_clr.RData")
rd.traits <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d_traits_clr.RData")

# --------------------#
# Make variance partitioning plots

rd.traits <- preprocess_varpart_rd.traits(rd.traits = rd.traits)
result <- plot_varpartbars_compareOTUprep(dbrda.intermed = rd.traits, 
                                          subsettype = "d", trait_tf = T)
result$p.r2
result$p.r2.raw


# --------------------#
# Combine trait r2 values into dataframe with species r2s for Fig 2

# # overlap
# ro <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_o.RData")
# 
# # deployment
# rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d.RData")
# # use traits instead of wood species
# rd.traits <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d_traits.RData")
# rd.traits <- preprocess_varpart_rd.traits(rd.traits = rd.traits)
# 
# # filter for just fungi and bacteria, 2012 and 2014, otuprep = t100 or raw
# ro$varpart.df %>%
#   filter(is.na(group)) %>%
#   filter(otuprep == "t100") %>%
#   mutate(model = "ro") %>%
#   select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> ro.varp
# rd$varpart.df %>%
#   filter(!grepl('sapro', subset)) %>%
#   filter(!grepl('basidio', subset)) %>%
#   filter(otuprep == "t100") %>%
#   mutate(model = "rd") %>%
#   select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> rd.varp
# rd.traits$varpart.df %>%
#   filter(otuprep == "t100") %>%
#   mutate(model = "rd.traits") %>%
#   mutate(pval = NA) %>%
#   select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> rd.traits.varp
# varp.df <- rbind(ro.varp, rd.varp, rd.traits.varp)
# 
# # check that Adj.R.square is (approx) the same for topo and logLoc between "rd" and "rd.traits" models -- these are "indfract".. individual fraction attributable to specific term
# varp.df %>%
#   filter(model %in% c("rd","rd.traits")) %>%
#   filter(terms  %in% c("topo","logLoc")) %>%
#   arrange(subset, terms)

