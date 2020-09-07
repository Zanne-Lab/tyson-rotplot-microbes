# rp6_roleOf_coclade.R

# Co-occurring microbial clade as predictors of microbial composition

# -------------------------------------------------------------------#
### *Hyp (deploy)* Fungal and bacterial composition will be highly correlated because (i) fungi can house bacteria, (ii) both clades engineer common habitat, (iii) both clades use common resources
# Only investigating this hypothesis with the "deployment subsets"
overlap.otuPreplist <- readRDS(file = paste0(intermed.path, "overlap_otuPreplist.RData"))
deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))


# --------------------#
# Run mantel tests
mantel.list <- list()
for(i in 1:length(deploy.otuPreplist)){
  mantel.list[[i]] <- batch_write_manteltests_rotplots(data.overlap.choosePrep = overlap.otuPreplist[[i]],
                                                       data.deploy.choosePrep = deploy.otuPreplist[[i]], 
                                                       presAbs = FALSE)
}
names(mantel.list) <- names(deploy.otuPreplist)
mantel.df <- list_to_df(mantel.list)

mantel.df
write.csv(mantel.df, 
          file = paste0(output.path, "mantel_od_clr.csv"))


# params <- plottingParams_rotplotNobags()
# params$yearSubset$levels
# 
# mantel.df$yearSubset <- factor(mantel.df$yearSubset, 
#                                levels = c("2012_o","2014_o", 
#                                           params$yearSubset$levels))
# mantel.df$source <- factor(mantel.df$source, 
#                            levels = c("r100pa","r500pa","r1000pa","t100pa"))
# ggplot(mantel.df, aes(x = source, y = r, fill = yearSubset)) +
#   geom_bar(stat= "identity", position = "dodge")



# -------------------------------------------------------------------#
### *Hyp 1 (deploy)* Fungal composition explains bacterial composition because (i) fungi can house bacteria, (ii) fungi engineer wood habitat, (iii) overlapping resource use
### *Hyp 2 (deploy)* Bacterial composition explains fungal composition because (i) bacteria are critical for importing N, a limiting nutrient, (ii) bacteria engineer wood habitat to a limited extent, (iii) overlapping resource use
### *Hyp 3 (deploy)* Fungal composition explains bacterial composition to a greater extent than the reverse because (i) fungi engineer wood habitat to a greater extent, (ii) fungi have a bigger body size

# first, make a dataframe for all the models -- don't include sapros and basidios
data <- deploy.otuPreplist$c100
tab.names <- names(data)
condition <- !grepl("sapro", tab.names) & !grepl("basidio", tab.names)
tab.names <- tab.names[condition]
pair.df <- create_fung_bact_pairs_rotplot(tab.names)
pair.df
curr.data <- data[tab.names]

# second, determine which axes to keep for fungal and bacterial compositions
# this takes forever, it needs to happen on the cluster
# for each unique commX in pair.df....
#modelSelectionMDS_jaccard(commx.otu)
# make a list of the signif axes to keep -- for now, this is just the first 2 axes
commX.signifaxes.list <- rep(list(c(1,2)), length(pair.df$commX.tab))
names(commX.signifaxes.list) <- pair.df$commX.tab

# --------------------#
# third, fit a dbrda with 2 predictor matrices: select composition axes of co-occurring clade and study factors
# names(deploy.otuPreplist)
# presAbs.vec <- c(FALSE, TRUE, TRUE, TRUE)
# r.list <- list()
# for(i in 1:length(deploy.otuPreplist)){
#   r.list[[i]] <- batch_varpart_commx_rotplot(pair.df, 
#                                              commX.signifaxes.list,
#                                              datasets = deploy.otuPreplist[[i]],
#                                              presAbs = presAbs.vec[i])
#   print(i)
# }
# names(r.list) <- names(deploy.otuPreplist)
# saveRDS(r.list, file = "data_Rsynth/rotplot_intermediates/dbrda_commx_d_clr.RData")
r.list <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_commx_d_clr.RData")
r.df <- list_to_df(r.list)

# --------------------#
# Make variance partitioning plots

params <- plottingParams_rotplotNobags()
color.vec <- params$topo$colors
names(color.vec) <- NULL

# isolate just the infract parts and tidy factors
r.df %>%
  filter(fractType == "ind") %>%
  rename('otuprep' = 'source') %>%
  mutate(model = "rd.comm") %>%
  rename('gene'='commY.name') %>%
  separate(yearSubset, into = c("yrharv", "yrdep"), remove = F) %>%
  mutate(subset = factor(yearSubset,
                         levels = rev(params$yearSubset$levels.new), 
                         labels = rev(params$yearSubset$labels.new))) %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.squared, pval) %>%
  mutate(terms = recode_factor(terms, 
                               commx = "Fungal-bacterial\ninteractions",
                               env = "Shared environ.")) %>%
  mutate(gene = recode_factor(gene, `ITS`="Fungi",`16S`="Bacteria", .ordered = T)) %>%
  mutate(otuprep = recode_factor(otuprep, 
                                 c100 = "CLR",
                                 r100pa = "r100",
                                 r500pa = "r500",
                                 r1000pa = "r1000", .ordered = T)) -> rd.comm.varp

# ggplot
mytheme <- make_ggplot_theme()
p <- ggplot(data = rd.comm.varp, aes(x = subset, 
                                     y = Adj.R.squared, 
                                     fill = terms)) +
  facet_grid(otuprep~gene) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  scale_fill_manual(name = "Predictor type", values = color.vec) +
  xlab("Deployment subset") + ylab("Adjusted R squared") 
p
ggsave(p, filename = paste0(output.path, "dbrda_commx_d_allrare.pdf"), 
       width = 5, height = 5)

####
#only raw data 

rd.comm.varp %>%
  filter(otuprep == "c100") %>%
  mutate(subset = factor(subset, 
                         levels = rev(levels(subset)))) -> rd.comm.varp.raw

p <- ggplot(data = rd.comm.varp.raw, aes(x = terms, 
                                         y = Adj.R.squared)) +
  facet_grid(subset~gene) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  xlab("Predictor type") + ylab("Adjusted R squared")
p
ggsave(p, filename = paste0(output.path, "dbrda_commx_d.pdf"), 
       width = 5, height = 5)

