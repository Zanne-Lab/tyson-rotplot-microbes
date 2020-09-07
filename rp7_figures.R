# rp7_figures.R

# objects needed for each part...

#-----------#
# make r.df
r.list <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_commx_d_clr.RData")
names(r.list)
r.df <- list_to_df(r.list)

#-----------#
# make varp.df
ro <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_o_clr.RData")
rd <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d.RData") # this has c100 in it
rd.traits <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_d_traits_clr.RData")
rd.traits <- preprocess_varpart_rd.traits(rd.traits = rd.traits)
ro$varpart.df %>%
  filter(is.na(group)) %>%
  filter(otuprep == "c100") %>%
  mutate(model = "ro") %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> ro.varp
rd$varpart.df %>%
  filter(!grepl('sapro', subset)) %>%
  filter(!grepl('basidio', subset)) %>%
  filter(otuprep == "c100") %>%
  mutate(model = "rd") %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> rd.varp
rd.traits$varpart.df %>%
  filter(otuprep == "c100") %>%
  mutate(model = "rd.traits") %>%
  mutate(pval = NA) %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.square, pval) -> rd.traits.varp
varp.df <- rbind(ro.varp, rd.varp, rd.traits.varp)

#-----------#
# load deploy.otuPreplist
deploy.otuPreplist <- readRDS(file = paste0(intermed.path, "deploy_otuPreplist.RData"))

#-----------#
# load rd.ag
rd.ag <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_ag_d_clr.RData")


#-----------#
# Figure 1 = none
# Figure 2 = r.df, varp.df
# Figure 3 = r.df, varp.df
# Figure 4 = deploy.otuPreplist[['c100']]
# Figure 5 = r.df
# In-text calculations = r.df, varp.df, rd.ag (rp4X.R)


# -------------------------------------------------------------------#
# Figure 1 is a conceptual figure.

# -------------------------------------------------------------------#
# Make Figure 2 -- Species, Watershed, Stem; r2 from each model for non-rarified data

r.df %>%
  filter(fractType == "ind") %>%
  rename('otuprep' = 'source') %>%
  mutate(model = "rd.comm") %>%
  rename('gene'='commY.name') %>%
  rename('subset'='yearSubset') %>%
  separate(subset, into = c("yrharv", "dropthing"), remove = F) %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.squared, pval) %>%
  mutate(terms = recode_factor(terms, 
                               env = "Environ.", commx = "Co-clade")) %>%
  rename('Adj.R.square'='Adj.R.squared') %>%
  filter(otuprep == "c100") -> rd.comm.varp.raw
varp.df1 <- rbind(rd.comm.varp.raw, varp.df)


# just ro
varp.df1 %>%
  filter(model == "ro") %>%
  mutate(gene = recode_factor(gene, ITS = "Fungi", 
                              `16S` = "Bacteria", .ordered = T)) %>%
  mutate(terms = recode_factor(terms, 
                               yrdeploy = "Year deployed",
                               logLoc = "Within stem",
                               topo = "Within watershed",
                               species4 = "Wood species")) -> varp.df1.ro

mytheme <- make_ggplot_theme()

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

p <- ggplot(varp.df1.ro, aes(x = terms, y = Adj.R.square, fill = yrharv)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  mytheme +
  facet_grid(yrharv~gene) +
  ylab("Adjusted R squared") + xlab("") +
  scale_fill_manual(values = c(light.purple, dark.purple)) + 
  guides(fill=F)
p
#no color
p <- ggplot(varp.df1.ro, aes(x = terms, y = Adj.R.square)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  mytheme +
  facet_grid(yrharv~gene) +
  ylab("Adjusted R squared") + xlab("") 
p
ggsave(p, filename = paste0(output.path, "fig2_o_clr.pdf"), 
       width = 4, height = 2.25)


# just rd
varp.df1 %>%
  filter(model != "ro") %>%
  filter(!terms %in% c("logLoc","topo")) -> tmp
varp.df1 %>%
  filter(model == "rd") %>%
  filter(terms %in% c("logLoc","topo")) -> tmp2
varp.df2 <- rbind(tmp, tmp2)

varp.df2 %>%
  mutate(subset = recode_factor(subset, 
                                `201216S_d11` = "1 yr (2011-12)",
                                `2012ITS_d11` = "1 yr (2011-12)", 
                                `2012_d11` = "1 yr (2011-12)",
                                
                                `201416S_d11` = "3 yrs (2011-14)",
                                `2014ITS_d11` = "3 yrs (2011-14)",
                                `2014_d11` = "3 yrs (2011-14)",
                                
                                `201216S_d09` = "3 yrs (2009-12)",
                                `2012ITS_d09` = "3 yrs (2009-12)",
                                `2012_d09` = "3 yrs (2009-12)",
                                
                                `201416S_d09` = "5 yrs (2009-14)",
                                `2014ITS_d09` = "5 yrs (2009-14)",
                                `2014_d09` = "5 yrs (2009-14)",
                                
                                .ordered = T)) %>%
  filter(terms %in% c("species4", "topo", "logLoc")) %>%
  mutate(terms = recode_factor(terms, 
                               logLoc = "Within stem",
                               topo = "Within watershed",
                               species4 = "Wood species")) %>%
  mutate(gene = recode_factor(gene, ITS = "Fungi", 
                              `16S` = "Bacteria", .ordered = T)) -> varp.df2

p <- ggplot(varp.df2, aes(x = terms, y = Adj.R.square, fill = subset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("") +
  scale_fill_manual(values = c(light.red, light.blue, dark.red, dark.blue)) +
  guides(fill=F)
p
#no color
p <- ggplot(varp.df2, aes(x = terms, y = Adj.R.square)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("")
p
ggsave(p, filename = paste0(output.path, "fig2_d_clr.pdf"), 
       width = 4, height = 4.25)


# -------------------------------------------------------------------#
# Make Figure 3 -- Species, Traits; r2 from each model for non-rarified data

varp.df1 %>%
  filter(model %in% c("rd","rd.traits")) %>%
  filter(terms %in% c("species4","woodtraits")) %>%
  mutate(subset = recode_factor(subset, 
                                `201216S_d11` = "1 yr (2011-12)",
                                `2012ITS_d11` = "1 yr (2011-12)", 
                                `2012_d11` = "1 yr (2011-12)",
                                
                                `201216S_d09` = "3 yrs (2009-12)",
                                `2012ITS_d09` = "3 yrs (2009-12)",
                                `2012_d09` = "3 yrs (2009-12)",
                                
                                `201416S_d11` = "3 yrs (2011-14)",
                                `2014ITS_d11` = "3 yrs (2011-14)",
                                `2014_d11` = "3 yrs (2011-14)",
                                
                                `201416S_d09` = "5 yrs (2009-14)",
                                `2014ITS_d09` = "5 yrs (2009-14)",
                                `2014_d09` = "5 yrs (2009-14)",
                                
                                .ordered = T)) %>%
  mutate(terms = recode_factor(terms, 
                               woodtraits = "Wood traits",
                               species4 = "Wood species")) %>%
  mutate(gene = recode_factor(gene, ITS = "Fungi", 
                              `16S` = "Bacteria", .ordered = T)) -> varp.df3

varp.df3$subset <- factor(varp.df3$subset, 
                          levels = c('1 yr (2011-12)',
                                     '3 yrs (2011-14)',
                                     '3 yrs (2009-12)',
                                     '5 yrs (2009-14)'))

p <- ggplot(varp.df3, aes(x = terms, y = Adj.R.square, fill = subset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("") +
  scale_fill_manual(values = c(light.red, light.blue, dark.red, dark.blue)) +
  guides(fill=F)
p

#no color
p <- ggplot(varp.df3, aes(x = terms, y = Adj.R.square)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("") 
p
ggsave(p, filename = paste0(output.path, "fig3_clr.pdf"), 
       width = 4, height = 4.5)


# -------------------------------------------------------------------#
# Make Figure 4 -- Microbial compositions differ based on stem and watershed position

data.deploy.c100 <- deploy.otuPreplist[['c100']]

#"Within stem position"
#logLoc - composition plot
plot.data.deploy.logLoc <- lapply(data.deploy.c100, function(x){
  makeDF_ordplot_deploy(data = x, presAbs = FALSE, GRPname = "logLoc")
})

ordplot_deploy(plot.data.deploy.x = plot.data.deploy.logLoc, 
               output.path = output.path, 
               GRPname = "logLoc",
               no_groups = TRUE,
               xlims = c(-2,2), ylims = c(-2,2)) #produces 1 plot

#"Within log position"
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
# Make Figure 5 -- Shared environment, Co-clade; r2 from each model for non-rarified data

r.df %>%
  filter(fractType == "ind") %>%
  rename('otuprep' = 'source') %>%
  mutate(model = "rd.comm") %>%
  rename('gene'='commY.name') %>%
  rename('subset'='yearSubset') %>%
  separate(subset, into = c("yrharv", "dropthing"), remove = F) %>%
  select(model, subset, gene, yrharv, otuprep, terms, Adj.R.squared, pval) %>%
  rename('Adj.R.square'='Adj.R.squared') %>%
  filter(otuprep == "c100") -> rd.comm.varp.raw

# just rd
rd.comm.varp.raw %>%
  mutate(subset = recode_factor(subset, 
                                `201216S_d11` = "1 yr (2011-12)",
                                `2012ITS_d11` = "1 yr (2011-12)", 
                                `2012_d11` = "1 yr (2011-12)",
                                
                                `201216S_d09` = "3 yrs (2009-12)",
                                `2012ITS_d09` = "3 yrs (2009-12)",
                                `2012_d09` = "3 yrs (2009-12)",
                                
                                `201416S_d11` = "3 yrs (2011-14)",
                                `2014ITS_d11` = "3 yrs (2011-14)",
                                `2014_d11` = "3 yrs (2011-14)",
                                
                                `201416S_d09` = "5 yrs (2009-14)",
                                `2014ITS_d09` = "5 yrs (2009-14)",
                                `2014_d09` = "5 yrs (2009-14)",
                                
                                .ordered = T)) %>%
  mutate(terms = recode_factor(terms, 
                               env = "Shared environ.", 
                               commx = "Fungal-bacterial\ninteractions")) %>%
  mutate(gene = recode_factor(gene, ITS = "Fungi", 
                              `16S` = "Bacteria", .ordered = T)) -> rd.comm.varp2

rd.comm.varp2$terms <- factor(rd.comm.varp2$terms, levels =
                                c('Fungal-bacterial\ninteractions','Shared environ.'))
rd.comm.varp2$subset <- factor(rd.comm.varp2$subset, levels = 
                                 c("1 yr (2011-12)","3 yrs (2011-14)",
                                   "3 yrs (2009-12)","5 yrs (2009-14)"))

p <- ggplot(rd.comm.varp2, aes(x = terms, y = Adj.R.square, fill = subset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("") +
  scale_fill_manual(values = c(light.red,  
                               dark.red, 
                               light.blue, dark.blue)) +
  guides(fill=F)
p
ggsave(p, filename = paste0(output.path, "fig5_color_clr.pdf"), 
       width = 4, height = 4.5)



# no color
p <- ggplot(rd.comm.varp2, aes(x = terms, y = Adj.R.square)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  mytheme +
  facet_grid(subset~gene) +
  ylab("Adjusted R squared") + xlab("") +
  scale_x_discrete(limits=c("Fungal-bacterial\ninteractions","Shared environ."))
p
ggsave(p, filename = paste0(output.path, "fig5_clr.pdf"), 
       width = 4, height = 4.5)


# -------------------------------------------------------------------#
# Do in-text calculations based on Fig2 and Fig3 results
varp.df1

#1
# The position of a fungal community within the stem and watershed each explained about X-X% of variation in community composition;
varp.df1 %>%
  filter(gene == "ITS") %>%
  filter(model %in% c("ro","rd")) %>%
  filter(otuprep == "c100") %>%
  filter(terms %in% c("logLoc","topo")) %>%
  mutate(perc.var.expl = Adj.R.square * 100) %>%
  summarize(n.min = min(perc.var.expl),
            n.max = max(perc.var.expl))

# ...whereas, wood species identity explained X-X% (Fig 2).
varp.df1 %>%
  filter(gene == "ITS") %>%
  filter(model %in% c("ro","rd")) %>%
  filter(otuprep == "c100") %>%
  filter(terms %in% c("species4")) %>%
  mutate(perc.var.expl = Adj.R.square * 100) %>%
  summarize(n.min = min(perc.var.expl),
            n.max = max(perc.var.expl))

#2
# For bacterial community composition, wood species identity also tended to explain the largest amount of variation, but this factor explained X-X% less in bacterial relative to fungal communities  (Fig 1, Fig 2)
varp.df1 %>%
  filter(gene == "16S") %>%
  filter(model %in% c("ro","rd")) %>%
  filter(otuprep == "c100") %>%
  filter(terms %in% c("species4")) %>%
  mutate(perc.var.expl = Adj.R.square * 100) -> bact.tmp
varp.df1 %>%
  filter(gene == "ITS") %>%
  filter(model %in% c("ro","rd")) %>%
  filter(otuprep == "c100") %>%
  filter(terms %in% c("species4")) %>%
  mutate(perc.var.expl = Adj.R.square * 100) -> fung.tmp
tmp.df <- data.frame(bact.percexpl = bact.tmp$perc.var.expl,
                     fung.percexpl = fung.tmp$perc.var.expl)
tmp.df %>%
  mutate(bf.perc = (bact.percexpl/fung.percexpl)*100) %>%
  mutate(less.bf.perc = 100 - bf.perc) %>%
  summarize(n.min = min(less.bf.perc),
            n.max = max(less.bf.perc))

#3
# For example, within stem position and wood species identity explained approximately equivalent amounts of variation (~ X%) in bacterial composition in two data subsets (overlap 2014 Fig 2a, deployment 5 yrs Fig 2b).
varp.df1 %>%
  filter(gene == "16S") %>%
  filter(model == "ro") %>%
  filter(yrharv == 2014) %>%
  filter(terms %in% c("species4","logLoc")) %>%
  mutate(perc.var.expl = Adj.R.square * 100)
varp.df1 %>%
  filter(model == "rd") %>%
  filter(subset == "201416S_d09") %>%
  filter(terms %in% c("species4","logLoc")) %>%
  mutate(perc.var.expl = Adj.R.square * 100)


#4
# Whereas wood species accounted for X% of variation in fungi after 1 year (2011-12), it accounted for X-X% after 3 years, and just X% after 5 years (2009-14) (Fig 2).
varp.df1 %>%
  filter(gene == "ITS") %>%
  filter(model == "rd") %>%
  filter(terms %in% c("species4")) %>%
  mutate(perc.var.expl = Adj.R.square * 100)

# This pattern was less pronounced for bacterial communities.
varp.df1 %>%
  filter(gene == "16S") %>%
  filter(model == "rd") %>%
  filter(terms %in% c("species4")) %>%
  mutate(perc.var.expl = Adj.R.square * 100)


#5
# Relative ability of traits to explain community variation in response to wood species
varp.df1 %>%
  filter(terms %in% c('woodtraits','species4')) %>%
  filter(model != "ro") %>%
  select(subset, terms, Adj.R.square) %>%
  spread(key = terms, value = Adj.R.square) %>%
  mutate(perc = (woodtraits/species4)*100)

#6
# how much did fungi structure bacteria?
varp.df1 %>%
  filter(terms %in% c("Co-clade")) %>%
  select(subset, gene, Adj.R.square) %>%
  spread(key = gene, value = Adj.R.square) %>%
  mutate()


# -------------------------------------------------------------------#
# Do in-text calculations based on angio/gymno varpart model

rd.ag <- readRDS(file = "data_Rsynth/rotplot_intermediates/dbrda_ag_d.RData")

#1
# Identification of wood species as either angio- or gymnosperms explained X-X% of variation in fungal communities and X-X% of variation in bacterial communities
var.fract <- list_to_df(lapply(rd.ag, function(x){x$var.fract}))
var.fract  %>%
  filter(!grepl("sapro", source)) %>%
  filter(!grepl("basidio", source)) %>%
  filter(terms == "angio.gymno") %>%
  mutate(perc.var.expl = Adj.R.square * 100) %>%
  separate(source, into = c("yrharv","other"), sep = 4, remove = F) %>%
  separate(other, into = c("gene","yrdep"), sep = "_") %>%
  group_by(gene) %>%
  summarize(nmin = min(perc.var.expl),
            nmax = max(perc.var.expl))

#2
# Differences in fungal composition between angio- and gymnosperm hosts were relatively consistent after 1 year, 3 years, and 5 years of decay; whereas, bacterial communities on angio- and gymnosperms became more similar with time (Fig. S7).
dist.df <- list_to_df(lapply(rd.ag, function(x){x$dist.df}))
dist.df %>%
  filter(!grepl("sapro", source)) %>%
  filter(!grepl("basidio", source)) %>%
  filter(terms == "angio.gymno") %>%
  separate(source, into = c("yrharv","other"), sep = 4, remove = F) %>%
  separate(other, into = c("gene","yrdep"), sep = "_")

