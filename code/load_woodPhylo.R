
fix_problem_species <- function(tree, prob_species, dontchoose = wood_names){
  
  require(taxonlookup)
  for(i in 1:length(prob_species)){
    
    if(prob_species[i] == "Quercus_velutina"){
      replace <- "Quercus_rubra"
    }
    
    if(!prob_species[i] %in% c("Quercus_velutina")){
      genus <- taxonlookup:::split_genus(prob_species[i])
      replace <- sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    }
    
    while (replace%in%dontchoose) replace<-sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    tree$tip.label[tree$tip.label==replace]<-prob_species[i]
  }
  return(tree)
}

load_zanne_tree <- function(){
  
  require(phytools)
  require(diversitree)
  
  zae <- read.tree("data/zanne1.1.tre")
  logTrt <- load_logTrt() # this fxn is in code/rotplot_microbes/load_rotplotMetaData.R
  wood_names <- unique(logTrt$binomial[!is.na(logTrt$binomial)])
  set.seed(42)
  
  length(wood_names)
  zae_mod <- fix_problem_species(zae, 
                               prob_species = c("Carya_tomentosa","Quercus_velutina"), 
                               dontchoose = wood_names)
  
  zanneTree <- diversitree:::drop.tip.fixed(phy = zae_mod,
                                          zae_mod$tip.label[!zae_mod$tip.label%in%wood_names])

  return(zanneTree)
  
}
