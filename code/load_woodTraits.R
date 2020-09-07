#wood.traits <- load_fungalEnvirTraits()
#write.csv(wood.traits, file="data_Rsynth/fungEnvWoodTraits.csv")

###############

addTracheid <- function(traits){
  
  #Integrate Oyomoare's tracheid data into full dataset
  
  #Oyomoare's data
  #o.trait.data<-read.csv("data/traits/OyomoareAnatomyData/tisprop_21sp_2015Sept5.csv")
  #o.trachied1<-read.csv("data/traits/OyomoareAnatomyData/trachieddiameter_2015Sept5earlylatewood.csv")
  o.trachied2<-read.csv("data/traits/OyomoareAnatomyData/trachieddiameter_2015Sept6.csv")
  
  #make a species column
  o.trachied2$species4<-NA
  o.trachied2[grepl("JUVI", o.trachied2$Species_ind),"species4"]<-"JUVI"
  o.trachied2[grepl("PIEC", o.trachied2$Species_ind),"species4"]<-"PIEC"
  o.trachied2[grepl("PIST", o.trachied2$Species_ind),"species4"]<-"PIST"
  
  #average by species
  o.trachied2 %>%
    group_by(species4) %>%
    summarize(Conduit.diam=mean(FeretsDiameter..um.)) -> o.trachied.summ
  
  #insert into column "Conduit.D.um."
  traits[traits$Species == "JUVI","Conduit.D.um."]<-
    o.trachied.summ[o.trachied.summ$species4=="JUVI","Conduit.diam"]
  traits[traits$Species == "PIEC","Conduit.D.um."]<-
    o.trachied.summ[o.trachied.summ$species4=="PIEC","Conduit.diam"]
  traits[traits$Species == "PIST","Conduit.D.um."]<-
    o.trachied.summ[o.trachied.summ$species4=="PIST","Conduit.diam"]
  
  return(traits)
}

addparenchyma <- function(traits){
  
  #Integrate O's parenchyma data into full dataset
  
  #Oyomoare's data
  o.trait.data<-read.csv("data/traits/OyomoareAnatomyData/tisprop_21sp_2015Sept5.csv")
  #o.trachied1<-read.csv("data/traits/OyomoareAnatomyData/trachieddiameter_2015Sept5earlylatewood.csv")
  #o.trachied2<-read.csv("data/traits/OyomoareAnatomyData/trachieddiameter_2015Sept6.csv")
  
  o.trait.data %>%
    separate(Image, into=c("Species","rep"), sep= -1) %>%
    group_by(Species) %>%
    summarize(parenchymaFrac=mean(Parenchymafrac)) %>%
    left_join(traits) -> traits
  
  return(traits)
}

load_traits <- function(){
  
  traits <- read.csv("data/traits/PlantTraitsDatabase0.1.3.b.csv")
  
  # add species means from Oyomoare's data
  traits <- addTracheid(traits)
  traits <- addparenchyma(traits)
  
  # calculate C fractions
  traits$cellulose<-traits$WADFPerDM-traits$WLigPerDM
  traits$hemicellulose<-traits$WNDFPerDM-traits$WADFPerDM
  traits$fibre_CP<-traits$WNDFPerDM/traits$WCProtPerDM
  traits$lignin_CP<-traits$WLigPerDM/traits$WCProtPerDM
  
  return(traits)
  
}

load_fungalEnvirTraits <- function(){
  
  require('tidyverse')
  
  traits <- load_traits()
  traits.meta<-read.csv("data/traits/readme2__PlantTraitsDatabase0.1.3.b.csv")
  
  ############
  # Subset wood traits that describe the fungal environment
  
  #find relevant column names used the meta data
  wood.traits<-traits.meta[traits.meta$tissue=='wood' & !is.na(traits.meta$tissue),]
  
  #select key dimensions of fungal environment
  #start with anything in units of DM
  wood.traits.DM<-wood.traits[grepl("%DM", wood.traits$unit),]
  exclude<-c("Wd.Per_N","WPerDM")
  wood.traits.DM<-wood.traits.DM[!wood.traits.DM$colName %in% exclude,]
  #add density, conduit length and diameter, parenchymaFrac
  add<-c("Conduit.D(um)","Conduit.L(m)","StemTopDensity","StemBottomDensity","WMnPPM")
  wood.traits.add<-wood.traits[wood.traits$colName %in% add,]
  wood.traits.select<-rbind(wood.traits.DM, wood.traits.add)
  
  #fix the change in colnames when they become real column names
  colName1<-gsub("\\(","\\.",wood.traits.select$colName)
  colName2<-gsub("\\)","\\.", colName1)
  
  #add the newly calculated tissue fractions
  colName3<-c(colName2, "cellulose","hemicellulose","fibre_CP","lignin_CP", "parenchymaFrac")
  
  #select the columns in the trait data based on the selected traits
  select.traits<-data.frame(Species=traits[,"Species"],
                            traits[,colnames(traits) %in% colName3])
  
  ############
  # Check out trait correlations
  select.traits1<-select.traits[complete.cases(select.traits),]
  
  #look at the correlation among these traits
  corobj<-cor(select.traits1[,-1])
  #library(corrplot)
  #par(mfrow=c(1,1))
  #corrplot(corobj, method="pie",type="upper")
  
  #### correlations ####
  # StemTopDensity and StemBottomDensity (+)
  # WCProtPerDM and WAdjProtPerDM
  # WLigPerDM and WLigPerNDF
  # Wd.Per_N.2. and WCProtPerDM and WCAdjProtPerDM
  # Wd.Pre_C.2. and WLigPerDM and WLigPerNDF
  select.traits2<-select.traits1[,!colnames(select.traits1) %in%
                                   c("StemTopDensity","WAdjProtPerDM","WLigPerNDF","Wd.Per_N.2.","Wd.Per_C.2.")]
  
  #look at the correlation among these traits
  corobj<-cor(select.traits2[,-1])
  #corrplot(corobj, method="pie",type="upper")
  
  #### correlations ####
  # WCProtPerDM and WSolProtPerDM
  # WSolProtPerDM and WSolProtPerCP
  # WADFPerDM and WADFPerNDF
  # WADFPerDM and WNDFPerDM
  # WAshPerDM and WCaPerDM
  # WPPerDM and WKPerDM
  select.traits3<-select.traits2[,!colnames(select.traits2) %in%
                                   c("WSolProtPerDM","WADFPerNDF","WKPerDM")]
  
  #look at the correlation among these traits
  corobj<-cor(select.traits3[,-1])
  #corrplot(corobj, method="pie",type="upper")
  
  select.traits4<-select.traits3[,!colnames(select.traits3) %in%
                                   c("fibre_CP","lignin_CP","WSolProtPerCP","WADFPerDM","WAshPerDM","WNaPerDM","WNDFPerDM",
                                     "WMoistPerDM","WMgPerDM")]
  #look at the correlation among these traits
  corobj<-cor(select.traits4[,-1])
  #corrplot(corobj, method="pie",type="upper")
  
  
  ############
  #Tidy up
  
  #put back in the incomplete cases
  wood.traits<-traits[,colnames(traits) %in% colnames(select.traits4)]
  
  #simplify trait names
  wood.traits %>%
    rename("conduitDiam" = `Conduit.D.um.`,
           "conduitLeng" = `Conduit.L.m.`, 
           "density" = `StemBottomDensity`, 
           "crudeProt" = `WCProtPerDM`,
           "lignin" = `WLigPerDM`,
           "Ca" = `WCaPerDM`, 
           "P" = `WPPerDM`, 
           "Mn" = `WMnPPM`) -> wood.traits
  
  #standardize species codes
  
  #replace PRSE with PRSE/PRVI 
  wood.traits$Species<-as.character(wood.traits$Species)
  wood.traits[wood.traits$Species == "PRSE","Species"]<-"PRSE/PRVI" #rename the problem species
  wood.traits %>%
    rename("species4" = `Species`) -> wood.traits
  
  return(wood.traits)
}

