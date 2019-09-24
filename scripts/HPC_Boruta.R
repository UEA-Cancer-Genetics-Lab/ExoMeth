#This is the HPC-frioendly version of the NanoMeth Rmarkdown to generate the 1,000 resamples for feature selection runs using
# Boruta. This is a _very_ computationally heavy aspect of the work and takes several hours to complete. The results are saved as 
# a RDS object for local computation afterwards. 
### UNDERTAKEN WITH R 3.5 (important for random numbers) ###

#Catch if on the HPC and set working directory (quick hack):
if(getwd() == "/gpfs/home/xjz17vzu")
{setwd("MovemberMeth")}

######################## LIBRARIES ########################
#Cleaning and munging:
library(readxl)
library(tidyverse)
library(magrittr)
library(lubridate)

#Modelling:
library(Boruta)
library(randomForest)
library(glmnet)
library(caret)
library(yardstick)
library(rsample)

#Computation
library(furrr)

#Making things pretty:
library(formattable)
library(qwraps2)
library(cowplot)
library(ggpubr)
library(ggsci)
library(dabestr)
library(UpSetR)
library(cowplot)

#Settings and stuff:
options(qwraps2_markup = "markdown")
plan(multiprocess)

#knitr setup
knitr::opts_chunk$set(echo = FALSE)
set.seed(2903)

######################## DATA IMPORT ########################
# This chunk loads all the data required and sets some parameters for filtering and general munging of data
# The end products are the cleaned datasets ready for analysis and modelling.

#Set somne parameters that can be tweaked if necessary:
PSA_cutoff <- 100
Time_cutoff <- 39 #Number of weeks

#### NANOSTRING DATA ####
#Grab all the nanoString data:
NanoString <- read_csv("data/raw/NanoString/Nano167_Log2Pos_Norm.csv",
                       col_names = TRUE,
                       col_types = cols(
                         .default = col_double(),
                         Sample_ID = col_character()
                       )) %>% 
  rename("ERG3 exons 4-5" = "ERG 3 ex 4-5",
         "ERG3 exons 6-7" = "ERG3 ex 6-7",
         TIMP4 = "Timp4")
### DUBLIN DATA ####
Perry_Meth  <- read_csv("data/processed/Perry_Methylation_Clinical_Clean.csv",
                        #Ugly but safe way to specify columns:
                        col_types = c("cnnnnnnffnnffnfDDfffffff")) %>% 
  filter(PSA < PSA_cutoff) %>% 
  #mutate the dates of urine and biopsy into date formats
  mutate(UrineDate = ymd(UrineDate),
         BxDate    = ymd(BxDate),
         #calculate the difference in weeks
         TimeDiff  = difftime(UrineDate, BxDate, units = "weeks")) %>% 
  mutate(TimeDiff = replace_na(TimeDiff, 0)) %>% 
  #Filter out all samples more than 9 months (39 weeks) different
  filter(abs(TimeDiff) <= Time_cutoff) %>% 
  select(-TimeDiff, -UrineDate, -BxDate) %>% 
  #set all the Perry Probes to their REAL NAMES
  rename(mGSTP1  = GSTP1_AP,
         mAPC    = APC_AP,
         mSFRP2  = PG3,
         mIGFBP3 = PG4,
         mIGFBP7 = PG5,
         mPTSG2  = PG6)

### TORONTO DATA ####
Bapat_Meth  <- read_csv("data/processed/Bapat_Methylation_Clinical_Clean.csv",
                        #Ugly but safe way to specify columns:
                        col_types = c("cnnnnnnffnnffnfDDfffffff")) %>% 
  filter(PSA < PSA_cutoff) %>% 
  #mutate the dates of urine and biopsy into date formats
  mutate(UrineDate = ymd(UrineDate),
         BxDate    = ymd(BxDate),
         #calculate the difference in weeks
         TimeDiff  = difftime(UrineDate, BxDate, units = "weeks")) %>% 
  mutate(TimeDiff = replace_na(TimeDiff, 0)) %>% 
  #Filter out all samples more than 9 months (39 weeks) different
  filter(abs(TimeDiff) <= Time_cutoff) %>% 
  select(-TimeDiff, -UrineDate, -BxDate)

### Overlap of the two ####
OverlapMeth <- read_csv("data/processed/Overlapping_Methylation_Clinical_Clean.csv",
                        #Ugly but safe way to specify columns:
                        col_types = c("cnnnnnnnnnnnnffnnffnfDDfffffff")) %>% 
  filter(PSA < PSA_cutoff) %>% 
  #mutate the dates of urine and biopsy into date formats
  mutate(UrineDate = ymd(UrineDate),
         BxDate    = ymd(BxDate),
         #calculate the difference in weeks
         TimeDiff  = difftime(UrineDate, BxDate, units = "weeks")) %>% 
  mutate(TimeDiff = replace_na(TimeDiff, 0)) %>% 
  #Filter out all samples more than 9 months (39 weeks) different
  filter(abs(TimeDiff) <= Time_cutoff) %>% 
  select(-TimeDiff, -UrineDate, -BxDate) %>% 
  #set all the Perry Probes to their REAL NAMES
  rename(mGSTP1  = GSTP1_AP,
         mAPC    = APC_AP,
         mSFRP2  = PG3,
         mIGFBP3 = PG4,
         mIGFBP7 = PG5,
         mPTSG2  = PG6)

#Define the probes used by each lab:
Perry_Probes <- c("mGSTP1", "mAPC", "mSFRP2", "mIGFBP3", "mIGFBP7", "mPTSG2")
Bapat_Probes <- c("GSTP1_BB", "APC_BB", "BG3", "BG4", "BG5", "BG6")
outcomes <- c("Cat", "TriSig", "GleaSig", "LowGSig", "ClinSig", "is_C", "is_HC")
clin_vars <- c("PSA", "UrineVol", "DRESize", "Age")
NanoGenes <- NanoString %>% 
  select(-Sample_ID) %>% 
  colnames

#### NanoString Integration ####
#Create the NanoMeth dataset:
NanoPerry <- inner_join(Perry_Meth, NanoString, by = "Sample_ID") %>% 
  select(Sample_ID, Cat:is_HC, Perry_Probes, NanoGenes) %>% 
  filter(!is.na(TriSig)) %>% 
  mutate(TriSig = factor(TriSig, ordered = TRUE),
         DRESize = replace_na(DRESize, "unknown"),
         UrineVol = replace_na(UrineVol, median(.$UrineVol, na.rm = TRUE)))

#Also produce the full overlapping dataset with NanoString, Bapat and Perry probes:
NanoFull <- inner_join(OverlapMeth, NanoString, by = "Sample_ID") %>% 
  select(Sample_ID, Cat:is_HC, Perry_Probes, Bapat_Probes, NanoGenes) %>% 
  filter(!is.na(TriSig)) %>% 
  mutate(TriSig = factor(TriSig, ordered = TRUE),
         DRESize = replace_na(DRESize, "unknown"))

######################## FEATURE SELECTION RUNS ########################
set.seed(2903) #To be sure to be sure
#Produce the many resamples required to run boruta on:
Bor_samples <- bootstraps(NanoPerry, 
                          times = 1000, 
                          strata = "TriSig") 

# Create a helper function to apply boruta to resamples:
Reampled_Boruta <- function(Splits,
                            Variables,
                            Target){
  #Grab the resampled dataframe
  TrainData = Splits$data
  #Apply Boruta
  Bor = Boruta(x = TrainData[, Variables],
               y = TrainData[[Target]],
               maxRuns = 100)
  #Pull out the importances from Boruta and coerce to a tibble
  Importances = Bor[["ImpHistory"]] %>% 
    as_tibble() %>%
    #Get into a long format
    gather(key = "Variable", value = "Importance") %>% 
    #merge with the final decisions made by Boruta, coerced to a dataframe and all nice like:
    left_join(
      Bor[["finalDecision"]] %>% 
        as.data.frame %>% 
        rownames_to_column("Variable") %>% 
        set_colnames(c("Variable", "Decision")), 
      by = "Variable")
  return(Importances)
}

######################## BORUTA RUNS ########################
# This is the meat of the feature reduction; it maps through each comparator set of variables 
# and runs boruta on each of the resamples. It also takes forever. For the final runs this will
# be pushed to the HPC for computation and left to run, with the results saved.
ResampledBoruta <- 
  #for each of the different choices of variables, apply the boruta resample function above
  map(list(clin_vars, Perry_Probes, NanoGenes, 
           c(NanoGenes, Perry_Probes, clin_vars)), 
      function(Variables){
        #Apply and join the boruta function to each resample:
        map_dfr(Bor_samples$splits, function(Resample){
          
          Reampled_Boruta(Splits = Resample, Variables = Variables, Target = "TriSig")
        }) %>% 
          # Now the resampled boruta has run count the number of times each variable is confirmed:
          group_by(Variable) %>% 
          mutate(Proportion = length(which(Decision == "Confirmed"))/length(Decision)) %>% 
          ungroup() %>% 
          #Now give some nice names to the proportions:
          mutate(FinalDecision = case_when(
            str_detect(Variable, "shadow") ~ "Shadow",
            Proportion < 0.4 ~ "Rejected",
            Proportion < 0.80 ~ "Tentative",
            Proportion >= 0.80 ~ "Confirmed")
          )
      }) %>% 
  set_names(c("Clinical", "Methylation", "NanoString", "NanoMeth"))

#Save the results to an RDS for local work:
saveRDS(ResampledBoruta, "output/data_out/ResampledBoruta.RDS")