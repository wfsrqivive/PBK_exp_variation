library(tidyverse)
library(httk)


######Hepatic clearance######
clearance<-read.csv("./input/inputDataFile_intrinsic_clearance.csv")%>%
  filter(compound %in% c("Diazepam", "Caffeine","Verapamil", "Quinidine", "Midazolam", "Diltiazem")) %>%
  select(-X, -method.Fup, -method.phys.chem.Fup, 
         -Fup, -Fut, -reference.Fup, -original.reference.Fup) 

KaFa<-read.csv("./input/inputDataFile_intestinal_uptake.csv") %>%
  filter(compound %in% c("Diazepam", "Caffeine","Verapamil", "Quinidine", "Midazolam", "Diltiazem")) %>%
  select(compound, method.KaFa,Ka, Fa, Papp)

Fup<-read.csv("./input/invitro_Fup_input.csv") %>%
  filter(compound %in% c("Diazepam", "Caffeine","Verapamil", "Quinidine", "Midazolam", "Diltiazem")) %>%
  select(compound, method.Fup, Fup)

PBK.input.data<-full_join(clearance, KaFa, by = c("compound")) %>%
  full_join(., Fup, by = c("compound")) %>%
  filter(!is.na(Ka)) 

write.csv(PBK.input.data, "./input/inputDataFile_combined_individual.csv")
