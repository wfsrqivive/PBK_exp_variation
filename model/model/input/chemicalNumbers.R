library(tidyverse)
CLint_data<-read.csv("invitro_CLint_input.csv") %>%
  group_by(compound) %>%
  summarize(meanCLint = mean(Clint), 
            sdCLint = sd(Clint),
            cvCLint = sdCLint/meanCLint*100,
            nCLint = n()) %>%
  filter(nCLint>1)

Papp_data<-read.csv("invitro_Papp_input.csv") %>%
  group_by(compound) %>%
  summarize(meanPapp = mean(PappAB), 
            sdPapp = sd(PappAB),
            cvPapp = sdPapp/meanPapp*100,
            nPapp = n()) %>%
  mutate(meanPapp = ifelse(is.na(sdPapp), NA, meanPapp)) %>%
  filter(nPapp>1)

Fup_data<-read.csv("invitro_Fup_input.csv") %>%
  group_by(compound) %>%
  summarize(meanFup = mean(Fup), 
            sdFup = sd(Fup),
            cvFup = sdFup/meanFup*100,
            nFup = n()) %>%
  mutate(meanFup = ifelse(is.na(sdFup), NA, meanFup))%>%
  filter(nFup>1)

data_set<-full_join(CLint_data, Papp_data) %>%
  full_join(., Fup_data) %>%
  arrange(meanCLint,meanPapp, meanFup) %>%
  mutate_if(is.numeric, funs(signif(.,2))) %>% 
  mutate_at(vars(contains('cv')), funs(signif(.,2))) %>%
  tibble::rownames_to_column(., "number") %>%
  select(-sdCLint, -sdPapp,-sdFup)

write.csv(data_set, "data_set.csv")

  