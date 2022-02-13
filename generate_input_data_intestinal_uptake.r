library(tidyverse)
library(httk)

#f.input.data<-function(){
######Intestinal uptake######
R<-1 #radius small intestine
Tsi<-3.32 #small intestine transit time 
kt<-1/(Tsi/7)
#background:
#https://pubmed.ncbi.nlm.nih.gov/10486429/
#https://pub.iapchem.org/ojs/index.php/admet/article/view/638

invitro.KaFa<-read.csv("./input/invitro_Papp_input.csv") %>%
  mutate(Peff.cm.power.minus.four.per.s = (10^(0.4926*log10(PappAB)- 0.1454)))%>%
  mutate(Peff.cm.per.hr = Peff.cm.power.minus.four.per.s/10000*3600)%>%
  mutate(Ka= (Peff.cm.per.hr*2)/R) %>%
  mutate(Fa = 1-(1+Ka/kt)**-7) %>%
  mutate(method.KaFa = reference) %>%
  select(compound, CAS, Ka, Fa, method.KaFa, PappAB)


KaFa<-invitro.KaFa

write.csv(KaFa, "./input/inputDataFile_intestinal_uptake.csv")
