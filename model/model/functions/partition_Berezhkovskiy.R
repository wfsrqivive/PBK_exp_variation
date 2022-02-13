#calculate tissue:plasma partition coefficients based on BEREZHKOVSKIY, LEONID M. (2004) equations

calcKp_Berez <- function(logP, pKa1, pKa2, pKa3, Fup, BP, type, dat){
  
  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))
  
  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)
  
  Vwp <- dat$f_water[dat$tissue == "Plasma"]
  Vnlp <- dat$f_n_l[dat$tissue == "Plasma"]
  Vphp <- dat$f_pl[dat$tissue == "Plasma"]
  
  dat2 <- dat %>% filter(!tissue %in% c("Plasma","RBCs"))
  
  Vwt <- dat2$f_water[dat2$tissue != "Adipose"]
  Vwad <- dat2$f_water[dat2$tissue == "Adipose"]
  Vnlt <- dat2$f_n_l[dat2$tissue != "Adipose"]
  Vnlad <- dat2$f_n_l[dat2$tissue == "Adipose"]
  Vpht <- dat2$f_pl[dat2$tissue != "Adipose"]
  Vphad <- dat2$f_pl[dat2$tissue == "Adipose"]
  
  fut <- 1/(1+((1-Fup)/Fup)*0.5)
  
  pH <- dat$pH[dat$tissue == "Adipose"]
  #pH <- 7.4 # Use when comparing to PK-Sim Berez. method Kp predictions 
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species
  
  logD_star <- switch(type,
                      #1-neutral
                      logD,   
                      #2-monoprotic acid
                      logD-log10(1+10^(pH-pKa1)),
                      #3-monoprotic base
                      logD-log10(1+10^(pKa1-pH)), 
                      #4-diprotic acid
                      logD-log10(1+10^(2*pH-pKa1-pKa2)),
                      #5-diprotic base
                      logD-log10(1+10^(pKa1+pKa2-2*pH)), 
                      #6-monoprotic acid monoprotic base (acid comes first)
                      logD-log10(1+10^(pKa2-pKa1)),  
                      #7-triprotic acid
                      logD-log10(1+10^(3*pH-pKa1-pKa2-pKa3)),  
                      #8-triprotic base
                      logD-log10(1+10^(pKa1+pKa2+pKa3-3*pH)),  
                      #9-diprotic acid monoprotic base (first two are acid)
                      logD-log10(1+10^(pH-pKa1-pKa2+pKa3)), 
                      #10-diprotic base monoprotic acid (first one is acid)
                      logD-log10(1+10^(pKa2+pKa3-pKa1-pH)))       
  
  D_star <- 10^logD_star   
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+((Vwad/fut)+0.7*Vphad))/(D_star*(Vnlp+0.3*Vphp)+((Vwp/Fup)+0.7*Vphp)))

  
  P <- 10^logP
  Kpt <- ((P*(Vnlt+0.3*Vpht)+((Vwt/fut)+0.7*Vpht))/(P*(Vnlp+0.3*Vphp)+((Vwp/Fup)+0.7*Vphp))) 
  
  #Kp <- c(Kpad, Kpt)
  # name <- dat2$tissue %>% substr(1,2) %>% tolower()
  # name <- paste("Kp", name, sep="")
  # uParam <- split(Kp, name)
  # 
  # return(uParam)
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  # return(nms)
  Kp <- as.list(c(Kpad,Kpt))
  names(Kp) <- nms
  
  
  return(list(Kp))
}

#calcKp_SchmittOutcome
#tissue_comp<-read.csv("input_data/data/unified_tissue_comp.csv")
#test<-data.frame(logP =3.65, pKa1 =4.54, pKa2= NA, pKa3 = NA, Fup = 0.018, type = 2)
#test<-test%>% 
#  mutate(calcKp_Berez= calcKp_Berez(logP, pKa1, pKa2, pKa3, Fup, BP=1, type,  dat = tissue_comp)) %>%
#  unnest_wider(calcKp_Berez)

