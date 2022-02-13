#calculate tissue:plasma partition coefficients based on: Rodgers and Rowland http://jpharmsci.org/article/S0022-3549(16)31789-0/fulltext and http://jpharmsci.org/article/S0022-3549(16)32034-2/fulltext

calcKp_RR <- function(logP, pKa1, pKa2, pKa3,Fup, BP, type, dat){

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma
  
  
  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  Ka <- 10^(-pKa1)
  HCT <- 0.45 #hematocrit

  
  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT*Fup)
  
  X <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_IW-pKa1),
              #3-monoprotic base
              10^(pKa1-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa1)+10^(2*pH_IW-pKa1-pKa2),
              #5-diprotic base
              10^(pKa2-pH_IW)+10^(pKa1+pKa2-2*pH_IW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa2-pH_IW)+10^(pH_IW-pKa1),  
              #7-triprotic acid
              10^(pH_IW-pKa1)+10^(2*pH_IW-pKa1-pKa2)+10^(3*pH_IW-pKa1-pKa2-pKa3),  
              #8-triprotic base
              10^(pKa3-pH_IW)+10^(pKa3+pKa2-2*pH_IW)+10^(pKa1+pKa2+pKa3-3*pH_IW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa3-pH_IW)+10^(pH_IW-pKa1)+10^(2*pH_IW-pKa1-pKa2), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa1)+10^(pKa3-pH_IW)+10^(pKa2+pKa3-2*pH_IW))       
  
  Y <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_P-pKa1),
              #3-monoprotic base
              10^(pKa1-pH_P), 
              #4-diprotic acid
              10^(pH_P-pKa1)+10^(2*pH_P-pKa1-pKa2),
              #5-diprotic base
              10^(pKa2-pH_P)+10^(pKa1+pKa2-2*pH_P), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa2-pH_P)+10^(pH_P-pKa1),  
              #7-triprotic acid
              10^(pH_P-pKa1)+10^(2*pH_P-pKa1-pKa2)+10^(3*pH_P-pKa1-pKa2-pKa3),  
              #8-triprotic base
              10^(pKa3-pH_P)+10^(pKa3+pKa2-2*pH_P)+10^(pKa1+pKa2+pKa3-3*pH_P),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa3-pH_P)+10^(pH_P-pKa1)+10^(2*pH_P-pKa1-pKa2), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa1)+10^(pKa3-pH_P)+10^(pKa2+pKa3-2*pH_P))       
  
  Z <- switch(type,
              #1-neutral
              1,   
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa1-pH_RBC), 
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa2-pH_RBC)+10^(pKa1+pKa2-2*pH_RBC), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa2-pH_RBC)+10^(pH_RBC-pKa1),  
              #7-triprotic acid
              1,  
              #8-triprotic base
              10^(pKa3-pH_RBC)+10^(pKa3+pKa2-2*pH_RBC)+10^(pKa1+pKa2+pKa3-3*pH_RBC),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa3-pH_RBC)+10^(pH_RBC-pKa1)+10^(2*pH_RBC-pKa1-pKa2), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa1)+10^(pKa3-pH_RBC)+10^(pKa2+pKa3-2*pH_RBC)) 
  
  
  Ka_PR <- (1/Fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
  
  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa1>7) | (type==5 & pKa1 >7) | (type==6 & pKa2 > 7) | (type==8 & pKa1 > 7) | (type==9 & pKa3>7) | (type==10 & pKa2>7), 1,2)
  
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  
  
  # Multiply by Fup to get Kp rather than Kpu
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*Fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*Fup  #lipid
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*Fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*Fup #lipid
  }else{  #neutrals
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*Fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))*Fup  #lipid
  }
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms

  return(list(Kp))
  
}

#calcKp_RR(logP=2.61, pKa1=11.82, pKa2=7.29 , pKa3= NA,Fup = 0.3, BP=1, type =6, dat = tissue_comp_RR)






