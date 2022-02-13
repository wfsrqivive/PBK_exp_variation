model <- function (pars){#input required for the solver
  
  derivs <- function(t, state, parms) { # returns rate of change
    with(as.list(c(state, parms)), {
      
      Vad = BW*FVad	 
      Vbo = BW*FVbo        	 
      Vbr = BW*FVbr		
      Vgu = BW*FVgu          	 
      Vhe = BW*FVhe        	
      Vki = BW*FVki
      Vli = BW*FVli 
      Vlu = BW*FVlu
      Vmu = BW*FVmu
      Vsk = BW*FVsk 
      Vsp = BW*FVsp
      Vte = BW*FVte
      Vve = BW*FVve
      Var = BW*FVar
      Vpl = BW*FVpl 
      Vrb = BW*FVrb 
      Vre = BW*FVre
      
      Qad = QC*FQad            
      Qbo = QC*FQbo             
      Qbr = QC*FQbr             
      Qgu = QC*FQgu            
      Qhe = QC*FQhe             
      Qki = QC*FQki             
      Qh = QC*FQh              
      Qha = QC*FQh - QC*FQgu - QC*FQsp 
      Qlu = QC*FQlu		
      Qmu = QC*FQmu    
      Qsk = QC*FQsk    
      Qsp = QC*FQsp    
      Qte = QC*FQte    
      Qre = QC*FQre		
      
      Cadipose = Aad/Vad		
      Cbone = Abo/Vbo		 
      Cbrain = Abr/Vbr	
      Cgut = Agu/Vgu		
      Cheart = Ahe/Vhe	 
      Ckidney = Aki/Vki		
      Cliver = Ali/Vli	 
      Clung = Alu/Vlu			 
      Cmuscle = Amu/Vmu		
      Cskin = Ask/Vsk		 
      Cspleen = Asp/Vsp		 
      Ctestes = Ate/Vte		 
      Cvenous = Ave/Vve		
      Carterial = Aar/Var	
      Crest = Are/Vre 		 
      Cplasmavenous = Cvenous/BP
      CLmet = (CLint/fuhep)*SF*Vli*60/1000
      
      
      Cliverfree = Cliver*fup #/Kpli		
      Ckidneyfree = Ckidney*fup
      Cplasmavenousfree = Cplasmavenous*fup
      
      
      
      #differential equations body
      dD = - Ka*D
      dAad = Qad*(Carterial - Cadipose/Kpad*BP)										
      dAbo = Qbo*(Carterial - Cbone/Kpbo*BP)										
      dAbr = Qbr*(Carterial - Cbrain/Kpbr*BP)										
      dAgu = Ka*D + Qgu*(Carterial - Cgut/Kpgu*BP) #			
      dAhe = Qhe*(Carterial - Cheart/Kphe*BP)										
      dAki = Qki*(Carterial - Ckidney/Kpki*BP) - CLrenal*Cplasmavenousfree							
      dAli=  Qha*Carterial + Qgu*(Cgut/Kpgu*BP) + Qsp*(Cspleen/Kpsp*BP) - Qh*(Cliver/Kpli*BP) - Cliverfree*CLmet	 #sumUptake+
      dAlu = Qlu*Cvenous - Qlu*(Clung/Kplu*BP)
      dAmu = Qmu*(Carterial - Cmuscle/Kpmu*BP)
      dAsk = Qsk*(Carterial - Cskin/Kpsk*BP)	
      dAsp = Qsp*(Carterial - Cspleen/Kpsp*BP)
      dAte = Qte*(Carterial - Ctestes/Kpte*BP)
      dAve = Qad*(Cadipose/Kpad*BP) + Qbo*(Cbone/Kpbo*BP) + Qbr*(Cbrain/Kpbr*BP) +  Qhe*(Cheart/Kphe*BP) + Qki*(Ckidney/Kpki*BP) + Qh*(Cliver/Kpli*BP)  + Qmu*(Cmuscle/Kpmu*BP) + Qsk*(Cskin/Kpsk*BP) + Qte*(Ctestes/Kpte*BP) + Qre*(Crest/Kpre*BP) - Qlu*Cvenous 						
      dAar = Qlu*(Clung/Kplu*BP) - Qlu*Carterial
      dAre = Qre*(Carterial - Crest/Kpre*BP)		
  
      
      #{Defining amount metabolitezed and cleared by the kidney for the mass balance equation}
      dAliClearance =   Cliverfree*CLmet
      dAkiClearance =   CLrenal*Cplasmavenousfree
      
      #{Mass Balance}
      Total = dose
      Calculated = Aad+Abo+Abr+Agu+Ahe+Aki+Ali+Alu+
        Amu+Ask+Asp+Ate+Ave+Aar+Are+
        D+AliClearance+AkiClearance
      ERROR=((Total-Calculated)/Total+1E-30)*100
      MASSBBAL=Total-Calculated + 1 
      
      
      
      res<-c(dD, dAad, dAbo, dAbr, dAgu, 
             dAhe, dAki, dAli, dAlu, 
             dAmu, dAsk, dAsp, dAte, 
             dAve, dAar, dAre, 
             dAliClearance, dAkiClearance) 
      
      return(list(res, ERROR = ERROR, 
                  MASSBBAL = MASSBBAL,
                  Calculated = Calculated,
                  Cplasmavenous = Cplasmavenous,#model of Jones and Rowland is in mg/L
                  Cadipose = Cadipose, 
                  Cbone = Cbone,	 
                  Cbrain = Cbrain,	
                  Cgut = Cgut,		
                  Cheart = Cheart,	 
                  Ckidney = Ckidney,		
                  Cliver = Cliver,	 
                  Clung = Clung,			 
                  Cmuscle = Cmuscle,		
                  Cskin = Cskin,	 
                  Cspleen = Cspleen,		 
                  Cgonads = Ctestes,		 
                  Cvenous = Cvenous,		
                  Carterial = Carterial 
      ))
      
    }) #end with
  } #end derivs
  
  solvePBPKCmax<- as.data.frame(ode(y = c(D= 0.1*70*1, Aad=0, Abo=0, Abr=0, Agu=0, 
                                          Ahe=0, Aki=0, Ali=0, Alu=0, 
                                          Amu=0, Ask=0, Asp=0, Ate=0, 
                                          Ave=0, Aar=0, Are=0, 
                                          AliClearance = 0, AkiClearance = 0),
                                    times  = seq(0, 24, by = 0.1),
                                    parms = pars,
                                    atol=1e-50,
                                    func = derivs)) #%>%
    #mutate(AUC =  trapz(time,Cplasmavenous)) %>%
    #filter(AUC == max(AUC))
    #filter(Cplasmavenous == max(Cplasmavenous)) 
  return(solvePBPKCmax)
  
}

