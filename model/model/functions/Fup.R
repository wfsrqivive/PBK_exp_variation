##### fraction unbound in the in vitro hepatic clearance incubation according to Kilford et al. 2008

#Method according to Lobell & Sivarajah, 2003. Molecular Diversity, 7: 69-87.
#the equation for Fup should be used in combination with "ionizationCalculations.R" function

f.Fup<-function(ionization, frac.unionized, frac.ionized.acid, frac.ionized.base, logP, quat.nitrogen ){
  #permanently positively charged (contains quaternary nitrogen)
    if (quat.nitrogen == "yes") {
    Fup<-0.3978*logP-2.0965
    #uncharged
  } else if (frac.unionized >= 0.90){
    Fup<-0.4485*logP-0.4782  #if 90% is unionized at pH 7.4 than treated as neutral (logD and logP are then the same)
    #negatively charged
  } else if(ionization == "monoproticAcid" || ionization == "diproticAcid" & frac.unionized < 0.90){
    Fup<-frac.unionized*(0.4485*logP-0.4782)+frac.ionized.acid*(0.3649*logP+0.4162)
    #positively charged and AlogP98 > 0.2
  } else if(ionization == "monoproticBase" || ionization == "diproticBase" & frac.unionized < 0.90 & logP>0.2) {
    Fup<-frac.unionized*(0.4485*logP-0.4782)+frac.ionized.base*(0.4628*logP-1.0971)
    #positively charged and AlogP98 < 0.2
  } else if(ionization == "monoproticBase" || ionization == "diproticBase" & frac.unionized < 0.90 & logP<=0.2) {
    Fup<-frac.unionized*(0.4485*logP-0.4782)+frac.ionized.base*-0.954
    
    #zwitterion: if for more than 10% is ionized (acid) with less than 10%  base, then treat as acid
  } else if(ionization == "zwitterionic" & frac.ionized.acid >0.1 & frac.ionized.base<0.1) {
    Fup<-frac.unionized*(0.4485*logP-0.4782)+frac.ionized.acid*(0.3649*logP+0.4162)
    #zwitterion: if for more than 10% is ionized (base) with less than 10% acid, then treat as base
  } else if(ionization == "zwitterionic" & frac.ionized.acid <0.1 & frac.ionized.base>0.1) {
    Fup<-frac.unionized*(0.4485*logP-0.4782)+frac.ionized.base*(0.4628*logP-1.0971)
    #zwitterion: if acid and base are both present in a range between 50%-90%, then treat as true zitterionic
  } else if(ionization == "zwitterionic" & frac.unionized <0.9 & frac.ionized.acid >0.5 || frac.ionized.base>0.5) {
    Fup<--0.477
  } else {
    Fup<-1000000
  }
  return(1/(10^Fup+1))
}

#f.Fup("monoproticAcid",0.6964,0.3036,0, 2.9,"quat.nitrogen")
