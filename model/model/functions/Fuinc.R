##### fraction unbound in the in vitro hepatic clearance incubation according to Kilford et al. 2008

f.Fuinc<-function(method.Clint, frac.unionized, frac.ionized.acid, frac.ionized.base, logD, logP, conc.enzymes.Clint){
  #hepatocytes (Kilford et al., 2008)
  if (method.Clint == "hep" & frac.unionized >= 0.95 || method.Clint == "hep" & frac.ionized.acid >= 0.95){
    Fuinc<-1/(1+125*(conc.enzymes.Clint/1*0.005)*10^(0.072*logD^2+0.067*logD-1.126))
  } else if(method.Clint == "hep" & frac.ionized.base >= 0.95 ){
    Fuinc<-1/(1+125*(conc.enzymes.Clint/1*0.005)*10^(0.072*logP^2+0.067*logP-1.126))
  } else if(method.Clint == "hep" & frac.unionized < 0.95 || method.Clint == "hep" &frac.ionized.acid < 0.95 || method.Clint == "hep" & frac.ionized.base < 0.95){
    Fuinc<-(frac.ionized.acid+frac.unionized)*(1/(1+125*(conc.enzymes.Clint/1*0.005)*10^(0.072*logD^2+0.067*logD-1.126))) + frac.ionized.base*(1/(1+125*(conc.enzymes.Clint/1*0.005)*10^(0.072*logP^2+0.067*logP-1.126))) 
  #microsomes or S9 (Kilford et al., 2008 based on Halifax and Houston 2006)
  } else if (method.Clint == "mic" & frac.unionized >= 0.95 || method.Clint == "mic" & frac.ionized.acid >= 0.95||method.Clint == "S9" & frac.unionized >= 0.95 || method.Clint == "S9" & frac.ionized.acid >= 0.95){
    Fuinc<-1/(1+conc.enzymes.Clint*10^(0.072*logD^2+0.067*logD-1.126)) 
  } else if(method.Clint == "mic" & frac.ionized.base >= 0.95 ||method.Clint == "S9" & frac.ionized.base >= 0.95 ){
    Fuinc<-1/(1+conc.enzymes.Clint*10^(0.072*logP^2+0.067*logP-1.126))
  } else if(method.Clint == "mic" & frac.unionized < 0.95 || method.Clint == "mic" &frac.ionized.acid < 0.95 || method.Clint == "mic" & frac.ionized.base < 0.95 || method.Clint == "S9" & frac.unionized < 0.95 || method.Clint == "S9" &frac.ionized.acid < 0.95 || method.Clint == "S9" & frac.ionized.base < 0.95){
    Fuinc<-(frac.ionized.acid+frac.unionized)*(1/(1+conc.enzymes.Clint*10^(0.072*logD^2+0.067*logD-1.126))) + frac.ionized.base*(1/(1+conc.enzymes.Clint*10^(0.072*logP^2+0.067*logP-1.126))) 
  #else  
  } else {
    Fuinc<-1
  }
  return(Fuinc)
}

#f.Fuinc(method.Clint ="S9", frac.unionized =100, frac.ionized.acid = 0, frac.ionized.base = 0, logD =3, logP = 3, conc.enzymes.Clint =3)
## calculations microsomes need to be added, checke with Simcyp app ##
## calculations microsomes need to be added, checke with Simcyp app ##