#####log D


#reference: Lipophilicity in Drug Action and Toxicology Applications of a Solvation Equation to Drug Transport Properties M. H. Abraham H. S. Chadha
f.logD <- function(ionization, pH, pKa1, pKa2, logP){
  if (ionization=="neutral") {
    logD<-logP
  } else if (ionization =="monoproticAcid") {
    logD<-logP-log10(1+10^(pH-pKa1))
  } else if (ionization =="monoproticBase") {
    logD<-logP-log10(1+10^(pKa1-pH))
  } else if (ionization == "diproticAcid") {
    logD<-logP-log10(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else if (ionization == "diproticBase") {
    logD<-logP-log10(1+10^(pKa2-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    logD<-logP-log10(1+10^(pH-pKa1)+10^(pKa2-pH))
  }
  return(logD) 
}

#####fractions Ionized (acid/base) vs unionized

#reference: Lipophilicity in Drug Action and Toxicology Applications of a Solvation Equation to Drug Transport Properties M. H. Abraham H. S. Chadha
f.frac.unionized <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="neutral") {
    frac.unionized<-1
  } else if (ionization =="monoproticAcid") {
    frac.unionized<-1/(1+10^(pH-pKa1))
  } else if (ionization =="monoproticBase") {
    frac.unionized<-1/(1+10^(pKa1-pH))
  } else if (ionization == "diproticAcid") {
    frac.unionized<-1/(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else if (ionization == "diproticBase") {
    frac.unionized<-1/(1+10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    frac.unionized<-1/(1+10^(pH-pKa1)+10^(pKa2-pH))
  }
  return(frac.unionized) 
}


f.frac.ionized.acid <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="zwitterionic"||ionization=="monoproticAcid") {
    fracIonized.acid<-1-1/(1+10^(pH-pKa1))
  } else if (ionization=="diproticAcid"){
    fracIonized.acid<-1-1/(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else {
    fracIonized.acid<-0
  }
  return(fracIonized.acid) 
}

f.frac.ionized.base <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="zwitterionic") {
    frac.ionized.base<-1-1/(1+10^(pKa2-pH))
  } else if (ionization=="monoproticBase") {
    frac.ionized.base<-1-1/(1+10^(pKa1-pH))
  } else if (ionization=="diproticBase") {
    frac.ionized.base<-1-1/(1+10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    frac.ionized.base<-0
  }
  return(frac.ionized.base) 
}

#f.frac.ionized.acid("diproticAcid", pH=6.5, 7.68, 9.3)

