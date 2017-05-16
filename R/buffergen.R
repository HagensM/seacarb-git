# Copyright (C) 2016 Mathilde Hagens (M.Hagens@uu.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

buffergen <- function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="u74", 
                      gas="potential", warn="y", eos="eos80", long=1.e20, lat=1.e20, NH4t=0, HSt=0) {
  
  n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), 
           length(ks), length(b), length(gas), length(NH4t), length(HSt))
  if(length(flag)!=n){flag <- rep(flag[1],n)}
  if(length(var1)!=n){var1 <- rep(var1[1],n)}
  if(length(var2)!=n){var2 <- rep(var2[1],n)}
  if(length(S)!=n){S <- rep(S[1],n)}
  if(length(T)!=n){T <- rep(T[1],n)}
  if(length(Patm)!=n){Patm <- rep(Patm[1],n)}
  if(length(P)!=n){P <- rep(P[1],n)}
  if(length(Pt)!=n){Pt <- rep(Pt[1],n)}
  if(length(Sit)!=n){Sit <- rep(Sit[1],n)}
  if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
  if(length(kf)!=n){kf <- rep(kf[1],n)}
  if(length(ks)!=n){ks <- rep(ks[1],n)}
  if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
  if(length(b)!=n){b <- rep(b[1],n)}
  #if(length(gas)!=n){gas <- rep(gas[1],n)}
  if(length(NH4t)!=n){NH4t <- rep(NH4t[1],n)}
  if(length(HSt)!=n){HSt <- rep(HSt[1],n)}
  
  # if the concentrations of total silicate, total phosphate, total nitrite,
  # total ammonium, and total hydrogen sulfide are NA, they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  NH4t[is.na(NH4t)] <- 0
  HSt[is.na(HSt)] <- 0
  
  # Calculate acid-base speciation
  Carbfull <- carbfull(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b,
                   gas=gas, warn=warn, eos=eos, long=long, lat=lat, NH4t=NH4t, HSt=HSt)

  #---------------------------------------------------------------------
  #--------------------    buffer effects    ---------------------------
  #---------------------------------------------------------------------
  
  # Define a function to calculate the buffer factors
  bufferfunc <- function(Carbfull, species=species){
  Nb_data <- dim(Carbfull)[1]
  Nb_species <- length(species)
    
  with(as.list(Carbfull),{

    # Calculation of acid-base fractions
    n1  <- NH4 / NH4t
    n2  <- NH3 / NH4t
    b1  <- BOH3 / BOR
    b2  <- BOH4 / BOR
    c1  <- CO2 / DIC
    c2  <- HCO3 / DIC
    c3  <- CO3 / DIC
    p1  <- H3PO4 / Pt
    p2  <- H2PO4 / Pt
    p3  <- HPO4 / Pt
    p4  <- PO4 / Pt
    s1  <- H2S / HSt
    s2  <- HS / HSt
    si1 <- SiOH4 / Sit
    si2 <- SiOOH3 / Sit
    si3 <- SiO2OH2 / Sit
    f1  <- HF / FLUO
    f2  <- F / FLUO
    so1 <- HSO4 / ST
    so2 <- SO4 / ST
    
    # Creation of output data-frames
    dALK.dH <-   data.frame(matrix(vector(), Nb_data, Nb_species))
    dtotX.dH <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dALK.dX <-   data.frame(matrix(vector(), Nb_data, Nb_species))
    dtotX.dX <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dALK.dpH <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dtotX.dpH <- data.frame(matrix(vector(), Nb_data, Nb_species))
    dH.dALK <-   data.frame(matrix(vector(), Nb_data, Nb_species))
    dX.dALK <-   data.frame(matrix(vector(), Nb_data, Nb_species))
    dX.dtotX <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dpH.dALK <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dH.dtotX <-  data.frame(matrix(vector(), Nb_data, Nb_species))
    dpH.dtotX <- data.frame(matrix(vector(), Nb_data, Nb_species))
    
    # Assign names to output data-frames
    colnames(dALK.dH) <- colnames(dtotX.dH) <- colnames(dALK.dX) <- 
      colnames(dtotX.dX) <- colnames(dALK.dpH) <- colnames(dtotX.dpH) <-
      colnames(dH.dALK) <- colnames(dX.dALK) <- colnames(dX.dtotX) <-
      colnames(dpH.dALK) <- colnames(dH.dtotX) <- colnames(dpH.dtotX) <- species
    
    # Define function for calculating ALK.X
    ALK.species <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
         X.name == "DIC") {HCO3 + 2*CO3} else
           if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
             NH3} else
               if(X.name == "H2PO4" || X.name == "H3PO4" || 
                  X.name == "HPO4" || X.name == "PO4" || 
                  X.name == "Pt") {-H3PO4 + HPO4 + 2*PO4} else
                    if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {HS} else
                        if(X.name == "SiOH4" || X.name == "SiOOH3" || 
                           X.name == "SiO2OH2" || X.name == "Sit") {
                          SiOOH3 + 2*SiO2OH2} else
                            if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
                              BOH4} else
                                 if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
                                   -HF} else
                                     if(X.name == "SO4" || X.name == "HSO4" || 
                                        X.name == "ST") {-HSO4} else
                                         NULL
    }
    
    # Define function for determining N
    N.function <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "BOH3" || X.name == "NH4" || 
         X.name == "H2PO4" || X.name == "H2S" || X.name == "SiOH4" || X.name == "F" || 
         X.name == "SO4" || X.name == "DIC" || X.name == "BOR" || 
         X.name == "NH4t" || X.name == "Pt" || 
         X.name == "HSt" || X.name == "Sit" || 
         X.name == "FLUO" || X.name == "ST") {0} else
           
           if(X.name == "HCO3" || X.name == "NH3" || 
              X.name == "HPO4" || X.name == "HS" ||
              X.name == "SiOOH3" || X.name == "BOH4") {1} else
                
                if(X.name == "CO3" || X.name == "PO4" || X.name == "SiO2OH2") {2} else
                     
                     if(X.name == "H3PO4" || X.name == "HF" || X.name == "HSO4") {-1} else
                         
                         NULL
    }
    
    # Define function for calculating dALKdH  
    dALKdH.function <- function(X.name)
    {
      N <- N.function(X.name)
      ALK.X <- ALK.species(X.name)
            
      dALKdH.CO2 <- rep(0.0, Nb_data)
      isDIC <- DIC!=0
        if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
           X.name == "DIC") {
          dALKdH.CO2[isDIC] <- (-1/H[isDIC])*(-N*ALK.X + HCO3[isDIC] + 4*CO3[isDIC])}  else {
            dALKdH.CO2[isDIC] <- (-1/H[isDIC])*(HCO3[isDIC]*(c1[isDIC]-c3[isDIC]) + 2*CO3[isDIC]*(2*c1[isDIC]+c2[isDIC]))
        }

      dALKdH.NH4 <- rep(0.0, Nb_data)
      isNH4t <- NH4t!=0
        if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
          dALKdH.NH4[isNH4t] <- (-1/H[isNH4t])*(-N*ALK.X + NH3[isNH4t])} else {
            dALKdH.NH4[isNH4t] <- (-1/H[isNH4t])*(NH3[isNH4t]*n1[isNH4t])
        }

      dALKdH.H2PO4 <- rep(0.0, Nb_data)
      isPt <- Pt!=0
        if(X.name == "H2PO4" || X.name == "H3PO4" || X.name == "HPO4" || 
           X.name == "PO4" || X.name == "Pt") {
          dALKdH.H2PO4[isPt] <- (-1/H[isPt])*(-N*ALK.X + H3PO4[isPt] + HPO4[isPt] + 4*PO4[isPt])}  else {
            dALKdH.H2PO4[isPt] <- (-1/H[isPt])*(-H3PO4[isPt]*(-p2[isPt]-2*p3[isPt]-3*p4[isPt]) + 
                                     HPO4[isPt]*(2*p1[isPt]+p2[isPt]-p4[isPt]) + 2*PO4[isPt]*(3*p1[isPt]+2*p2[isPt]+p3[isPt]))
        }

      dALKdH.H2S <- rep(0.0, Nb_data)
      isHSt <- HSt!=0
        if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {
          dALKdH.H2S[isHSt] <- (-1/H[isHSt])*(-N*ALK.X + HS[isHSt])}  else {
            dALKdH.H2S[isHSt] <- (-1/H[isHSt])*(HS*s1[isHSt])
        } 
      
      dALKdH.SiOH4 <- rep(0.0, Nb_data)
      isSit <- Sit!=0
        if(X.name == "SiOH4" || X.name == "SiOOH3" || 
           X.name == "SiO2OH2" || X.name == "Sit") {
          dALKdH.SiOH4[isSit] <- (-1/H[isSit])*(-N*ALK.X + SiOOH3[isSit] + 4*SiO2OH2[isSit])} else {
            dALKdH.SiOH4[isSit] <- (-1/H[isSit])*(SiOOH3[isSit]*(si1[isSit]-si3[isSit]) + 
                                     2*SiO2OH2[isSit]*(2*si1[isSit]+si2[isSit]))
        } 
      
      dALKdH.BOH3 <- rep(0.0, Nb_data)
      isBOR <- BOR!=0
        if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
          dALKdH.BOH3[isBOR] <- (-1/H[isBOR])*(-N*ALK.X + BOH4[isBOR])} else {
            dALKdH.BOH3[isBOR] <-(-1/H[isBOR])*(BOH4[isBOR]*b1[isBOR])
        }

      dALKdH.F <- rep(0.0, Nb_data)
      isFLUO <- FLUO!=0
        if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
          dALKdH.F[isFLUO] <- (-1/H[isFLUO])*(-N*ALK.X + HF[isFLUO])} else {
            dALKdH.F[isFLUO] <- (-1/H[isFLUO])*(-HF[isFLUO]*f2[isFLUO])
        }

      dALKdH.SO4 <- rep(0.0, Nb_data)
      isST <- ST!=0
        if(X.name == "SO4" || X.name == "HSO4" || X.name == "ST") {
          dALKdH.SO4[isST] <- (-1/H[isST])*(-N*ALK.X + HSO4[isST])} else {
            dALKdH.SO4[isST] <- (-1/H[isST])*(-HSO4[isST]*so2[isST])
        }

      dALKdH.H <- (-1/H)*(OH+H)       # Internal enhancement of buffering
      return(dALKdH.CO2 + dALKdH.NH4 + dALKdH.H2PO4 + dALKdH.H2S + dALKdH.SiOH4 + dALKdH.BOH3 + 
               dALKdH.F + dALKdH.SO4 + dALKdH.H)
    }

    # Define function for determining totX
    totX.function <- function(X.name)
    {
      if(X.name == "CO2" || X.name == "HCO3" || X.name == "CO3" || 
         X.name == "DIC") {DIC} else
           if(X.name == "NH3" || X.name == "NH4" || X.name == "NH4t") {
             NH4t} else
               if(X.name == "H2PO4" || X.name == "H3PO4" || 
                  X.name == "HPO4" || X.name == "PO4" || 
                  X.name == "Pt") {Pt} else
                    if(X.name == "H2S" || X.name == "HS" || X.name == "HSt") {HSt} else
                      if(X.name == "SiOH4" || X.name == "SiOOH3" || X.name == "SiO2OH2" || 
                         X.name == "Sit") {Sit} else
                           if(X.name == "BOH3" || X.name == "BOH4" || X.name == "BOR") {
                             BOR} else
                               if(X.name == "F" || X.name == "HF" || X.name == "FLUO") {
                                 FLUO} else
                                   if(X.name == "SO4" || X.name == "HSO4" || 
                                      X.name == "ST") {ST} else
                                        NULL
    }
    
    # Define function to calculate beta.H
    beta.H_func <- function()
    {
      dALKdH.CO2 <- rep(0.0, Nb_data)
      isDIC <- DIC!=0
      dALKdH.CO2[isDIC] <- (-1/H[isDIC])*(HCO3[isDIC]*(c1[isDIC]-c3[isDIC]) + 2*CO3[isDIC]*(2*c1[isDIC]+c2[isDIC]))  
      
      dALKdH.NH4 <- rep(0.0, Nb_data)
      isNH4t <- NH4t!=0
      dALKdH.NH4[isNH4t] <- (-1/H[isNH4t])*(NH3[isNH4t]*n1[isNH4t])
      
      dALKdH.H2PO4 <- rep(0.0, Nb_data)
      isPt <- Pt!=0
      dALKdH.H2PO4[isPt] <- (-1/H[isPt])*(-H3PO4[isPt]*(-p2[isPt]-2*p3[isPt]-3*p4[isPt]) + 
                                        HPO4[isPt]*(2*p1[isPt]+p2[isPt]-p4[isPt]) + 2*PO4[isPt]*(3*p1[isPt]+2*p2[isPt]+p3[isPt]))
      
      dALKdH.H2S <- rep(0.0, Nb_data)
      isHSt <- HSt!=0
      dALKdH.H2S[isHSt] <- (-1/H[isHSt])*(HS[isHSt]*s1[isHSt])
      
      dALKdH.SiOH4 <- rep(0.0, Nb_data)
      isSit <- Sit!=0
      dALKdH.SiOH4[isSit] <- (-1/H[isSit])*(SiOOH3[isSit]*(si1[isSit]-si3[isSit]) + 
                                         2*SiO2OH2[isSit]*(2*si1[isSit]+si2[isSit]))
      
      dALKdH.BOH3 <- rep(0.0, Nb_data)
      isBOR <- BOR!=0
      dALKdH.BOH3[isBOR] <-(-1/H[isBOR])*(BOH4[isBOR]*b1[isBOR])
      
      dALKdH.F <- rep(0.0, Nb_data)
      isFLUO <- FLUO!=0
      dALKdH.F[isFLUO] <- (-1/H[isFLUO])*(-HF[isFLUO]*f2[isFLUO])
      
      dALKdH.SO4 <- rep(0.0, Nb_data)
      isST <- ST!=0
      dALKdH.SO4[isST] <- (-1/H[isST])*(-HSO4[isST]*so2[isST])
      
      dALKdH.H <- (-1/H)*(OH+H)       # Internal enhancement of buffering
      return(dALKdH.CO2 + dALKdH.NH4 + dALKdH.H2PO4 + dALKdH.H2S + dALKdH.SiOH4 + dALKdH.BOH3 + 
               dALKdH.F + dALKdH.SO4 + dALKdH.H)
    }
    
    # Define a function to calculate the Revelle factor 
    Revelle_func <- function()
    {
      X.name   <- "CO2"
      X        <- get(X.name)
      N        <- N.function(X.name)
      ALK.X    <- ALK.species(X.name)
      totX     <- totX.function(X.name)
      
      dALK.dH   <- dALKdH.function(X.name)
      dX.dtotX <- X * H * dALK.dH / (ALK.X^2 + (H*dALK.dH - N*ALK.X)*totX)
      RF <- as.numeric(DIC / CO2 * dX.dtotX)
      
      return(RF)
    }
    
    # Loop that calculates output species per variable
    for (i in 1:length(species))
    {
      
      # In case a change in the total concentration is specified, 
      # change X to the reference species for ALK
      X.name <- species[i]
      
      if(X.name == "DIC") X <- get("CO2") else                      
        if(X.name == "NH4t") X <- get("NH4") else                     
          if(X.name == "Pt") X <- get("H2PO4") else
            if(X.name == "HSt") X <- get("H2S") else 
              if(X.name == "Sit") X <- get("SiOH4") else
                if(X.name == "BOR") X <- get("BOH3") else
                  if(X.name == "FLUO") X <- get("F") else
                    if(X.name == "ST") X <- get("SO4") else
                      X <- get(species[i])
                        
      N              <- N.function(X.name)
      ALK.X          <- ALK.species(X.name)
      totX           <- totX.function(X.name)
      
      isTotX      <- totX!=0
      isNotSTFLUO <- !X.name %in% c("HF","F","FLUO","HSO4","SO4","ST")
      # Calculation of factors in Table 3a
      dALK.dH[[i]]     <- dALKdH.function(X.name)
      dtotX.dH[[i]][isTotX]    <- (N*totX[isTotX] - ALK.X[isTotX]) / H[isTotX]
      dALK.dX[[i]][isTotX]     <- ALK.X[isTotX] / X[isTotX]
      dtotX.dX[[i]][isTotX]    <- totX[isTotX] / X[isTotX]
      
      # Transformation from dH to dpH
      dH.dpH         <- -log(10)*H
      
      dALK.dpH[[i]]    <- dH.dpH * dALK.dH[[i]]
      dtotX.dpH[[i]]   <- dH.dpH * dtotX.dH[[i]]
      
      # Calculation of factors in Table 3b
      denom <- ALK.X^2 + (H*dALK.dH[[i]] - N*ALK.X)*totX
      comput  <- totX*H  / denom
      dH.dALK[[i]][isTotX][isNotSTFLUO]     <- comput[isTotX][isNotSTFLUO]
      comput  <- -ALK.X*H / denom
      dH.dtotX[[i]][isTotX]     <- comput[isTotX]
      comput  <- -X * (N*totX - ALK.X)  / denom
      dX.dALK[[i]][isTotX]      <- comput[isTotX]
      comput  <- X * H * dALK.dH[[i]] / denom
      dX.dtotX[[i]][isTotX]     <- comput[isTotX]
      
      # Transformation from dH to dpH
      dpH.dH         <- 1/(-log(10)*H)
      
      dpH.dALK[[i]]    <- dpH.dH * dH.dALK[[i]]
      dpH.dtotX[[i]]   <- dpH.dH * dH.dtotX[[i]]
    }  
    
    # Calculate the Revelle factor
    RF <- Revelle_func()
    
    # Calculate beta.H
    beta.H <- as.numeric(beta.H_func())
      
    # Combine the results in a list
    result <- list(Carbfull=as.matrix(Carbfull), dALK.dH=dALK.dH, dtotX.dH=dtotX.dH, dALK.dX=dALK.dX, 
                   dtotX.dX=dtotX.dX, dALK.dpH=dALK.dpH, dtotX.dpH=dtotX.dpH, dH.dALK=dH.dALK, 
                   dH.dtotX=dH.dtotX, dX.dALK=dX.dALK, dX.dtotX=dX.dtotX, dpH.dALK=dpH.dALK, 
                   dpH.dtotX=dpH.dtotX, beta.H=beta.H, RF=RF)
    
    # Return result
    return(result)
  })
}
  
  # Define all species for which bufferfunc needs to be calculated
  species <- names(Carbfull)[c(7,14:15,20:36,39:40,16,41:45)]
  
  output <- bufferfunc(Carbfull, species)
  return(output)
  
}
