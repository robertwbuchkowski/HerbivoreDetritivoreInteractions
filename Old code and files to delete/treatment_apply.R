# Code for applying treatments to ODE simulations
runlong = F


# Sampling and treatment data frames --------------------------------------

hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309, 314)+365, c(99, 280, 323, 114)+365*2)

samp = unique(c(hopsamp, soilsamp, wormsamp))

timestosample = samp[order(samp)]

rm(hopsamp, soilsamp, wormsamp,samp)

# Creates dataframes that can be used as events in ODE simulations

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2_expt = 1200 + 365

tmax2 = ifelse(runlong, tmax2_expt + 365*20,tmax2_expt)

timelist2 = c(timelist2_expt, rep(seq(4,23, by=1), each=2)*365 + rep(c(121,305), 20))

eshock_WE <- data.frame(var = rep("W", length(timelist2)),
                        time =  timelist2,
                        value = rep(0.65, length(timelist2)),
                        method = rep("mult", length(timelist2)))


eshock <- data.frame(var = c("R", "P", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 4), timelist2),
                     value = c(0.1, 0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+4))

eadd <- data.frame(var = c("R", "P", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 4), timelist2),
                   value = c(0.1, 0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, 40)),
                   method = c(rep("mult",4),rep("add", length(timelist2))))

rm(timelist2_expt,tmax2_expt,timelist2)

# Parameter set run ------------------------------------

params<- c(Vlm_mod = 5e-6,
           Vsm_mod = 4e-06,
           Klm_mod = 0.143,
           Ksm_mod = 0.143,
           Vlw = 2.4e-05,#2.4e-06, #2.4e-05 before correction to Type I
           Vsw = 4.6e-04,#4.1e-05,#0.00462 before correction to Type I
           l = 0,
           q = 0.1,
           IN= 0.015, # Crowther et al. 2015
           tm = 0.005,
           tw = 0.0002,
           th = 0.5, # based on surivival from Schmitz lab experiments [MOD]
           fi= 0.6, #0.6,
           fo=0.0003, #0.003,
           fus=0.001,
           # Plant relationships
           alphaR_mod = 1e-6,
           alphaA_mod = 0.0001,
           VN_mod = 3e-6, #From JRS project
           KN_mod = 0.5, #From JRS project
           A_P_mod = 0.005*639.0611,
           VH_mod = 0.1, #From JRS project 0.0025 - 0.0029
           StU = 0.1, # Guess
           fstable = 0.6, # Guess!
           fstableworm = 0.8, # Guess!
           # Substrate use efficiencies (conversion efficiencies)
           SUEh = 0.50,
           SUE = 0.50,
           SUEws = 0.02,
           SUEwl = 0.04,
           SUEwm = 0.6,
           # Assimulation efficiencies
           AEh = 0.25,
           AEw = 0.5,
           # break in the model
           Ea = 0.25,
           Kappa = 8.62e-05,
           Tref_W = 288.15,
           Tref_P = 297.65,
           B = 2493,
           D = 26712,
           th_winter = 10,
           alphaA_winter = 0.1,
           Hmin = 0.005,
           Vslope = 0.063,
           Kslope = 0.007,
           Vint = 5.47,
           Kint = 3.19,
           Nplant = 1)

# Add the model and default starting conditions -----

yint_base= c(R=95, # WE with R:S ratio from Buchkowski et al. 2018
        P = 0.04, # WE plots
        L=26, # WE plots with C:N ratio from files
        M=7.2, # WE plots
        W=9.6, # WE plots
        N=0.1, # WE plots
        S=108, # proportion using my historial data of 135
        U= 27, # proportion using my historial data of 135
        H=0.01) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11

if(params["Nplant"] ==1){
  yint = yint_base
}else{
  temp <- rep(yint_base[c("R", "P")]/params["Nplant"], each=params["Nplant"])
  names(temp) <- c(paste0("R",seq(1, params["Nplant"],1)),paste0("P",seq(1, params["Nplant"],1)))
  
  yint = c(temp, yint_base[!(names(yint_base) %in% c("R", "P"))])
  
  rm(yint_base)
}

# Functions for model and temperature forcing -----

# See load_climate_data for the parameterization of this function
LTtemp = function(doy){

  -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846

}

#Doesn't work for multiple plants yet
# ...need to find a way to build the vectors in the model

plantmodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{

    # Model of earthworm growth. From ASA Johnston
    A_W = exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    A_P = exp(-B/LTtemp(t %% 365))/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/LTtemp(t %% 365))))
    
    # Modeled microbial dynamics MIMICS
    tempC = (LTtemp(t %% 365)-273.15)
    Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
    Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
    Klm = exp(Kslope*tempC + Kint)*Klm_mod
    Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
    
    # Plant root nitrogen uptake
    VR = exp(Vslope*tempC + Vint)*rep(VN_mod, Nplant)
    
    KN = rep(KN_mod,Nplant)
    
    # Cut off at freezing (want night, so add 5C to freezing cut off)
    if(LTtemp(t %% 365) > (273.15 + 5) ){
      AP = A_P*rep(A_P_mod, Nplant)
      VH = rep(VH_mod, Nplant)
      alphaR = rep(alphaR_mod, Nplant)
      alphaA = rep(alphaA_mod, Nplant)
      th2 = th
    }else{
      AP = rep(0, Nplant) # no AG growth when freezing
      VH = rep(0, Nplant) # no herbivory when freezing
      alphaR = rep(0, Nplant) # no root death when freezing
      alphaA = rep(alphaA_winter, Nplant) # rapid aboveground senescence
      th2 = ifelse(H < Hmin, 0, th_winter) # rapid herbivore death to min pop
    }
    
    #Equations
    dL = sum(alphaR*R*R) + sum(alphaA*P) + 
      (1-SUEh)*sum(VH*H*P) + th2*H*H + 
      A_W*tw*W*W - Vlm*L*M/(Klm + M) - 
      A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) - 
      tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = AEw*SUEwl*A_W*Vlw*L*W + 
      AEw*SUEws*A_W*Vsw*S*W +
      AEw*SUEws*A_W*Vsw*U*W + 
      AEw*SUEwm*A_W*Vsw*W*M - 
      A_W*tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + 
      (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) - 
      sum(VR*N*R/(KN+N)) + 
      (1-AEh)*SUEh*sum(VH*H*P) +
      (1-AEw)*SUEwl*A_W*Vlw*L*W + 
      (1-AEw)*SUEws*A_W*Vsw*S*W +
      (1-AEw)*SUEws*A_W*Vsw*U*W +
      (1-AEw)*SUEwm*A_W*Vsw*W*M
    
    dS = fstable*tm*M - 
      StU*Vsm*S*M/(Ksm + M) + 
      fstableworm*(1-SUEwl)*A_W*Vlw*L*W +
      fstableworm*(1-SUEws)*A_W*Vsw*S*W -
      A_W*Vsw*S*W +
      fi*N - fo*S +
      fus*U
    
    dU = (1-fstable)*tm*M + 
      StU*Vsm*S*M/(Ksm + M) + 
      (1-fstableworm)*(1-SUEwl)*A_W*Vlw*L*W + 
      (1-fstableworm)*(1-SUEws)*A_W*Vsw*S*W - 
      SUEws*A_W*Vsw*U*W - 
      Vsm*U*M/(Ksm + M) -
      fus*U
    
    dR = sum(VR*N*R/(KN+N)) - sum(alphaR*R*R) - AP*R
    
    dP = AP*R - sum(VH*H*P) - sum(alphaA*P)
    
    dH = AEh*SUEh*sum(VH*H*P) - th2*H*H
    
    Nmin = IN + fo*S + 
      (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) + 
      (1-AEh)*SUEh*sum(VH*H*P) +
      (1-AEw)*SUEwl*A_W*Vlw*L*W + 
      (1-AEw)*SUEws*A_W*Vsw*S*W +
      (1-AEw)*SUEws*A_W*Vsw*U*W +
      (1-AEw)*SUEwm*A_W*Vsw*W*M
    
    GPP = sum(VR*N*R/(KN+N))
    
    SU = S/(S+U)
    
    dy = c(dR, dP, dL, dM, dW, dN, dS, dU, dH)
    
    return(list(dy, Nmin=Nmin, GPP=GPP, percentStable=SU))
    
  }
  )
}

