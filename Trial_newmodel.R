# Trial

library(FME)
library(tidyverse)

trial <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    # Constants
    Ea = 0.25
    Kappa = 8.62e-05
    Tref_W = 288.15
    Tref_P = 297.65
    B = 2493
    D = 26712
    Vslope = 0.063
    Vint = 5.47
    
    # Model of earthworm growth. From ASA Johnston
    wgr = wgr_mod*exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    hgr = hgr_mod*exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    pgr = A_P_mod*exp(-B/LTtemp(t %% 365))/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/LTtemp(t %% 365))))
    
    # Microbial Growth
    Vmd = exp(Vslope*(LTtemp(t %% 365)-273.15) + Vint)*Vmd_mod
    Vmn = exp(Vslope*(LTtemp(t %% 365)-273.15) + Vint)*Vmn_mod
    
    # Overall Model Equations
    
    dH = eh*hgr*Vhp*H*P - th*H
    
    dP = pgr*P*min(ppn*Vpn*Ni,1) - hgr*Vhp*H*P - tp*P*P
    
    dW = eW*(wgr*Vwd*W*No + wgr*Vwm*W*M) - tw*W
    
    dNo = th*H + tp*P*P + tw*W + tm*M - wgr*Vwd*W*No - Vmd*M*No - l*No
    
    dM = eM*Vmd*M*No - wgr*Vwm*W*M - tm*M
    
    dNi = (1 - eM)*Vmd*M*No + (1 - eW)*(wgr*Vwd*W*No + wgr*Vwm*W*M) + (1 - eh)*
      hgr*Vhp*H*P - q*Ni + iN - pgr*P*min(ppn*Vpn*Ni,1)
    
    dy = c(dH, dP, dW, dNo, dM, dNi)
    
    LIM = ifelse(ppn*Vpn*Ni<1, 0,1) # 0 = N-limited, 1 = C-limited
    
    return(list(c(dy), hgr=hgr, pgr = pgr, LIM=LIM))
    
  }
  )
}

params<- c(Vmd_mod = 8e-6,
           Vmn_mod = 8e-4,
           Vwd = 2.4e-05,#2.4e-06, #2.4e-05 before correction to Type I
           Vwm = 2.4e-05,#2.4e-06, #2.4e-05 before correction to Type I
           A_P_mod = 639.0611,#0.008*639.0611,
           wgr_mod = 1,
           hgr_mod = 0.1,
           ppn = 1,
           Vpn = 1,
           pgr = 1,
           Vhp = 0.0025/1.2, #From JRS project 0.0025 - 0.0029
           eh = 0.7,
           eM = 0.50,
           eW = 0.02,
           q = 0.2,
           l = 0.01,
           iN= 0.05,
           tm = 0.05,
           tw = 0.000015,
           th = 0.01, # based on surivival from Schmitz lab experiments
           tp = 0.00001,
           fi=0.6, #0.6,
           fo=0.002 #0.003,
           )

yint= c(H=0.01, # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11
        P=95.238095, # WE with R:S ratio from Buchkowski et al. 2018
        W=9.639477, # WE plots
        No=25.526097+134.845515, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        Ni=0.100000 # WE plots
        ) 

(ystable = stode(y=yint, func=trial, parms=params)$y)           


plot(ode(y=yint, times = 1:1000, func=trial, parms=params))

# 
# Vmd = 2.2e-10
# Vmn = 2.2e-9
# Vwd = 2.4e-05#2.4e-06 #2.4e-05 before correction to Type I
# Vwm = 2.4e-05#2.4e-06 #2.4e-05 before correction to Type I
# A_P_mod = 639.0611#0.008*639.0611
# wgr_mod = 1
# hgr_mod = 0.1
# ppn = 1
# Vpn = 1
# pgr = 1
# Vhp = 0.0025/1.2 #From JRS project 0.0025 - 0.0029
# eh = 0.7
# eM = 0.50
# eW = 0.02
# q = 0.2
# l = 0.01
# iN= 0.05
# tm = 0.05
# tw = 0.000015
# th = 0.01 # based on surivival from Schmitz lab experiments
# tp = 0.0001
# wgr = 1
# pgr = 0.1
# hgr = 0.1
# 
# 
# 
# eh*hgr*Vhp*0.01*95.24 - th*0.01
# 
# pgr*95.24*min(ppn*Vpn*0.1,1) - hgr*Vhp*0.01*95.24 - tp*95.24*95.24
# 
# eW*(wgr*Vwd*9.64*160 + wgr*Vwm*9.64*7.16) - tw*9.64
# 
# th*0.01 + tp*95.24*95.24 + tw*9.64 + tm*7.16 - wgr*Vwd*9.64*160 - Vmd*7.16*160 - l*160
# 
# eM*Vmd*7.16*160 - wgr*Vwm*9.64*7.16 - tm*7.16
# 
# (1 - eM)*Vmd*7.16*160 + (1 - eW)*(wgr*Vwd*9.64*160 + wgr*Vwm*9.64*7.16) + (1 - eh)*
#   hgr*Vhp*0.01*95.24 - q*0.1 + iN - pgr*95.24*min(ppn*Vpn*0.1,1)


# Try just AG

trial <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    # Constants
    Ea = 0.25
    Kappa = 8.62e-05
    Tref_W = 288.15
    Tref_P = 297.65
    B = 2493
    D = 26712
    Vslope = 0.063
    Vint = 5.47
    
    # Model of earthworm growth. From ASA Johnston
    wgr = wgr_mod*exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    hgr = hgr_mod*exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    pgr = A_P_mod*exp(-B/LTtemp(t %% 365))/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/LTtemp(t %% 365))))
    
    tp2 = ifelse((LTtemp(t %% 365)-273.15)>0, tp, tp*1000)
    
    th2 = ifelse((LTtemp(t %% 365)-273.15)<0 & H > 0.01, th*100, th)
    
    # Overall Model Equations
    
    dH = eh*hgr*Vhp*H*P - th2*H
    
    dP = pgr*P*min(ppn*Vpn*Ni,1) - hgr*Vhp*H*P - tp2*P*P
    
    dNi = (1 - eh)*hgr*Vhp*H*P - q*Ni + iN - pgr*P*min(ppn*Vpn*Ni,1)
    
    dy = c(dH, dP, dNi)
    
    LIM = ifelse(ppn*Vpn*Ni<1, 0,1) # 0 = N-limited, 1 = C-limited
    
    return(list(c(dy), th2=th2, pgr = pgr, LIM=LIM))
    
  }
  )
}

params<- c(Vmd_mod = 8e-6,
           Vmn_mod = 8e-4,
           Vwd = 2.4e-05,#2.4e-06, #2.4e-05 before correction to Type I
           Vwm = 2.4e-05,#2.4e-06, #2.4e-05 before correction to Type I
           A_P_mod = 639.0611,
           wgr_mod = 1,
           hgr_mod = 0.1,
           ppn = 0.5,
           Vpn = 0.05,
           pgr = 1,
           Vhp = 0.1*0.025, #From JRS project 0.0025 - 0.0029
           eh = 0.7,
           eM = 0.50,
           eW = 0.02,
           q = 0.1,
           l = 0.01,
           iN= 0.05,
           tm = 0.05,
           tw = 0.000015,
           th = 0.001, # based on surivival from Schmitz lab experiments
           tp = 0.00000001,
           fi=0.6, #0.6,
           fo=0.002 #0.003,
)

yint= c(H=0.01,#.01, # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11
        P=95.238095, # WE with R:S ratio from Buchkowski et al. 2018
        Ni=0.100000 # WE plots
) 

plot(ode(y=yint, times = 1:1000, func=trial, parms=params))
yint[1] = 0
plot(ode(y=yint, times = 1:1000, func=trial, parms=params))

# why not try a smooth annual growth function that relies on the results.

## Try just AG once per year

trial <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    dH = eh*Vhp*H*P - th*H
    
    dP = pgr*P - Vhp*H*P - tp*P*P
    
    return(list(c(dH, dP)))
    
  }
  )
}

stode(y = c(H = 0.01, P = 10), func = trial, parms = exp(out$par))

params= c(
          Vhp = 0.00025,
          th = 0.01,
          pgr = 75,
          tp = 0.01)

eh = 0.7

fitfcn <- function(params){
  with(as.list(exp(params)),{
    aa = (th/(eh*Vhp) - 10)^2 + 
      (pgr/tp - 15)^2 +
      100*(-((th*tp - eh*pgr*Vhp)/(eh*Vhp^2))-0.01)^2
    
    return(aa)
  })
}
fitfcn(log(params))

out = optim(log(params), fitfcn)

exp(out$par)

with(as.list(exp(out$par)),{
 c(th/(eh*Vhp), pgr/tp, (-((th*tp - eh*pgr*Vhp)/(eh*Vhp^2))))
  # P2, P1, H2
 
})

# With DD herbivore death

params= c(
          Vhp = 0.00025,
          th = 0.01,
          pgr = 75,
          tp = 0.01)

eh = 0.7

fitfcn <- function(params){
  with(as.list(exp(params)),{
    aa = ((pgr*th)/(th*tp + eh*Vhp^2) - 10)^2 + 
      (pgr/tp - 15)^2 +
      1000*((eh*pgr*Vhp)/(th*tp + eh*Vhp^2)-0.01)^2
    
    return(aa)
  })
}
fitfcn(log(params))

out = optim(log(params), fitfcn)

exp(out$par)

with(as.list(exp(out$par)),{
  c(pgr/tp,
    (pgr*th)/(th*tp + eh*Vhp^2),
    (eh*pgr*Vhp)/(th*tp + eh*Vhp^2)
    )
  # P2, P1, H2
  
})

# Hard to get the behaviour right. When the herbivores don't have DD death they destroy the plants, when they do, they don't impact the plants...