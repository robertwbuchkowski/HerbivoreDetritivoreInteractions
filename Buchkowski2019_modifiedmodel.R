Buchkowski2019 <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    if(AAT == 1){
      Temp = 281.9906
    }else{
      Temp = -12.8244*cos(2*3.14/365*(t %% 365)-0.3666)+281.9846
    }
    
    # Temperature-dependent functions
    A_W = exp(-0.25/8.62e-05*(1/Temp-1/288.15)) # Based on Johnston et al. 2018 
    A_P = exp(-2493/Temp)/(1 + (2493/(26712-2493))*exp(26712*(1/(297.65) - 1/Temp))) # Based on Feng
    
    V0 = exp(0.063*(Temp - 273.15) + 5.47) # Based on Wieder et al. 2013
    K0 = exp(0.007*(Temp - 273.15) + 3.19) # Based on Wieder et al. 2013
    
    mfr = (V0*Vlm*L*M/(K0*Klm + M))/((V0*Vsm*S*M/(K0*Ksm + M)) + V0*Vlm*L*M/(Klm + M))
    
    dL = tp*P + tw*W*W + th*H*H + # litter inputs from dead bodies
      (1-a_hp)*A_W*Vhp*H*P/(Khp + N) - # sloppy herbivore feeding/ urine loss
      V0*Vlm*L*M/(K0*Klm + M) - A_W*Vlw*L*W/(Klw + L) - # litter consumed by worms and microbes
      q2*L
    
    dM = a_m*(V0*Vlm*L*M/(K0*Klm + M) + V0*Vsm*S*M/(K0*Ksm + M)) - # microbes eat litter and soil
      A_W*Vlw*mfr*M*W/(Klw + L) - A_W*Vsw*(1-mfr)*M*W/(Ksw + S) - # worms eat microbes in litter and soil
      tm*M # microbes die
    
    dW = a_wl*(A_W*Vlw*L*W/(Klw + L)) + # worms eat litter
      A_W*Vsw*S*W/(Ksw + S) + # worms eat soil (unassimlated soil stays put)
      A_W*Vlw*mfr*M*W/(Klw + L) + A_W*Vsw*(1-mfr)*M*W/(Ksw + S) - # worm eat microbes (unassimulated microbes survive gut passage)
      tw*W*W # worms die
    
    dN = IN - q*N - # system gains and losses
      fi*N + fo*S + # exchange with SOM and IORG
      (1-a_m)*(V0*Vlm*L*M/(K0*Klm + M) + V0*Vsm*S*M/(K0*Ksm + M)) - # microbial mineralization
      A_P*Vnp*N*P/(Knp + N) # plant uptake of nitrogen 
    
    dS = tm*M - # input from dead microbes
      V0*Vsm*S*M/(K0*Ksm + M) - # loss to microbes
      A_W*Vsw*S*W/(Ksw + S) + # loss to worms
      (1-a_wl)*(Vlw*L*W/(Klw + L)) + # unassimulated worm faeces
      fi*N - fo*S # exchange with IORG
    
    dP = A_P*Vnp*N*P/(Knp + N) - # nitrogen gain
      tp*P - # density-dependent death
      A_W*Vhp*H*P/(Khp + N) # herbivory
    
    dP2 = A_P*Vnp*1.5*N*P2/(Knp + N) - # nitrogen gain
      tp*P2 - # density-dependent death
      A_W*Vhp*2*H*P2/(Khp + N) # herbivory
    
    dH = a_hp*A_W*Vhp*H*P/(Khp + N) + # herbivory normal plant
      a_hp*A_W*Vhp*H*P2/(Khp + N) - # herbivory second plant
      th*H*H # death
    
    return(list(c(dL, dM, dW, dN, dS, dP, dP2, dH)))
    
  }
  )
}