# Cluster model equilibrium 
# Robert W. Buchkowski
# April 10, 2020

require(deSolve) # version 1.21

# How many replicates do you want to run in this code file?
NTOT = 50

# Load in the required functions ----
# A temperature function that takes day of the year (doy) and returns the temperature in Kelvin
LTtemp = function(doy){
  
  -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
  
}

# The complex model
singlemodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    TEMP = LTtemp(t %% 365) # Calculate temperature via day of the year
    
    # Model of earthworm growth. From ASA Johnston
    A_W = exp(-Ea/Kappa*(1/TEMP-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    A_P = exp(-B/TEMP)/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/TEMP)))
    
    # Modeled microbial dynamics MIMICS
    tempC = (TEMP-273.15)
    Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
    Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
    Klm = exp(Kslope*tempC + Kint)*Klm_mod
    Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
    
    dL = tp*P*P + (1-SUEh)*A_W*Vhp*H*P + th*H*H + tw*W*W - Vlm*L*M/(Klm + M) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - A_P*Vpf*N*P/(Kpf+N)
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = A_P*Vpf*N*P/(Kpf+N) - tp*P*P - A_W*Vhp*H*P
    
    dH = SUEh*A_W*Vhp*H*P - th*H*H
    
    dy = c(dP, dL, dM, dW, dN, dS, dH)
    
    return(list(c(dy)))
    
  }
  )
}

source("Scripts/parameters.R")

# Loop over different parameter sets ----
# A function that loops over different parameter sets
singlerun <- function(idx){
  ID = round(runif(1)*1e8,0) # set round ID
  paramscur = params # imported params
  yint2 = yint
  
  # Modify the parameter values
  for(i in 1:(length(paramscur)-2)){
    paramscur[i] = rlnorm(1, meanlog = log(params[i]), sdlog = 0.3536)
  }
  
  paramscur[25] = params[25] # Replace the original value of Kappa
  
  paramscur2 = paramscur
  
  yts = 2000 # Years to simulate
  ytsNE = 4 # Years to sumulate for non-equilibrium
  
  # Simulate each treatment to equilibrium and judge it stable if the change is less than 1e-4 for any state variable
  
  stablerunHW = ode(y=yint2,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur2)
  if(dim(stablerunHW)[1] == 365*yts){
    ystableHWsave2 = ystableHWsave = ystableHW = stablerunHW[(365*yts),-1]
    
    HWs = ifelse(max(abs(stablerunHW[(365*yts),-1] - stablerunHW[(365*(yts-1)),-1])) > 1e-4, 0,1)
    
    ystableHW["Stable"] = HWs
  }else{
    ystableHW = yint2
    ystableHW["Stable"] = 2
  }
   
  stablerunHW = ode(y=ystableHWsave,times = seq(1, 365*ytsNE,1), func=singlemodel, parms=paramscur2)
  
  
  # Each step after the first model with both animals sets their parameters or and popualtion to zero before running, then resets to default before the next run. E.g. ->
  paramscur2[c("Vlw", "Vsw")] = 0
  yint2["W"] = 0
  ystableHWsave["W"] = 0
  
  stablerunH = ode(y=yint2,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur2)
  if(dim(stablerunH)[1] == 365*yts){
    ystableH = stablerunH[(365*yts),-1]
    
    Hs = ifelse(max(abs(stablerunH[(365*yts),-1] - stablerunH[(365*(yts-1)),-1])) > 1e-4, 0,1)
    
    ystableH["Stable"] = Hs
  }else{
    ystableH = yint2
    ystableH["Stable"] = 2
  } 
  
  stablerunH = ode(y=ystableHWsave,times = seq(1, 365*ytsNE,1), func=singlemodel, parms=paramscur2)
  
  yint2 = yint
  ystableHWsave = ystableHWsave2
  paramscur2 = paramscur
  
  paramscur2[c("Vhp")] = 0
  yint2["H"] = 0
  ystableHWsave["H"] = 0
  
  stablerunW = ode(y=yint2,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur2)
  if(dim(stablerunW)[1] == 365*yts){
    ystableW = stablerunW[(365*yts),-1]
    
    Ws = ifelse(max(abs(stablerunW[(365*yts),-1] - stablerunW[(365*(yts-1)),-1])) > 1e-4, 0,1)
    
    ystableW["Stable"] = Ws
  }else{
    ystableW = yint2
    ystableW["Stable"] = 2
  } 
    
  stablerunW = ode(y=ystableHWsave,times = seq(1, 365*ytsNE,1), func=singlemodel, parms=paramscur2)
  
  paramscur2[c("Vlw", "Vsw")] = 0
  yint2["W"] = 0
  ystableHWsave["W"] = 0
  
  stablerunN = ode(y=yint2,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur2)
  if(dim(stablerunN)[1] == 365*yts){
    ystableNsave = ystableN = stablerunN[(365*yts),-1]
    
    Ns = ifelse(max(abs(stablerunN[(365*yts),-1] - stablerunN[(365*(yts-1)),-1])) > 1e-4, 0,1)
    
    ystableN["Stable"] = Ns
    
  }else{
    ystableN = yint2
    ystableN["Stable"] = 2
  }
  
  stablerunN = ode(y=ystableHWsave,times = seq(1, 365*ytsNE,1), func=singlemodel, parms=paramscur2)
  
  # Done the four year runs
  
  out3 = list(
    cbind(rbind(ystableHW,ystableH,ystableW,ystableN), ID = rep(ID, 4)), 
              c(paramscur, ID = ID),
              cbind(rbind(stablerunHW[365*ytsNE,-1],stablerunH[365*ytsNE,-1],stablerunW[365*ytsNE,-1],stablerunN[365*ytsNE,-1]), ID = rep(ID, 4), Treatment = c("HW", "H", "W", "N")))
  
  return(out3)
  
}

repseq = seq(1, NTOT, 1)

# Run the models: This takes a long time (days) on a single core!!
out1 = lapply(repseq, FUN=singlerun)

# Clean up the data and output
out2 = vector(mode = "list", length = NTOT)
out3 = vector(mode = "list", length = NTOT)
out4 = vector(mode = "list", length = NTOT)

for(jj in 1:NTOT){
  out2[[jj]] = out1[[jj]][[1]]
  
  out3[[jj]] = out1[[jj]][[2]]
  
  out4[[jj]] = out1[[jj]][[3]]
}

outf2 <- do.call("rbind", out2)
outf3 <- do.call("rbind", out3)
outf4 <- do.call("rbind", out4)

fname = paste0(round(runif(1), 7)*10000000,round(runif(1), 7)*10000000)

fname2 = paste0("[YOUR UNIQUE PATH TO SAVE THE DATA]/model_", fname, "_eqm.csv")

fname3 = paste0("[YOUR UNIQUE PATH TO SAVE THE DATA]/model_", fname, "_param.csv")

fname4 = paste0("[YOUR UNIQUE PATH TO SAVE THE DATA]/model_", fname, "_noneqm.csv")

write.csv(outf2, fname2, row.names = T)
write.csv(outf3, fname3, row.names = F)
write.csv(outf4, fname4, row.names = F)

# The compiled version of the data are saved as complex_model_10000_eqm.rds