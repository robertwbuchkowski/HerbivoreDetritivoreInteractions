# Model for the cluster

require(deSolve) # version 1.21

# How many replicates do you want to run in this code file?
NTOT = 100

# How many years do you want to simulate the experiments?
simyear = 100

# Load in the required functions ----

# A temperature function that takes day of the year (doy) and returns the temperature in Kelvin (calculated in extract_data_model_analysis.R)

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

# Complex model sampling and treatment data frames --------

# A calculation of when the sample the model output in order to match it with the empirical results. If you save all the model data from >100,000 simulations the data file is prohibitively large

hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309, 314)+365, c(99, 280, 323, 114)+365*2)

samp = unique(c(hopsamp, soilsamp, wormsamp))

timestosample = samp[order(samp)]

rm(hopsamp, soilsamp, wormsamp,samp)

# Creates dataframes that can be used as events in ODE simulations -----

# The events are used the remove earthworms, plant biomass, or litter at the start or throughout the simulated experiments

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2_expt = 1200 + 365

runlong = T
tmax2 = ifelse(runlong, tmax2_expt + 365*simyear,tmax2_expt)

timelist2 = c(timelist2_expt, rep(seq(4,simyear+3, by=1), each=2)*365 + rep(c(121,305), simyear))

eshock_WE <- data.frame(var = rep("W", length(timelist2)),
                        time =  timelist2,
                        value = rep(0.65, length(timelist2)),
                        method = rep("mult", length(timelist2)))


eshock <- data.frame(var = c("P", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 3), timelist2),
                     value = c(0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+3))

eadd <- data.frame(var = c("P", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 3), timelist2),
                   value = c(0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, simyear*2)),
                   method = c(rep("mult",3),rep("add", length(timelist2))))

rm(timelist2_expt,timelist2)

# Loop over different parameter sets ----

# A function that simulates the model for a randomly drawn parameter set

singlerun <- function(idx){
  ID = round(runif(1)*1e8,0) # set round ID
  ID2 = round(runif(1)*1e8,0) # set backup round ID
  paramscur = params # imported params
  ERRRRR = T # set loop
  ystable = yint
  
  while(ERRRRR){ # A loop that draws parameter sets until one is found that generates a stable equilibrium 
    
    for(i in 1:(length(paramscur)-2)){
      paramscur[i] = rlnorm(1, meanlog = log(params[i]), sdlog = 0.3536)
    }
    
    paramscur[25] = params[25] # Replace the original value of Kappa
    
    yts = 2000
    
    stablerun = ode(y=ystable,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur)
    if(dim(stablerun)[1] == 365*yts){
      ystable = stablerun[(365*yts),-1]
      
      ERRRRR = max(abs(stablerun[(365*yts),-1] - stablerun[(365*(yts-1)),-1])) > 1e-4
    }

  }

  ystable2 = ystable # Save the stable model output
  
  # Simulate each treatment in turn 
  
  # Returning worms to the 1-m^2 plots
  output2_WE_Return = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur)
  
  # Removing worms from the 1-m^2 plots, using an event
  output2_WE_Remove = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                          events = list(data=eshock_WE))
  
  # Adding worms to the experimental plots, using an event
  output2_HW = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                   events = list(data=eadd))
  
  # Removing worms to the experimental plots, using an event
  output2_H = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eshock))
  
  ystable["H"] = 0
  
  # Adding worms to the experimental plots, using an event
  output2_W = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eadd))
  
  # Removing worms to the experimental plots, using an event
  output2_0 = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eshock))
  
  # Simulating plots in old-fields without grasshoppers
  output2_WE_Hopper = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur)
  
  # Compile all the simulations
  # Note: WE_Return is listed twice because the same simulation is used to model the 1-m^2 earthworm plots (with natural herbivory) and the herbivore expeirments with herbivores added.
  out = rbind(output2_WE_Return, output2_WE_Hopper,
              output2_WE_Return,output2_WE_Remove,
              output2_HW, output2_H, 
              output2_W, output2_0)
  
  out1 = cbind(out, Treatment = rep(c(7,6,5,4,3,2,1,0), each=tmax2),
               Expt = rep(c(1,1,1,1,2,2,2,2), each=tmax2),
               Run = rep(ID, dim(out)[1]),
               Run2 = rep(ID2, dim(out)[1]))
  
  out2 = out1[out1[,"time"] %in% timestosample,]
  
  lg <- dim(out2)[1]
  
  out3 <- cbind(out2,
                PARS = c(unname(paramscur),i, rep(-2, lg -1 - length(paramscur))),
                YSTABLE = c(unname(ystable2), rep(-2, lg - length(ystable2))))
  
  return(out3)
  
}

repseq = seq(1, NTOT, 1)

# Run the analysis. This takes a long time on a single core (days for 100 samples)
out1 = lapply(repseq, FUN=singlerun)

outf <- do.call("rbind", out1)

fname = paste0("[YOUR UNIQUE PATH TO SAVE THE DATA]/modelcluster_", round(runif(1), 7)*10000000,round(runif(1), 7)*10000000, ".csv")

write.csv(outf, fname, row.names = F)

# The results of this analysis are saved in the file: "fullmodeloutput_2020-04-07.rds"