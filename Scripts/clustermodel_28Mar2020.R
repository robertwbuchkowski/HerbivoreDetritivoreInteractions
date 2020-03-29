# Model for cluster: Version 4.0

require(deSolve)

# How many replicates do you want to run in this code file?
NTOT = 100

# How many years do you want to simulate the experiments?
simyear = 100

# Load in the required functions ----
LTtemp = function(doy){
  
  -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
  
}

singlemodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    TEMP = LTtemp(t %% 365)
    
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

params<- c(Vlm_mod = 8e-6,
           Vsm_mod = 4e-06,
           Klm_mod = 0.143,
           Ksm_mod = 0.143,
           Vlw = 2.4e-06,#2.4e-06, #2.4e-05 before correction to Type I
           Vsw = 4.1e-05,#4.1e-05,#0.00462 before correction to Type I
           Vpf = 0.001/0.00018, #From JRS project 0.03
           Kpf = 0.08, #From JRS project 0.006
           Vhp = 0.01, #From JRS project 0.0025 - 0.0029
           SUEh = 0.7,
           SUE = 0.50,
           SUEws = 0.01,
           SUEwl = 0.02,
           SUEwm = 0.3,
           q = 0.1,
           IN= 0.02,
           l = 0.0001,
           tm = 0.01,
           tw = 0.00001,
           th = 20, # based on surivival from Schmitz lab experiments
           fi=0.6, #0.6,
           fo=0.001, #0.003,
           tp = 0.000005, #0.00008,
           Ea = 0.25,
           Kappa = 8.62e-05,
           B = 2493,
           D = 26712,
           Vslope = 0.063,
           Kslope = 0.007,
           Vint = 5.47,
           Kint = 3.19,
           Tref_W = 288.15,
           Tref_P = 297.65)

yint= c(P=95.238095, # WE with R:S ratio from Buchkowski et al. 2018
        L=25.526097, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        W=9.639477, # WE plots
        N=0.100000, # WE plots
        S=134.845515, # my historical data
        H=0.009) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11

# Single model sampling and treatment data frames --------

hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309, 314)+365, c(99, 280, 323, 114)+365*2)

samp = unique(c(hopsamp, soilsamp, wormsamp))

timestosample = samp[order(samp)]

rm(hopsamp, soilsamp, wormsamp,samp)

# Creates dataframes that can be used as events in ODE simulations

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

singlerun <- function(idx){
  ID = round(runif(1)*1e8,0) # set round ID
  ID2 = round(runif(1)*1e8,0) # set backup round ID
  paramscur = params # imported params
  ERRRRR = T # set loop
  ystable = yint
  
  while(ERRRRR){
    
    for(i in 1:(length(paramscur)-2)){
      paramscur[i] = rlnorm(1, meanlog = log(params[i]), sdlog = 0.3536)
    }
    
    yts = 2000
    
    stablerun = ode(y=ystable,times = seq(1, 365*yts,1), func=singlemodel, parms=paramscur)
    if(dim(stablerun)[1] == 365*yts){
      ystable = stablerun[(365*yts),-1]
      
      ERRRRR = max(abs(stablerun[(365*yts),-1] - stablerun[(365*(yts-1)),-1])) > 1e-4
    }

  }

  ystable2 = ystable
  
  output2_WE_Return = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur)
  
  output2_WE_Remove = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                          events = list(data=eshock_WE))
  
  output2_HW = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                   events = list(data=eadd))
  
  output2_H = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eshock))
  
  ystable["H"] = 0
  
  output2_W = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eadd))
  
  output2_0 = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur,
                  events = list(data=eshock))
  
  output2_WE_Hopper = ode(y=ystable,times = 1:tmax2, func=singlemodel, parms=paramscur)
  
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

out1 = lapply(repseq, FUN=singlerun)

outf <- do.call("rbind", out1)

fname = paste0("/gpfs/loomis/home.grace/fas/schmitz/rwb45/Model_parameter_reps/modelcluster_", round(runif(1), 7)*10000000,round(runif(1), 7)*10000000, ".csv")

write.csv(outf, fname, row.names = F)