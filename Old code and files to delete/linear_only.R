# Just a nice and simple linear model

# Model set up ----

NSLM <- function(t, y, pars){
  
  with(as.list(c(pars,y)),{
    
    # Temperature 
    Temp = -12.8244*cos(2*3.14/365*(t %% 365)-0.3666)+281.9846
    
    # Temperature-dependent functions (scaled to mean of 1)
    Aw = 1/(0.838083) * exp(-0.25/8.62e-05*(1/Temp-1/288.15)) # Based on Johnston et al. 2018 
    Ap = 1/(0.0001474841) * exp(-2493/Temp)/(1 + (2493/(26712-2493))*exp(26712*(1/(297.65) - 1/Temp))) # Based on Feng
    Am = 1/(484.8603) * exp(0.063*(Temp - 273.15) + 5.47) # Based on Wieder et al. 2013
    
    dH = Aw*khp0*P - klh*H
    dP = Ap*kpx0*X - Aw*khp0*P - klp*P
    dX = Ix - q*X - Ap*kpx0*X - ksx*X - Am*kmx0*X + kxw*W + kxs*S + kxm*M + kxl*L
    dW = Aw*kws0*S + Aw*kwm0*M + Aw*kwl0*L - kxw*W - klw*W
    dS = ksx*X + ksm*M + ksl*L - kxs*S - Aw*kws0*S - kms0*S
    dM = Am*kmx0*X + Am*kms0*S + Am*kml0*L - kxm*M - Aw*kwm0*M - ksm*M
    dL = klh*H + klp*P + klw*W - kxl*L - Aw*kwl0*L - ksl*L - kml0*L- q*L
    
    return(list(c(dH, dP, dX, dW, dS, dM, dL)))
    
  })
}

EQM <- function(pars){
  with(as.list(c(pars)),{
    
    # Temperature 
    Temp = mean(-12.8244*cos(2*3.14/365*(seq(1,365) %% 365)-0.3666)+281.9846)
    
    # Temperature-dependent functions (scaled to mean of 1)
    Aw = 1/(0.838083) * exp(-0.25/8.62e-05*(1/Temp-1/288.15)) # Based on Johnston et al. 2018 
    Ap = 1/(0.0001474841) * exp(-2493/Temp)/(1 + (2493/(26712-2493))*exp(26712*(1/(297.65) - 1/Temp))) # Based on Feng
    Am = 1/(484.8603) * exp(0.063*(Temp - 273.15) + 5.47) # Based on Wieder et al. 2013
    
    # calculate average parameters
    
    khp = Aw*khp0
    kpx = Ap*kpx0
    kmx = Am*kmx0
    kws = Aw*kws0
    kwm = Aw*kwm0
    kwl = Aw*kwl0
    kms = Am*kms0
    kml = Am*kml0
    
    output = c((Ix*khp*kpx*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + ksl + kwl + kxl + q) + 
                              klw*(kwm*kws*kxl + kml*kws*kxm + kws*kxl*kxm + ksl*kwm*kxs + kwm*kxl*kxs + kml*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + 
                                     (kwm + kxm)*(kws + kxs)*q + ksm*((kml + ksl)*kxs + (kws + kxs)*(kxl + q)) + kms*(kwm*(kxl + q) + kxm*(kml + ksl + kxl + q)))))/
                 (klh*(khp + klp)*q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*kpx*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + ksl + kwl + kxl + q) + 
                          klw*(kwm*kws*kxl + kml*kws*kxm + kws*kxl*kxm + ksl*kwm*kxs + kwm*kxl*kxs + kml*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + 
                                 (kwm + kxm)*(kws + kxs)*q + ksm*((kml + ksl)*kxs + (kws + kxs)*(kxl + q)) + kms*(kwm*(kxl + q) + kxm*(kml + ksl + kxl + q)))))/
                 ((khp + klp)*q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                                   klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                                          ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                                          kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                                          kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + ksl + kwl + kxl + q) + 
                      klw*(kwm*kws*kxl + kml*kws*kxm + kws*kxl*kxm + ksl*kwm*kxs + kwm*kxl*kxs + kml*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + 
                             (kwm + kxm)*(kws + kxs)*q + ksm*((kml + ksl)*kxs + (kws + kxs)*(kxl + q)) + kms*(kwm*(kxl + q) + kxm*(kml + ksl + kxl + q)))))/
                 (q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*(kmx*ksl*ksm*kws + kpx*ksl*ksm*kws + ksl*ksm*ksx*kws + kmx*ksm*kwl*kws + kpx*ksm*kwl*kws + ksm*ksx*kwl*kws + kmx*ksl*kwm*kws + 
                      kpx*ksl*kwm*kws + ksl*ksx*kwm*kws + kmx*kwl*kwm*kws + kpx*kwl*kwm*kws + ksx*kwl*kwm*kws + kmx*ksm*kws*kxl + ksm*ksx*kws*kxl + 
                      kmx*kwm*kws*kxl + ksx*kwm*kws*kxl + kpx*ksl*kws*kxm + ksl*ksx*kws*kxm + kpx*kwl*kws*kxm + ksx*kwl*kws*kxm + ksx*kws*kxl*kxm + 
                      kpx*ksm*kwl*kxs + kmx*ksl*kwm*kxs + kmx*kwl*kwm*kxs + kpx*kwl*kwm*kxs + kmx*kwm*kxl*kxs + kpx*kwl*kxm*kxs + 
                      kml*(kms*(kmx + kpx + ksx)*kwm + (kmx + kpx + ksx)*(ksm + kwm)*kws + ksx*kws*kxm + (kmx + kpx)*kwm*kxs) + 
                      ((kmx + ksx)*(ksm + kwm)*kws + ksx*kws*kxm + kmx*kwm*kxs)*q + 
                      kms*(kpx*(ksl + kwl)*kwm + kpx*kwl*kxm + kmx*kwm*(ksl + kwl + kxl + q) + ksx*kwm*(ksl + kwl + kxl + q))))/
                 (q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*(klw*(kml*ksm*(kmx + kpx + ksx) + kml*ksx*kxm + kmx*(ksl*(ksm + kwm) + ksm*(kxl + q)) + (ksm + kwm + kxm)*(kpx*ksl + ksx*(ksl + kxl + q))) + 
                      kxw*(kml*(kmx*ksm + kpx*ksm + ksx*(ksm + kwm + kxm)) + kmx*ksm*(ksl + kwl + kxl + q) + 
                             (ksm + kwm + kxm)*(kpx*ksl + ksx*(ksl + kwl + kxl + q)))))/
                 (q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*(klw*(kml*(kmx + kpx + ksx)*(kms + kws) + kml*(kmx + kpx)*kxs + kmx*(ksl*kxs + (kws + kxs)*(kxl + q)) + 
                           kms*(kpx*ksl + kmx*(ksl + kxl + q) + ksx*(ksl + kxl + q))) + 
                      kxw*(kml*(kms*(kmx + kpx + ksx) + (kmx + kpx)*(kws + kxs)) + kmx*(kws + kxs)*(ksl + kwl + kxl + q) + 
                             kms*(kpx*ksl + kmx*(ksl + kwl + kxl + q) + ksx*(ksl + kwl + kxl + q)))))/
                 (q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))),
               (Ix*(klw*(kms*(kmx + kpx + ksx)*kwm + kmx*(ksm + kwm)*kws + kms*kpx*kxm + kmx*kwm*kxs + (ksm + kwm + kxm)*((kpx + ksx)*kws + kpx*kxs)) + 
                      kpx*(kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw))/
                 (q*((kms*(kwm + kxm) + (ksm + kwm + kxm)*(kws + kxs))*kxw*(kml + kpx + ksl + kwl + kxl + q) + 
                       klw*(kpx*ksm*kws + ksm*ksx*kws + kpx*kwm*kws + ksx*kwm*kws + kmx*(ksm + kwm)*kws + ksm*kws*kxl + kwm*kws*kxl + kml*kws*kxm + kpx*kws*kxm + 
                              ksx*kws*kxm + kws*kxl*kxm + kml*ksm*kxs + kpx*ksm*kxs + ksl*ksm*kxs + kmx*kwm*kxs + kpx*kwm*kxs + ksl*kwm*kxs + ksm*kxl*kxs + 
                              kwm*kxl*kxs + kml*kxm*kxs + kpx*kxm*kxs + ksl*kxm*kxs + kxl*kxm*kxs + (ksm + kwm + kxm)*(kws + kxs)*q + 
                              kms*(kmx*kwm + ksx*kwm + kwm*kxl + (kml + ksl)*kxm + kxl*kxm + kpx*(kwm + kxm) + (kwm + kxm)*q)))))
    
    names(output) = c("H","P","X","W","S","M","L")
    return(output)
  })
  
  
}

params <-c(khp0 = 0.1,
           klh = 0.1,
           kpx0 = 0.1,
           klp = 0.1,
           Ix = 1,
           q = 0.1,
           ksx = 0.1,
           kmx0 = 0.1,
           kxw = 0.1,
           kxs = 0.1,
           kxm = 0.1,
           kxl = 0.1,
           kws0 = 0.1,
           kwm0 = 0.1,
           kwl0 = 0.1,
           klw = 0.1,
           ksm = 0.1,
           ksl = 0.1,
           kms0 = 0.1,
           kml0 = 0.1)

plot(ode(y = EQM(params), times = 1:365*5, func = NSLM, parms = params))

# Run a perliminary model fit to get the parameters in a reasonable range ----

EQMfit <- function(pars){
  PARS = exp(pars)
  return(sum(abs(c(H=0.0149,P=10.7,X=0.192,W=7.18,S=135,M=6.73,L=26) - 
            EQM(PARS))/
        c(H=0.0149,P=10.7,X=0.192,W=7.18,S=135,M=6.73,L=26)))
}

length(params)

N = 100
startpar = matrix(runif(N*20, min = 0, max = 1), nrow=N, ncol = 20)
endpar = startpar
fit = rep(0, N)
conv = rep(-1, N)

for(i in 1:N){
  params[] = log(startpar[i,])
  eqmfit = optim(params, fn = EQMfit, control=list(maxit = 10000))
  endpar[i,] = exp(eqmfit$par)
  fit[i] = eqmfit$value
  conv[i] = eqmfit$convergence
  if(i%%10 == 0) print(i)
}
  
colnames(endpar) = names(params)

hist(fit)

endpar %>%
  as_tibble() %>%
  mutate(Run = seq(1,100)) %>%
  gather(-Run, key= param,value = Value) %>%
  ggplot(aes(x=log(Value))) +
  geom_histogram() +
  facet_wrap(.~param) +
  theme_classic()

which(fit == min(fit))

EQM(endpar[68,])

saveRDS(endpar[68,],"Data/parameters.Rds")

# ... Functions for replicating the treatments ----

hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309, 314)+365, c(99, 280, 323, 114)+365*2)

samp = unique(c(hopsamp, soilsamp, wormsamp))

timestosample = samp[order(samp)]

rm(hopsamp, soilsamp, wormsamp,samp)

# Creates dataframes that can be used as events in ODE simulations

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2 = 1200 + 365

timelist2 = c(timelist2_expt, rep(seq(4,23, by=1), each=2)*365 + rep(c(121,305), 20))

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
                   value = c(0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, 40)),
                   method = c(rep("mult",3),rep("add", length(timelist2))))

rm(timelist2_expt,timelist2)

# Run simulations ----

guess = readRDS("Data/parameters.Rds")

datatomodel2 = readRDS("Data/datatomodel2.Rds")

Objective <- function(PAR){
  ystable = EQM(PAR)
  ok = 1
  while(ok > 1e-4){
    saf = ode(y = ystable, times = 1:100*365, 
              func = NSLM, parms = PAR)
    ystable = saf[99,-1]
    ok = sum(abs(saf[99,-1] - saf[98,-1]))
  }
  
  outRt = ode(y=ystable,times = 1:tmax2, 
              func=NSLM, 
              parms=PAR)
  
  COST = modCost(model = outRt,
                 obs = datatomodel2 %>%
                   filter(Treatment == "Rt") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std")
  
  outRmW = ode(y=ystable,times = 1:tmax2, 
               func=NSLM, parms=PAR,
               events = list(data=eshock_WE))
  
  COST = modCost(model = outRmW,
                 obs = datatomodel2 %>%
                   filter(Treatment == "RmW") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  
  outHW = ode(y=ystable,times = 1:tmax2, 
              func=NSLM, parms=PAR,
              events = list(data=eadd))
  
  COST = modCost(model = outHW,
                 obs = datatomodel2 %>%
                   filter(Treatment == "HW") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  
  outH = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR,
             events = list(data=eshock))
  
  COST = modCost(model = outH,
                 obs = datatomodel2 %>%
                   filter(Treatment == "H") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  
  ystable["H"] = 0
  PAR_noH = PAR
  
  PAR_noH[c("khp0", "klh")] = 0
  
  outW = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR_noH,
             events = list(data=eadd))
  
  COST = modCost(model = outW,
                 obs = datatomodel2 %>%
                   filter(Treatment == "W") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   filter(StateVar != "H") %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  
  outN = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR_noH,
             events = list(data=eshock))
  
  COST = modCost(model = outN,
                 obs = datatomodel2 %>%
                   filter(Treatment == "N") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   filter(StateVar != "H") %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  
  outRmH = ode(y=ystable,times = 1:tmax2, 
               func=NSLM, parms=PAR_noH)
  
  COST = modCost(model = outRmH,
                 obs = datatomodel2 %>%
                   filter(Treatment == "RmH") %>%
                   mutate(StateVar = ifelse(StateVar == "N", "X", StateVar)) %>%
                   filter(StateVar != "H") %>%
                   rename(name = StateVar, val = Biomass) %>%
                   select(name,time, val) %>%
                   data.frame(), y= "val", weight = "std", cost=COST)
  return(COST)
  
}

system.time(
  cost1 <- Objective(guess)
)

system.time(
  fit1 <- modFit(f = Objective, p = guess, lower = 0)
)

# Time: 1626.891

postfit <- function(PAR){
  ystable = EQM(PAR)
  ok = 1
  while(ok > 1e-4){
    saf = ode(y = ystable, times = 1:100*365, 
              func = NSLM, parms = PAR)
    ystable = saf[99,-1]
    ok = sum(abs(saf[99,-1] - saf[98,-1]))
  }
  
  outRt = ode(y=ystable,times = 1:tmax2, 
              func=NSLM, 
              parms=PAR)
  
  outRmW = ode(y=ystable,times = 1:tmax2, 
               func=NSLM, parms=PAR,
               events = list(data=eshock_WE))
  
  outHW = ode(y=ystable,times = 1:tmax2, 
              func=NSLM, parms=PAR,
              events = list(data=eadd))
  
  outH = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR,
             events = list(data=eshock))
  
  ystable["H"] = 0
  PAR_noH = PAR
  
  PAR_noH[c("khp0", "klh")] = 0
  
  outW = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR_noH,
             events = list(data=eadd))
  
  outN = ode(y=ystable,times = 1:tmax2, 
             func=NSLM, parms=PAR_noH,
             events = list(data=eshock))
  
  outRmH = ode(y=ystable,times = 1:tmax2, 
               func=NSLM, parms=PAR_noH)

  return(rbind(outRt, outRmW, outRmH, outHW, outH,
               outW,outN))
  
}

pf1 <-postfit(fit1$par)

pf1 <- as.data.frame(pf1)

pf1[,"Treatment"] = rep(c("Rt", "RmW", "RmH", "HW", "H",
                          "W","N"), each = dim(pf1)[1]/7)

datatomodel2 %>%
  mutate(Treatment = ifelse(Treatment =="RtH", "Rt", Treatment)) %>%
  mutate(StateVar = ifelse(StateVar =="N", "X", StateVar)) %>%
  group_by(time, Treatment, StateVar) %>%
  summarize(Expt = mean(Biomass)) %>%
  left_join(
    pf1 %>% as_tibble() %>%
      gather(-time, -Treatment, key = StateVar, value = Model)
  ) %>%
  ggplot(aes(x=Model, y=Expt, color=Treatment)) +
  geom_point() + facet_wrap(.~StateVar, scales="free")

pf1 %>% as_tibble() %>%
  gather(-time, -Treatment, key = StateVar, value = Model) %>%
  ggplot(aes(x=time, y=Model, color=Treatment)) +
  geom_line() + facet_wrap(.~StateVar, scales="free")

# Try a second guess to see how robust the fit is

guess2 = guess*(1 + runif(length(guess))) 

system.time(
  fit2 <- modFit(f = Objective, p = guess2, lower = 0)
)

# Time: 3168.886

round(100*(fit1$par-fit2$par)/fit1$par)

pf2 <-postfit(fit2$par)

pf2 <- as.data.frame(pf2)

pf2[,"Treatment"] = rep(c("Rt", "RmW", "RmH", "HW", "H",
                          "W","N"), each = dim(pf2)[1]/7)

datatomodel2 %>%
  mutate(Treatment = ifelse(Treatment =="RtH", "Rt", Treatment)) %>%
  mutate(StateVar = ifelse(StateVar =="N", "X", StateVar)) %>%
  group_by(time, Treatment, StateVar) %>%
  summarize(Expt = mean(Biomass)) %>%
  left_join(
    pf2 %>% as_tibble() %>%
      gather(-time, -Treatment, key = StateVar, value = Model)
  ) %>%
  ggplot(aes(x=Model, y=Expt, color=Treatment)) +
  geom_point() + facet_wrap(.~StateVar, scales="free")

pf2 %>% as_tibble() %>%
  gather(-time, -Treatment, key = StateVar, value = Model) %>%
  ggplot(aes(x=time, y=Model, color=Treatment)) +
  geom_line() + facet_wrap(.~StateVar, scales="free")


# Compare two models (not doing exactly the same things...)
pf1 %>% as_tibble() %>%
  gather(-time, -Treatment, key = StateVar, value = Model1) %>%
  left_join(
    pf2 %>% as_tibble() %>%
      gather(-time, -Treatment, key = StateVar, value = Model2)
  ) %>%
  mutate(time2 = time %% 100) %>%
  filter(time2 == 0) %>%
  ggplot(aes(x=Model1, y=Model2, color=Treatment)) +
  geom_point() + facet_wrap(.~StateVar, scales="free") +
  theme_classic()
