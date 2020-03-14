# Just a nice and simple linear model

require(FME)
require(tidyverse)

# Model set up ----

NSLM <- function(t, y, pars){
  
  with(as.list(c(pars,y)),{
    
    # Temperature 
    if(AAT == 1){
      Temp = 281.9906
    }else{
      Temp = -12.8244*cos(2*3.14/365*(t %% 365)-0.3666)+281.9846
    }
    
    # Temperature-dependent functions (scaled to mean of 1)
    Aw = 1/(0.838083) * exp(-0.25/8.62e-05*(1/Temp-1/288.15)) # Based on Johnston et al. 2018 
    Ap = 1/(0.0001474841) * exp(-2493/Temp)/(1 + (2493/(26712-2493))*exp(26712*(1/(297.65) - 1/Temp))) # Based on Feng
    Am = 1/(484.8603) * exp(0.063*(Temp - 273.15) + 5.47) # Based on Wieder et al. 2013
    
    dH = eh*Aw*khp0*H*P - th*H*H
    dP = Ap*kpx0*P*X - tp*P - Aw*khp0*H*P
    dX = X0 - qx*X - Ap*kpx0*P*X + jm*M
    dM = Am*rm*M + Am*kml0*M*L - tm*M*M - qm*M - jm*M
    dW = Aw*rw*W + Aw*kwl0*W*L - tw*W*W
    dL = tw*W*W + tm*M*M + th*H*H + tp*P - Am*kml0*M*L - Aw*kwl0*W*L - ql*L
    
    return(list(c(dH, dP, dX, dW, dM, dL)))
    
  })
}

params <-c(eh = 0.5,
           khp0 = 0.1,
           kpx0 = 5,
           X0 = 1,
           qx = 0.1,
           jm = 0.1,
           kwl0 = 0.1,
           rw = 0.01,
           kml0 = 0.1,
           rm = 0.05,
           tm = 0.1,
           qm = 0.1,
           ql = 0.1,
           tp = 0.1,#1/0.0510,
           tw = 1/14.36,
           th = 1/0.0644,
           AAT = 1)

(RS = runsteady(y = y1, func = NSLM, parms = params)$y)

y4 = y3 = y2 = y1 = RS

y2["H"] = 0
y3["W"] = 0
y4[c("H", "W")] = 0

params4 = params3 = params2 = params

params2[c("khp0","th")] = 0
params3[c("kwl0", "rw", "tw")] = 0
params4[c("khp0", "kwl0", "rw", "tw", "th")] = 0

out1 = ode(y = y1, times = 1:365*5, func = NSLM, parms = params)
out2 = ode(y = y2, times = 1:365*5, func = NSLM, parms = params2)
out3 = ode(y = y3, times = 1:365*5, func = NSLM, parms = params3)
out4 = ode(y = y4, times = 1:365*5, func = NSLM, parms = params4)

rbind(out1, out2, out3, out4) %>%
  as_tibble() %>%
  mutate(Treatment = rep(c("HW", "W", "H", "N"), each= 365)) %>%
  gather(-time, -Treatment, key = StateVar, value=Model) %>%
  ggplot(aes(x=time, y=Model, color = Treatment)) +
  geom_line() + 
  facet_wrap(.~StateVar, scales = "free") + theme_classic()

rbind(runsteady(y = y1, func = NSLM, parms = params)$y,
      runsteady(y = y2, func = NSLM, parms = params2)$y,
      runsteady(y = y3, func = NSLM, parms = params3)$y,
      runsteady(y = y4, func = NSLM, parms = params4)$y)

sort(unique(datatomodel2$time))/365


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

datatomodel2 = readRDS("Data/datatomodel2.Rds")

Objective <- function(PAR){
  ystable = runsteady(y = y1, func = NSLM, parms = PAR)$y
  ok = 1
  while(ok > 1e-4){
    saf = ode(y = ystable, times = 1:100*365, 
              func = NSLM, parms = PAR)
    ystable = saf[99,-1]
    ok = sum(abs(saf[99,-1] - saf[98,-1]))
  }
  
  PAR["AAT"] = 0
  
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
  
  PAR_noH[c("khp0")] = 0
  
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

postfit <- function(PAR){
  ystable = runsteady(y = y1, func = NSLM, parms = PAR)$y
  ok = 1
  while(ok > 1e-4){
    saf = ode(y = ystable, times = 1:100*365, 
              func = NSLM, parms = PAR)
    ystable = saf[99,-1]
    ok = sum(abs(saf[99,-1] - saf[98,-1]))
  }
  
  PAR["AAT"] = 0
  
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
  
  PAR_noH[c("khp0")] = 0
  
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

Objective(PAR = params)

system.time(
  cost1 <- Objective(PAR)
)

system.time(
  fit1 <- modFit(f = Objective, p = params, lower = 0)
)

# Time: 3.23 hours

fit1$par

# eh        khp0          th        kpx0          tp          X0 
# 0.347167954 0.016594590 0.192679816 5.756071331 1.262005583 0.335001737 
# qx          jm        kwl0          rw          tw        kml0 
# 0.476678861 0.013110934 0.065783030 0.051693531 0.095565338 0.258821928 
# rm          tm          qm          ql         AAT 
# 0.003578889 0.031925778 0.062507329 0.163406645 1.000012440 

system.time(
  cost2 <- Objective(fit1$par)
)

cost2$residuals %>%
  as_tibble() %>%
  ggplot(aes(x = mod, y = obs)) +
  geom_point() + 
  facet_wrap(.~name, scales="free") + theme_classic()

# Almost no variation in the model output, but large empirical variation

cost2$residuals %>%
  as_tibble() %>%
  ggplot(aes(x = x, y = obs)) +
  geom_point() + 
  geom_point(aes(x = x, y = mod), col="red") + 
  facet_wrap(.~name, scales="free") + theme_classic()


pf1 <-postfit(fit1$par)

pf1 <- as.data.frame(pf1)

pf1[,"Treatment"] = rep(c("Rt", "RmW", "RmH", "HW", "H",
                          "W","N"), each = 1565)

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

# Time: 

(fit1$par-fit2$par)/fit1$par