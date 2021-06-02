# More complex model simulations:
# Author: Robert W. Buchkowski
# Date created: Dec. 2020
# Date modified: June 1/2021

# This script compares the single model with and without reverse MM dynamics for the microbial pool. 

require(FME) # version 1.3.5
require(lubridate) # version 1.7.4
require(tidyverse) # version 1.2.1

# How many years do you want to simulate the experiments?
simyear = 100
# 0.0 Create directory if necessary -------
if(!dir.exists(paste0("simplemodel_",Sys.Date()))){
  dir.create(paste0("simplemodel_",Sys.Date()))
}

# .. 0.1 Load in data and functions ----------------------------------------------------

# A temperature function that takes day of the year (doy) and returns the temperature in Kelvin
LTtemp = function(doy){
  
  -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
  
}

# Load in the data from field experiments
datatomodel2 = read_csv("Data/datatomodel2.csv")

# ------- EXPANATION FOR MAJOR CODE CHUNKS ---------#
# Sections 1 to 3 below run the different model versions. Section 1 runs the complex model with a single plant species. Section 2 runs the Direct effects model. Section 3 runs the multiple species model. In Sections 2 and 3, there are two version of each model with slight differences in their parameterization. Sections 1.0.0.1 and 3.1.0.1 load in the sampling frame and the event functions to simulate treatments. These are analagous to the ones presented in the script complex_model_non-eqm_to_cluster.R

# These models can take a while to run. You can look only at the results by skipping to section 5 and loading in the outputs provided in the data files.

# The section 6 plots the model over a longer time frame so that dynamics can be observed after initial transience has past.

# 1.0 Single Model (i.e the Complex Model) -----------------------------------------

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

# Load in the parameters:
source("Scripts/parameters.R")

if(F){ 
  yts = 100
  A0 = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel, parms=params)
  (yint0 = A0[365*(yts-1),-1])
  
  yts = 20
  A1 = ode(y=yint0,times = seq(1,365*yts,1), func=singlemodel, parms=params)
  yint0["H"] = 0; params["Vhp"] = 0 
  A2 = ode(y=yint0,times = seq(1,365*yts,1), func=singlemodel, parms=params)
  
  Adiff = 100*(A2 - A1)/A1 ; plot(Adiff[,"P"], type = "l")
  
  A1 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "Yes") %>%
    bind_rows(
      A2 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "No")
    ) %>% 
    gather(-time, - Herb, key = StateVar, value = Biomass) %>%
    ggplot(aes(x= time, y= Biomass, color = Herb)) + geom_line() + 
    facet_wrap(.~StateVar, scale = "free") + theme_classic()
}

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .... 1.0.0.1 Single model sampling and treatment data frames --------------------------------------

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

# Test one example
if(F){
  output2_HW = ode(y=yint3,times = 1:tmax2_expt, func=singlemodel, parms=params,
                   events = list(data=eadd))
  
  plot(output2_HW)
}

# .. 1.0.1 Simulations -------------------------------------------------

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Convert plant biomass to aboveground rather than belowground
head(output)
output[,"P1"] = output$P
output$P = output$P*0.1 # 10% of plant biomass is aboveground

if(F){
  
  
  datatomodel2a = datatomodel2 %>% 
    group_by(time, Treatment, StateVar) %>%
    summarize(std = sd(Biomass), Biomass = mean(Biomass)) %>%
    mutate(lower = Biomass - std, upper = Biomass + std) %>% 
    filter(Treatment %in% c("Rt","RmW","HW","H","W","N")) %>%
    ungroup()
  
  pdf(paste0("simplemodel_",Sys.Date(),"/Grahipcs",
             round(100*(hour(Sys.time()) + (minute(Sys.time()))/60)),
             ".pdf"), width=7, height=7)
  
  output %>% filter(time < tmax2_expt) %>% gather(-time, -Treatment, key = StateVar, value = Biomass) %>% ggplot() + geom_line(aes(x=time, y=Biomass), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y")  + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Return")) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)")+ geom_errorbar(data = datatomodel2a, aes(x=time, ymin = lower, ymax = upper, color=Treatment)) + geom_point(data = datatomodel2a, aes(x=time, y= Biomass, color=Treatment), size=2)
  
  outdataCOMP = datatomodel2a %>% 
    left_join(as_tibble(output) %>% gather(-time, -Treatment, key = StateVar, value = Model))
  
  SVorder = unique(outdataCOMP$StateVar)
  
  vec = rep(NA, length(SVorder))
  for(i in 1:length(SVorder)){
    ms =lm(Biomass~Model + 0, data=outdataCOMP %>% filter(StateVar==SVorder[i]))
    mss = ifelse(ms$coefficients[1]>=0, "", "-")
    vec[i] = paste0("R[Full]^2 == ",mss,round(summary(ms)$adj.r.squared,2))
  }
  rm(ms, mss)
  vec2 = rep(NA, length(SVorder))
  for(i in 1:length(SVorder)){
    ms =lm(Biomass~Model + 0, data=outdataCOMP %>% 
             filter(StateVar==SVorder[i] & Treatment %in% c("N","H","W","HW")))
    mss = ifelse(ms$coefficients[1]>=0, "", "")
    vec2[i] = paste0("R[Expt]^2 == ",mss, round(summary(ms)$adj.r.squared,2))
  }
  
  
  maxes = outdataCOMP %>% 
    group_by(StateVar) %>%
    summarize(BMmax = max(max(Biomass + upper), max(Model)))
  
  ann_text <- tibble(SVorder, vec, vec2) %>%
    gather(-SVorder, key = type, value = R2) %>%
    rename(StateVar = SVorder) %>%
    left_join(maxes) %>%
    mutate(Model = BMmax*0.2) %>%
    mutate(Biomass = ifelse(type == "vec", BMmax*0.9, BMmax*0.8))
  
  rm(vec, vec2)
  
  ggplot(outdataCOMP, aes(x=Model, y=Biomass, color=Treatment)) + geom_abline(intercept=0, slope=1, lty=2) + geom_point() + geom_errorbar(data = outdataCOMP, aes(ymin = lower, ymax = upper, color=Treatment)) +  facet_wrap(c("StateVar"), scales="free") + theme_classic() + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Addition")) + geom_text(data = ann_text,label = ann_text$R2, color="black", parse=T) + geom_blank(data = maxes %>% mutate(Model = BMmax, Biomass = BMmax, Treatment = "N") %>% bind_rows(maxes %>%  mutate(Model = 0, Biomass = 0, Treatment = "N")))
  
  output2 = output
  
  output2["Treatment"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
  
  as_tibble(output2) %>% gather(-time, -Treatment, key = StateVar, value = Biomass) %>% spread(key=Treatment, value=Biomass) %>% mutate(WE = T5-T4, Both = T3-T0, H = T2-T0, W = T1-T0) %>% mutate(Interaction = Both-(H+W)) %>% select(time, StateVar, WE, Both, H, W,Interaction) %>% gather(-time, -StateVar, key=Treatment, value=Biomass) %>% ggplot() + geom_line(aes(x=time, y=Biomass), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2)
  
  dev.off()
}

# Save the single model output
singleoutput = output
rm(output)

# 2.0 Model with DD microbial death and forward MM ----

singlemodel_fwdMM <-function(t, y,pars){
  
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
    
    dL = tp*P*P + (1-SUEh)*A_W*Vhp*H*P + th*H*H + tw*W*W - Vlm*L*M/(Klm + L) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + L) + Vsm*S*M/(Ksm + L)) - tm*M*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + L) + Vsm*S*M/(Ksm + L)) - A_P*Vpf*N*P/(Kpf+N)
    
    dS = tm*M*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + L) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = A_P*Vpf*N*P/(Kpf+N) - tp*P*P - A_W*Vhp*H*P
    
    dH = SUEh*A_W*Vhp*H*P - th*H*H
    
    dy = c(dP, dL, dM, dW, dN, dS, dH)
    
    return(list(c(dy)))
    
  }
  )
}

# Load in the parameters:
source("Scripts/parameters.R")

# Adjust Parameters to match new functions
params["tm"] = params["tm"]/6.715
params["Klm_mod"] = params["Klm_mod"]*0.26
params["Ksm_mod"] = params["Ksm_mod"]*0.26

if(F){ 
  yts = 100
  A0 = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel_fwdMM, parms=params)
  (yint0 = A0[365*(yts-1),-1])
  
  yts = 20
  A1 = ode(y=yint0,times = seq(1,365*yts,1), func=singlemodel_fwdMM, parms=params)
  yint0["H"] = 0; params["Vhp"] = 0 
  A2 = ode(y=yint0,times = seq(1,365*yts,1), func=singlemodel_fwdMM, parms=params)
  
  Adiff = 100*(A2 - A1)/A1 ; plot(Adiff[,"P"], type = "l")
  
  A1 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "Yes") %>%
    bind_rows(
      A2 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "No")
    ) %>% 
    gather(-time, - Herb, key = StateVar, value = Biomass) %>%
    ggplot(aes(x= time, y= Biomass, color = Herb)) + geom_line() + 
    facet_wrap(.~StateVar, scale = "free") + theme_classic()
}

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel_fwdMM, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .. 2.0.1 Simulations -------------------------------------------------

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=singlemodel_fwdMM, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Convert plant biomass to aboveground rather than belowground
head(output)
output[,"P1"] = output$P
output$P = output$P*0.1 # 10% of plant biomass is aboveground

if(F){
  
  
  datatomodel2a = datatomodel2 %>% 
    group_by(time, Treatment, StateVar) %>%
    summarize(std = sd(Biomass), Biomass = mean(Biomass)) %>%
    mutate(lower = Biomass - std, upper = Biomass + std) %>% 
    filter(Treatment %in% c("Rt","RmW","HW","H","W","N")) %>%
    ungroup()
  
  pdf(paste0("simplemodel_",Sys.Date(),"/Grahipcs",
             round(100*(hour(Sys.time()) + (minute(Sys.time()))/60)),
             ".pdf"), width=7, height=7)
  
  output %>% filter(time < tmax2_expt) %>% gather(-time, -Treatment, key = StateVar, value = Biomass) %>% ggplot() + geom_line(aes(x=time, y=Biomass), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y")  + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Return")) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)")+ geom_errorbar(data = datatomodel2a, aes(x=time, ymin = lower, ymax = upper, color=Treatment)) + geom_point(data = datatomodel2a, aes(x=time, y= Biomass, color=Treatment), size=2)
  
  outdataCOMP = datatomodel2a %>% 
    left_join(as_tibble(output) %>% gather(-time, -Treatment, key = StateVar, value = Model))
  
  SVorder = unique(outdataCOMP$StateVar)
  
  vec = rep(NA, length(SVorder))
  for(i in 1:length(SVorder)){
    ms =lm(Biomass~Model + 0, data=outdataCOMP %>% filter(StateVar==SVorder[i]))
    mss = ifelse(ms$coefficients[1]>=0, "", "-")
    vec[i] = paste0("R[Full]^2 == ",mss,round(summary(ms)$adj.r.squared,2))
  }
  rm(ms, mss)
  vec2 = rep(NA, length(SVorder))
  for(i in 1:length(SVorder)){
    ms =lm(Biomass~Model + 0, data=outdataCOMP %>% 
             filter(StateVar==SVorder[i] & Treatment %in% c("N","H","W","HW")))
    mss = ifelse(ms$coefficients[1]>=0, "", "")
    vec2[i] = paste0("R[Expt]^2 == ",mss, round(summary(ms)$adj.r.squared,2))
  }
  
  
  maxes = outdataCOMP %>% 
    group_by(StateVar) %>%
    summarize(BMmax = max(max(Biomass + upper), max(Model)))
  
  ann_text <- tibble(SVorder, vec, vec2) %>%
    gather(-SVorder, key = type, value = R2) %>%
    rename(StateVar = SVorder) %>%
    left_join(maxes) %>%
    mutate(Model = BMmax*0.2) %>%
    mutate(Biomass = ifelse(type == "vec", BMmax*0.9, BMmax*0.8))
  
  rm(vec, vec2)
  
  ggplot(outdataCOMP, aes(x=Model, y=Biomass, color=Treatment)) + geom_abline(intercept=0, slope=1, lty=2) + geom_point() + geom_errorbar(data = outdataCOMP, aes(ymin = lower, ymax = upper, color=Treatment)) +  facet_wrap(c("StateVar"), scales="free") + theme_classic() + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Addition")) + geom_text(data = ann_text,label = ann_text$R2, color="black", parse=T) + geom_blank(data = maxes %>% mutate(Model = BMmax, Biomass = BMmax, Treatment = "N") %>% bind_rows(maxes %>%  mutate(Model = 0, Biomass = 0, Treatment = "N")))
  
  output2 = output
  
  output2["Treatment"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
  
  as_tibble(output2) %>% gather(-time, -Treatment, key = StateVar, value = Biomass) %>% spread(key=Treatment, value=Biomass) %>% mutate(WE = T5-T4, Both = T3-T0, H = T2-T0, W = T1-T0) %>% mutate(Interaction = Both-(H+W)) %>% select(time, StateVar, WE, Both, H, W,Interaction) %>% gather(-time, -StateVar, key=Treatment, value=Biomass) %>% ggplot() + geom_line(aes(x=time, y=Biomass), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2)
  
  dev.off()
}

# Save the single model output
forwardMMoutput = output
rm(output)

# 3.0 Compare the reverse and forward versions of the model -----



as_tibble(singleoutput) %>%
  filter(Treatment %in% c("N", "H", "W","HW") & time %in% seq(180, 38065, by = 365)) %>%
  select(- P1) %>%
  pivot_longer(P:H) %>%
  mutate(Type = "Rev") %>%
  bind_rows(
    as_tibble(forwardMMoutput)%>%
      filter(Treatment %in% c("N", "H", "W","HW") & time %in% seq(180, 38065, by = 365)) %>%
      select(- P1) %>%
      pivot_longer(P:H) %>%
      mutate(Type = "Fwd")
  ) %>%
  pivot_wider(names_from = Type) %>%
  ggplot(aes(x = Rev, y = Fwd, color = Treatment)) + geom_point(alpha = 0.5) + facet_wrap(.~name, scales = "free") + theme_classic() + scale_color_manual(breaks = c("N", "W", "H", "HW"), values=c("purple", "brown", "green", "blue"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)")) + ylab("Traditional density-dependent model") + xlab("Reverse Michaelis-Menten denisty-depednent model")

as_tibble(singleoutput) %>%
  filter(Treatment %in% c("N", "H", "W","HW") & time %in% seq(180, 38065, by = 365)) %>%
  select(- P1) %>%
  pivot_longer(P:H) %>%
  mutate(Type = "Rev") %>%
  bind_rows(
    as_tibble(forwardMMoutput)%>%
      filter(Treatment %in% c("N", "H", "W","HW") & time %in% seq(180, 38065, by = 365)) %>%
      select(- P1) %>%
      pivot_longer(P:H) %>%
      mutate(Type = "Fwd")
  ) %>%
  ggplot(aes(x = time, y = value, color = Treatment, linetype = Type)) + geom_line() + facet_wrap(.~name, scales = "free") + theme_classic() + scale_color_manual(breaks = c("N", "W", "H", "HW"), values=c("purple", "brown", "green", "blue"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)"))

# 3.0 Single Model: Run the model for a long time (1000 years) ------

simyear = 1000 # How long to run the model?

# Reload the model and parameters, for convenience
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
           Tref_W = 288.15,
           Tref_P = 297.65,
           B = 2493,
           D = 26712,
           Vslope = 0.063,
           Kslope = 0.007,
           Vint = 5.47,
           Kint = 3.19)

yint= c(P=95.238095, # WE with R:S ratio from Buchkowski et al. 2018
        L=25.526097, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        W=9.639477, # WE plots
        N=0.100000, # WE plots
        S=134.845515, # my historical data
        H=0.009) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11

yts = 1000 # Years to simulate to equilibrium 

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .... 3.0.0.1 Single model sampling and treatment data frames for 1000 years --------------------------------------

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

# .. 3.0.1 Simulations of the long-term simulation -----

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=singlemodel, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Convert plant biomass to aboveground rather than belowground
head(output)
output[,"P1"] = output$P
output$P = output$P*0.1 # 10% of plant biomass is aboveground

# Save the single model output
longtermoutput = output
write_rds(longtermoutput, "Data/longtermoutput.rds")
rm(output)

# .. 3.0.2 Plotting of the long-term simulation ----

longtermoutput = read_rds("Data/longtermoutput.rds")

longtermoutput["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)

toplot = as_tibble(longtermoutput) %>% 
  select(-Treatment) %>%
  rename(Treatment = Treatment2) %>%
  gather(-time, -Treatment, key = StateVar, value = Biomass) %>% 
  spread(key=Treatment, value=Biomass) %>% 
  mutate(eWE = (T5-T4)/T4, cBoth = (T3-T0)/T0, aH = (T2-T0)/T0, bW = (T1-T0)/T0) %>% 
  mutate(dInteraction = (T3-T2-T1+T0)/T0) %>% 
  select(time, StateVar, eWE, cBoth, aH, bW,dInteraction) %>% 
  gather(-time, -StateVar, key=Treatment, value=Biomass)

png(paste0("simplemodel_",Sys.Date(),"/longtime.png"), width=14, height=8, units = "in", res = 1000)
toplot %>% filter(time %in% seq(1, 1005*365, by = 365/10)) %>% 
  mutate(time = time/365) %>% 
  filter(!(StateVar %in% c("H", "W","P1", "P2"))) %>%
  filter(Treatment %in% c("aH", "bW", "dInteraction")) %>%
  ggplot() + 
  geom_line(aes(x=time, y=Biomass), alpha=0.5) + theme_classic() + 
  facet_wrap(Treatment~StateVar, scales="free_y",labeller=labeller(StateVar = StateVar_names, Treatment = effect_names), ncol = 5, nrow = 3) + ylab("Effect on Biomass (proporation of control)") + 
  xlab("Time (years)") + geom_hline(yintercept = 0, lty=2) + theme(legend.position = "top")
dev.off()
