require(FME)
require(lubridate)
require(tidyverse)

# How many years do you want to simulate the experiments?
simyear = 100
# 0.0 Create directory if necessary -------
if(!dir.exists(paste0("simplemodel_",Sys.Date()))){
  dir.create(paste0("simplemodel_",Sys.Date()))
}

# .. 0.1 Load in data and functions ----------------------------------------------------

LTtemp = function(doy){
  
  -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
  
}

datatomodel2 = read_csv("Data/datatomodel2.csv")

# 1.0 Single Model -----------------------------------------

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

  # P            L            M            W            N 
  # 34.609554901  8.570928208  5.429586815  9.495648285  0.160438844 
  # S            H 
  # 81.147450074  0.006371372 
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

# 2.0 Direct effects -----
# . 2.1 Direct effects #1 -----
singlemodel_direct <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    # Scaled so each effect changes parameter 5% based on expected population size
    tp_Hmod = tp - H * 1.852524e-05 + 2.5e-07 
    Vpf_Wmod = Vpf + W * 0.03 - 0.3
    
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
    
    dL = tp_Hmod*P*P + (1-SUEh)*A_W*Vhp*H*P + th*H*H + tw*W*W - Vlm*L*M/(Klm + M) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - A_P*Vpf_Wmod*N*P/(Kpf+N)
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = A_P*Vpf_Wmod*N*P/(Kpf+N) - tp_Hmod*P*P - A_W*Vhp*H*P
    
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

# P            L            M            W            N 
# 34.609554901  8.570928208  5.429586815  9.495648285  0.160438844 
# S            H 
# 81.147450074  0.006371372 

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel_direct, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .. 2.1.1 Simulations -------------------------------------------------

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Convert plant biomass to aboveground rather than belowground
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
directoutput = output
rm(output)

# . 2.2 Direct effects #2 -----

singlemodel_direct2 <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    # Scaled so each effect changes parameter 5% based on expected population size
    Vhp_Hmod = Vhp - H * 0.025
    Vpf_Wmod = Vpf + W * 0.03 - 0.3
    
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
    
    dL = tp*P*P + (1-SUEh)*A_W*Vhp_Hmod*H*P + th*H*H + tw*W*W - Vlm*L*M/(Klm + M) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - A_P*Vpf_Wmod*N*P/(Kpf+N)
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = A_P*Vpf_Wmod*N*P/(Kpf+N) - tp*P*P - A_W*Vhp_Hmod*H*P
    
    dH = SUEh*A_W*Vhp_Hmod*H*P - th*H*H
    
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

# P            L            M            W            N 
# 34.609554901  8.570928208  5.429586815  9.495648285  0.160438844 
# S            H 
# 81.147450074  0.006371372 

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel_direct2, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .. 2.2.1 Simulations -------------------------------------------------

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=singlemodel_direct2, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Convert plant biomass to aboveground rather than belowground
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
directoutput2 = output
rm(output)

# 3.0 Multiple Model -------

multiplemodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    TEMP = LTtemp(t %% 365)
    
    # Model of earthworm growth. From ASA Johnston
    A_W = exp(-Ea/Kappa*(1/TEMP-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    A_P = exp(-B/TEMP)/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/TEMP)))
    
    # Modeled Herbivore death
    # th = ifelse(TEMP > 273.15, th0, ifelse(H < Hmin, 0, th0*WF))
    # Vhp = ifelse(TEMP < 273.15, 0, Vhp0)
    
    # Modeled microbial dynamics MIMICS
    tempC = (TEMP-273.15)
    Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
    Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
    Klm = exp(Kslope*tempC + Kint)*Klm_mod
    Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
    
    dL = a21*P1*P2 + a12*P1*P2 + (tp1*P1*P1 + tp2*P2*P2) + (1-SUEh)*(A_W*Vhp1*H*P1 + A_W*Vhp2*H*P2) + th*H*H + tw*W*W - Vlm*L*M/(Klm + M) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - (A_P*Vpf1*N*P1/(Kpf1+N) + A_P*Vpf2*N*P2/(Kpf2+N))
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP1 = A_P*Vpf1*N*P1/(Kpf1+N) - tp1*P1*P1 - A_W*Vhp1*H*P1 - a12*P1*P2
    
    dP2 = A_P*Vpf2*N*P2/(Kpf2+N) - tp2*P2*P2 - A_W*Vhp2*H*P2 - a21*P1*P2
    
    dH = SUEh*(A_W*Vhp1*H*P1 + A_W*Vhp2*H*P2) - th*H*H

    dy = c(dP1, dP2, dL, dM, dW, dN, dS, dH)
    
    return(list(c(dy)))
    
  }
  )
}

# . 3.1 Multiple Model #1 ----

params<- c(Vlm_mod = 8e-6,
           Vsm_mod = 4e-06,
           Klm_mod = 0.143,
           Ksm_mod = 0.143,
           Vlw = 2.4e-06,#2.4e-06, #2.4e-05 before correction to Type I
           Vsw = 4.1e-05,#4.1e-05,#0.00462 before correction to Type I
           
           Vpf1 = 0.001/0.00018, #From JRS project 0.03
           Kpf1 = 0.08, #From JRS project 0.006
           Vhp1 = 0.01, #From JRS project 0.0025 - 0.0029
           
           Vpf2 = 1.1*0.001/0.00018, #From JRS project 0.03
           Kpf2 = 0.08, #From JRS project 0.006
           Vhp2 = 1.8*0.01, #From JRS project 0.0025 - 0.0029
           
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
           
           tp1 = 0.000005, #0.00008,
           
           tp2 = 0.000005, #0.00008,
           
           a21 = 0.000003,
           
           a12 = 0.0000045,
           
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

yint= c(P1=95.238095/2, # WE
        P2=95.238095/2, # WE
        L=25.526097, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        W=9.639477, # WE plots
        N=0.100000, # WE plots
        S=134.845515, # my historical data
        H=0.009) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,365), func=multiplemodel, parms=params)

max(initialrun[(dim(initialrun)[1]-1),-1] - initialrun[(dim(initialrun)[1]-2),-1])

(yint3 = initialrun[(dim(initialrun)[1]-1),-1])

# P            L            M            W            N 
# 34.609554901  8.570928208  5.429586815  9.495648285  0.160438844 
# S            H 
# 81.147450074  0.006371372 

if(F){
  # yts = 2000
  # A0 = ode(y=yint,times = seq(1,365*yts,1), func=multiplemodel, parms=params)
  # yint0 = A0[365*(yts-1),-1]
  yint0 = yint3
  
  yts = 20
  A1 = ode(y=yint0,times = seq(1,365*yts,1), func=multiplemodel, parms=params)
  yint0["H"] = 0; params["Vhp1"] = 0 ; params["Vhp2"] = 0 
  A2 = ode(y=yint0,times = seq(1,365*yts,1), func=multiplemodel, parms=params)
  
  par(mfrow=c(1,2)); Adiff = 100*(A2 - A1)/A1 ; plot(Adiff[,"P2"], type = "l", col = "black", ylab = "Removing herbivore % change"); points(Adiff[,"P1"], type = "l", col = "red"); rg = range(c(A1[,"P1"]/A1[,"P2"],A2[,"P1"]/A2[,"P2"])); legend("topleft", legend = c("P2", "P1"), col = c("black", "red"), lty = 1); plot(P1/P2~time, data = A1, type = "l", ylim = c(floor(rg[1]), ceiling(rg[2])), col = "blue"); points(P1/P2~time, data = A2, type = "l", col = "orange"); par(mfrow=c(1,1));legend("topright", legend = c("Herbivore", "No herbivore"), col = c("blue", "orange"), lty = 1); rm(rg)
  
  A1 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "Yes") %>%
    bind_rows(
      A2 %>% as.data.frame() %>% as_tibble() %>% mutate(Herb = "No")
    ) %>% 
    gather(-time, - Herb, key = StateVar, value = Biomass) %>%
    ggplot(aes(x= time, y= Biomass, color = Herb)) + geom_line() + 
    facet_wrap(.~StateVar, scale = "free") + theme_classic()
}

# .... 3.1.0.1 Multiple model sampling and treatment data frames --------------------------------------

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


eshock <- data.frame(var = c("P1","P2", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 4), timelist2),
                     value = c(0.1,0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+4))

eadd <- data.frame(var = c("P1","P2", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 4), timelist2),
                   value = c(0.1,0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, simyear*2)),
                   method = c(rep("mult",4),rep("add", length(timelist2))))

rm(timelist2_expt,timelist2)

if(F){
  # Test a similation
  
  output2_HW = ode(y=yint3,times = 1:tmax2_expt, func=multiplemodel, parms=params,
                   events = list(data=eadd))
  
  output2_HW %>% as.data.frame() %>% as_tibble() %>%
    gather(-time, key = StateVar, value = Biomass) %>%
    ggplot(aes(x= time, y= Biomass)) + geom_line() + 
    facet_wrap(.~StateVar, scale = "free") + theme_classic()
}

# .. 3.1.1 Simulations ----

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Joint plant biomass
output[,"P"] = output$P1 + output$P2

# Convert plant biomass to aboveground rather than belowground
output$P = output$P*0.1 # 10% of plant biomass is aboveground

# Save the single model output
multipleoutput = output

if(F){
  pdf(paste0("simplemodel_",Sys.Date(),"/Grahipcs_mult_",
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

rm(output)

# . 3.2 Multiple Model #2 ----

params<- c(Vlm_mod = 8e-6,
           Vsm_mod = 4e-06,
           Klm_mod = 0.143,
           Ksm_mod = 0.143,
           Vlw = 2.4e-06,#2.4e-06, #2.4e-05 before correction to Type I
           Vsw = 4.1e-05,#4.1e-05,#0.00462 before correction to Type I
           
           Vpf1 = 0.001/0.00018, #From JRS project 0.03
           Kpf1 = 0.08, #From JRS project 0.006
           Vhp1 = 0.01, #From JRS project 0.0025 - 0.0029
           tp1 = 5e-06, #0.00008,
           a12 = 3e-06,
           
           Vpf2 = 0.001/0.00018, #From JRS project 0.03
           Kpf2 = 0.5*0.08, #From JRS project 0.006
           Vhp2 = 1.8*0.01, #From JRS project 0.0025 - 0.0029
           tp2 = 5e-06, #0.00008,
           a21 = 3e-06,
           
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

yint= c(P1=95.238095/2, # WE
        P2=95.238095/2, # WE
        L=25.526097, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        W=9.639477, # WE plots
        N=0.100000, # WE plots
        S=134.845515, # my historical data
        H=0.009) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11


yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,365), func=multiplemodel, parms=params)

max(initialrun[(dim(initialrun)[1]-1),-1] - initialrun[(dim(initialrun)[1]-2),-1])

(yint3 = initialrun[(dim(initialrun)[1]-1),-1])

# P            L            M            W            N 
# 34.609554901  8.570928208  5.429586815  9.495648285  0.160438844 
# S            H 
# 81.147450074  0.006371372 
# .. 3.2.1 Simulations ----

output2_WE_Return = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params)

output2_WE_Remove = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                        events = list(data=eshock_WE))

output2_HW = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                 events = list(data=eadd))

output2_H = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eshock))

yint3["H"] = 0

output2_W = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eadd))

output2_0 = ode(y=yint3,times = 1:tmax2, func=multiplemodel, parms=params,
                events = list(data=eshock))

output = as.data.frame(rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0))

rm(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)

output["Treatment"] = as.factor(rep(c("Rt","RmW","HW","H","W","N"), each=tmax2))

# Joint plant biomass
output[,"P"] = output$P1 + output$P2

# Convert plant biomass to aboveground rather than belowground
output$P = output$P*0.1 # 10% of plant biomass is aboveground

# Save the single model output
multipleoutput2 = output

rm(output)

# 4.0 Compare the dynamics of single and multiple species outputs -----

datatomodel2a = datatomodel2 %>% 
  group_by(time, Treatment, StateVar) %>%
  summarize(std = sd(Biomass), Biomass = mean(Biomass)) %>%
  mutate(lower = Biomass - std, upper = Biomass + std) %>% 
  filter(Treatment %in% c("Rt","RmW","HW","H","W","N")) %>%
  ungroup() %>% left_join(
    data.frame(Treatment = c("Rt",  "RmW", "HW",  "H",   "W",   "N"),
               Treatment3 = c("eRt",  "fRmW", "dHW",  "bH",   "cW",   "aN"))
  )

singleoutput["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
singleoutput["Type"] = "Single"
singleoutput["P2"] = NA

directoutput["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
directoutput["Type"] = "Direct"
directoutput["P2"] = NA

directoutput2["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
directoutput2["Type"] = "Direct2"
directoutput2["P2"] = NA

multipleoutput["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
multipleoutput["Type"] = "Multiple"

multipleoutput2["Treatment2"] = rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2)
multipleoutput2["Type"] = "Multiple2"

# Write data files
write_rds(singleoutput, "Data/singleoutput.rds")
write_rds(directoutput, "Data/directoutput.rds")
write_rds(directoutput2, "Data/directoutput2.rds")
write_rds(multipleoutput, "Data/multipleoutput.rds")
write_rds(multipleoutput2, "Data/multipleoutput2.rds")

# Read data files
singleoutput = read_rds("Data/singleoutput.rds")
directoutput = read_rds("Data/directoutput.rds")
directoutput2 = read_rds("Data/directoutput2.rds")
multipleoutput = read_rds("Data/multipleoutput.rds")
multipleoutput2 = read_rds("Data/multipleoutput2.rds")

outputall = rbind(singleoutput, directoutput,directoutput2, multipleoutput,multipleoutput2) %>% left_join(
  data.frame(Treatment = c("Rt",  "RmW", "HW",  "H",   "W",   "N"),
             Treatment3 = c("eRt",  "fRmW", "dHW",  "bH",   "cW",   "aN"))
)

speciescomp = outputall %>% filter(Type %in% c("Multiple", "Multiple2")) %>%
  select(time, P1, P2, Treatment, Type)

# . 4.1 Plot the outputs into a single file for viewing ----

StateVar_names <- c(
  "H" = "Herbivore" ,
  "M" = "Microbial",
  "N" = "Inorganic N",
  "P" = "Plant",
  "W" = "Earthworm",
  "L" = "Litter",
  "S" = "Soil organic matter"
)

Treatment3_names <- c(
  "bH" = "Herbivores (Expt)" ,
  "cW" = "Earthworms (Expt)",
  "dHW" = "Both (Expt)",
  "fRmW" = "Remove worms",
  "eRt" = "Return worms",
  "aN" = "Neither (Expt)"
)

effect_names <- c(
  "aH" = "+ Herbivores (Expt)" ,
  "bW" = "+ Earthworms (Expt)",
  "cBoth" = "+ Both (Expt)",
  "dInteraction" = "Interaction (Expt)",
  "eWE" = "Returning worms"
)

toplot = as_tibble(outputall) %>% 
  select(-Treatment,-Treatment3) %>%
  rename(Treatment = Treatment2) %>%
  gather(-time, -Treatment, -Type, key = StateVar, value = Biomass) %>% 
  spread(key=Treatment, value=Biomass) %>% 
  mutate(eWE = T5-T4, cBoth = T3-T0, aH = T2-T0, bW = T1-T0) %>% 
  mutate(dInteraction = (T3-T2-T1+T0)/T0) %>% 
  select(time, StateVar, Type, eWE, cBoth, aH, bW,dInteraction) %>% 
  gather(-time, -StateVar, -Type, key=Treatment, value=Biomass)

pdf(paste0("simplemodel_",Sys.Date(),"/Grahipcs_compare_",
           round(100*(hour(Sys.time()) + (minute(Sys.time()))/60)),
           ".pdf"), width=9.5, height=7)

outputall %>% as_tibble() %>% filter(time < tmax2_expt) %>% select(-Treatment2,-Treatment, -P1, -P2) %>% gather(-Type, -time, -Treatment3, key = StateVar, value = Biomass) %>% ggplot() + geom_line(aes(x=time, y=Biomass, linetype = Type), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment3, scales="free_y",labeller=labeller(StateVar = StateVar_names, Treatment3 = Treatment3_names)) + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Return")) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)")+ geom_errorbar(data = datatomodel2a, aes(x=time, ymin = lower, ymax = upper, color=Treatment3)) + geom_point(data = datatomodel2a, aes(x=time, y= Biomass, color=Treatment3), size=2)

source("Scripts/plotcompare.R")

plotcompare(outputall, "Single")
plotcompare(outputall, "Direct")
plotcompare(outputall, "Direct2")
plotcompare(outputall, "Multiple")
plotcompare(outputall, "Multiple2")

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

toplot %>% filter(time %in% seq(1, 105*365, by = 365/10)) %>% 
  mutate(time = time/365) %>% 
  filter(!(StateVar %in% c("H", "W","P1", "P2"))) %>%
  ggplot() + geom_line(aes(x=time, y=Biomass, color = Type), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y",labeller=labeller(StateVar = StateVar_names, Treatment = effect_names)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (years)") + geom_hline(yintercept = 0, lty=2) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

toplot %>% filter(time %in% seq(1, 105*365, by = 365/10)) %>% mutate(time = time/365) %>% filter(!(StateVar %in% c("P1", "P2","H", "W"))) %>% filter(Treatment == "dInteraction") %>%
  ggplot() + geom_line(aes(x=time, y=Biomass, color = Type), alpha=0.5) + theme_classic() + facet_wrap(.~StateVar, scales="free",labeller=labeller(StateVar = StateVar_names)) + ylab("Interaction effect (proportion of control)") + xlab("Time (years)") + geom_hline(yintercept = 0, lty=2) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

speciescomp %>% as_tibble() %>% mutate(Prop = P1/(P1+P2))  %>% 
  filter(time %in% seq(1, 105*365, by = 365/10)) %>%
  mutate(Year = time/365) %>%
  filter(Treatment %in% c("N", "H", "W", "HW")) %>%
  mutate(LW = ifelse(Treatment %in% c("N", "H"), 2,1)) %>%
  ggplot(aes(x=Year, y = Prop, color = Treatment, size = LW, linetype = Type)) + geom_line() + theme_classic() + scale_color_manual(values=c("purple", "brown", "green", "blue"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)")) + ylab("Proportion of fast growing plant") + scale_size(guide = F, range = c(1,3)) + xlab("Time (years)")

dev.off()

png(paste0("simplemodel_",Sys.Date(),"/talk.png"), width=7, height=5, units = "in", res = 1000)
toplot %>% filter(time %in% seq(1, 105*365, by = 365/10)) %>% mutate(time = time/365) %>% filter(!(StateVar %in% c("P1", "P2","H", "W"))) %>% filter(Treatment == "dInteraction") %>% 
  left_join(
    tibble(Type = c("Single","Direct", "Direct2", "Multiple", "Multiple2"),
           Type2 = c("Complex model","Direct effect 1", "Direct effect 2", "Multiple species 1", "Multiple species 2"))
  ) %>%
  ggplot() + geom_line(aes(x=time, y=Biomass, color = Type2), alpha=0.5) + theme_classic() + facet_wrap(.~StateVar, scales="free",labeller=labeller(StateVar = StateVar_names)) + ylab("Interaction effect (proportion of control)") + xlab("Time (years)") + geom_hline(yintercept = 0, lty=2) + scale_color_manual(name = "Model", values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.box = "horizontal")
dev.off()

png(paste0("simplemodel_",Sys.Date(),"/model_expt.png"), width=7, height=5,units = "in", res = 1000)
plotcompare(outputall, "Single")
dev.off()


effect_names <- c(
  "aH" = "+ Herbivore" ,
  "bW" = "+ Detritivore"
)

png(paste0("simplemodel_",Sys.Date(),"/appendixS2F1.png"), width=8, height=5, units = "in", res = 1000)
toplot %>% filter(time %in% seq(1, 105*365, by = 365/10)) %>% 
  mutate(time = time/365) %>% 
  filter(!(StateVar %in% c("H", "W","P1", "P2"))) %>%
  left_join(
    tibble(Type = c("Single","Direct", "Direct2", "Multiple", "Multiple2"),
           Type2 = c("Complex model","Direct effect 1", "Direct effect 2", "Multiple species 1", "Multiple species 2"))
  ) %>%
  filter(Treatment %in% c("aH", "bW")) %>%
  ggplot() + geom_line(aes(x=time, y=Biomass, color = Type2), alpha=0.5) + theme_classic() + facet_wrap(Treatment~StateVar, scales="free_y",labeller=labeller(StateVar = StateVar_names, Treatment = effect_names), ncol =5, nrow = 2) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (years)") + geom_hline(yintercept = 0, lty=2) + scale_color_manual(name = "Model", values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + theme(legend.position = "top")
dev.off()

# 5.0 Single Model: Run for 1000 years -----------------------------------------

simyear = 1000

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

yts = 1000

initialrun = ode(y=yint,times = seq(1,365*yts,1), func=singlemodel, parms=params)

max(initialrun[(365*yts),-1] - initialrun[(365*yts-365),-1])

initialrun %>% as.data.frame() %>% as_tibble() %>%
  filter(time > 365*(yts-5)) %>%
  gather(-time, key = StateVar, value = Biomass) %>%
  ggplot(aes(x= time, y= Biomass)) + geom_line() + 
  facet_wrap(.~StateVar, scale = "free") + theme_classic()

yint3 = initialrun[365*(yts-1),-1]

# .... 5.0.0.1 Single model sampling and treatment data frames --------------------------------------

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
longtermoutput = output
write_rds(longtermoutput, "Data/longtermoutput.rds")
rm(output)

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
