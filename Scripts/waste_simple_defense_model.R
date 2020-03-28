# WASTE CODE FROM SIMPLE_DEFERENCE_MODEL

# ..............Simulate MULTIPLE sets of Treatments -------------------------------------------------

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2_expt = 1200 + 365

tmax2 = tmax2_expt + 365*20

timelist2 = c(timelist2_expt, rep(seq(4,23, by=1), each=2)*365 + rep(c(121,305), 20))

eshock_WE <- data.frame(var = rep("W", length(timelist2)),
                        time =  timelist2,
                        value = rep(0.65, length(timelist2)),
                        method = rep("mult", length(timelist2)))


eshock <- data.frame(var = c("P", "PA", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 4), timelist2),
                     value = c(0.1, 0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+4))

eadd <- data.frame(var = c("P", "PA", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 4), timelist2),
                   value = c(0.1, 0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, 40)),
                   method = c(rep("mult",4),rep("add", length(timelist2))))

# varparams = expand.grid(alpha_A = params["alpha_A"]*seq(0.5, 1.5, length=3),alpha_R = params["alpha_R"]*seq(0.5, 1.5, length=3), Vpf = params["Vpf"]*seq(0.5, 1.5, length=3))

varparams = data.frame(alpha_A = params["alpha_A"]*c(1,1,1,1,1,0.5,1.5),
                       alpha_R = params["alpha_R"]*c(1,1,1,0.5,1.5,1,1), 
                       Vpf = params["Vpf"]*        c(1,0.5,1.5,1,1,1,1))

reps= dim(varparams)[1]

outputall = matrix(NA, ncol=9, nrow=tmax2*6*reps)

params_save = params

for(i in 1: reps){
  
  params = params_save
  params["alpha_A"] = varparams[i,"alpha_A"]
  params["alpha_R"] = varparams[i,"alpha_R"]
  params["Vpf"] = varparams[i,"Vpf"]
  
  
  initialrun = ode(y=yint,times = 1:(365*100), func=singlemodel, parms=params)
  print(paste("Done Ecosystem ",i," of ", reps, "spin-up"))
  
  yint3= yint2 = initialrun[365*100,-1]
  
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
  
  
  outputall[(((i-1)*tmax2*6+1):(i*tmax2*6)),] = rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)
  print(paste("Done Ecosystem ",i," of ", reps))
}

outputall = as.data.frame(outputall)
names(outputall) = c("time", names(ystable))
outputall["Treatment"] = as.factor(rep(c(5,4,3, 2, 1, 0), each=tmax2, length=tmax2*6*reps))
outputall["Run"] = rep(seq(1,reps), each=tmax2*6)

outputall = as_tibble(outputall)

outputPARAMS = as_tibble(varparams) %>% mutate(Run = seq(1, reps))

outputall %>% left_join(outputPARAMS) %>% write_csv(paste ("output_ALL_PARAMS_",format(Sys.time(), "%Y-%m-%d_%H%M"), ".csv", sep=""))

outputall <- read_csv("output_ALL_PARAMS_2019-01-25_1251.csv") %>% select(time:Run) %>% mutate(Treatment = as.factor(Treatment))

outputALL = as_tibble(outputall) %>% rename(R = P) %>% 
  rename(P = PA) %>%
  gather(-time, -Treatment, -Run, key=StateVar, value=Biomass)

outputALL_usevar = outputALL

jpeg(paste ("model_results_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

outputALL %>% filter(time < 1500) %>% ggplot() + geom_line(aes(x=time, y=Biomass,group=Run), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y")  + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Return")) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)")+ geom_errorbar(data = datatomodel2, aes(x=time, ymin = lower, ymax = upper, color=Treatment)) + geom_point(data = datatomodel2, aes(x=time, y= Biomass, color=Treatment), size=2)

dev.off()

outputALL %>% group_by(time, Treatment, StateVar) %>% rename(Model = Biomass) %>% summarise(lower_x = mean(Model) - sd(Model), upper_x = mean(Model) - sd(Model), Model = mean(Model)) %>% ungroup()

outdataCOMP = datatomodel3 %>% select(-time2, -SD) %>% 
  left_join(outputALL %>% group_by(time, Treatment, StateVar) %>% 
              rename(Model = Biomass) %>% summarise(lower_x = min(Model), upper_x = max(Model), Model = median(Model)) %>% 
              ungroup()
  )

vec = rep(NA, length(unique(outdataCOMP$StateVar)))
for(i in 1:length(unique(outdataCOMP$StateVar))){
  ms =lm(Biomass~Model + 0, data=outdataCOMP %>% filter(StateVar==unique(outdataCOMP$StateVar)[i]))
  mss = ifelse(ms$coefficients[1]>=0, "", "-")
  vec[i] = paste0("R[Full]^2 == ",mss,round(summary(ms)$adj.r.squared,2))
}
rm(ms, mss)
vec2 = rep(NA, length(unique(outdataCOMP$StateVar)))
for(i in 1:length(unique(outdataCOMP$StateVar))){
  ms =lm(Biomass~Model + 0, data=outdataCOMP %>% filter(StateVar==unique(outdataCOMP$StateVar)[i] & Treatment %in% c(1,2,3,4)))
  mss = ifelse(ms$coefficients[1]>=0, "", "")
  vec2[i] = paste0("R[Expt]^2 == ",mss, round(summary(ms)$adj.r.squared,2))
}


facet_names <- c(
  `H` = "Grashoppers (H)",
  `M` = "Soil Microbes (M)",
  `N` = "Inorganic Nitrogen (N)",
  `W` = "Earthworms (W)",
  `P` = "Plant (P)"
)

ann_text <- data.frame(StateVar = rep(unique(outdataCOMP$StateVar),2),
                       R2 = c(vec, vec2),
                       Model = rep(c(1,0.22,2,5,0.12), 2),
                       Biomass = c(c(9,0.35,12,17,0.048),c(9.3,0.35,12,17,0.048)*0.9))

rm(vec, vec2)

jpeg(paste ("model_compare_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

ggplot(outdataCOMP, aes(x=Model, y=Biomass, color=Treatment)) + geom_abline(intercept=0, slope=1, lty=2) + geom_point() + geom_errorbar(data = outdataCOMP, aes(ymin = lower, ymax = upper, color=Treatment)) + geom_errorbarh(data = outdataCOMP, aes(xmin = lower_x, xmax = upper_x, color=Treatment)) + facet_wrap(c("StateVar"), scales="free", labeller= as_labeller(facet_names)) + theme_classic() + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Addition")) + geom_text(data = ann_text,label = ann_text$R2, color="black", parse=T)

dev.off()

outputALL2 = outputALL

outputALL2["Treatment"] = as.factor(rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2, length=tmax2*6*reps))

outputPARAMS

outputALL2_0 = outputALL2 %>% spread(key=Treatment, value=Biomass) %>% mutate(WE = T5-T4, Both = T3-T0, H = T2-T0, W = T1-T0) %>% mutate(Interaction = Both-(H+W)) %>% select(time, Run, StateVar, WE, Both, H, W,Interaction) %>% gather(-time,-Run, -StateVar, key=Treatment, value=Biomass) %>% left_join(outputPARAMS) 

jpeg(paste ("model_longterm_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

ggplot(outputALL2_0) + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .))

dev.off()

facet_names2 <- c(
  `H` = "Add Grashoppers (H)",
  `M` = "Soil Microbes (M)",
  `N` = "Inorganic N (N)",
  `W` = "Add Earthworms (W)",
  `R` = "Plant Roots (R)",
  `S` = "Soil N (S)",
  `Interaction` = "Interaction effect",
  `P` = "Plant (P)"
)

jpeg(paste ("longterm_main_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("H", "W", "Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .))

dev.off()

outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .))

# ..............Simulate EARTHWORM x Grasshopper plant change Treatments -------------------------------------------------

singlemodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    # Modify plant parameters based on grasshopper and earthworm effects
    
    alpha_R_pre = 0.00012 - H * 8e-4
    Vpf_pre = 0.015 + W * 0.003
    
    # Model of earthworm growth. From ASA Johnston
    A_W = exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    A_P1 = A_P_mod*exp(-B/LTtemp(t %% 365))/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/LTtemp(t %% 365))))
    
    # no AG growth when freezing
    A_P = ifelse(LTtemp(t %% 365) > 273.15, A_P1, 0)
    
    # AG rapid death when freezing
    alpha_A2 = ifelse(LTtemp(t %% 365) > 273.15, alpha_A, alpha_A*WF) 
    
    # Plant root dynamics slow in winter
    alpha_R2 = ifelse(LTtemp(t %% 365) > 273.15, alpha_R_pre, alpha_R_pre/WF) 
    
    Vpf_R = ifelse(LTtemp(t %% 365) > 273.15, Vpf_pre, Vpf_pre/WF) 
    
    # Modeled Herbivore death
    th2 = ifelse(LTtemp(t %% 365) > 273.15, 
                 ifelse(H < Hmin, 0, th*WF), th)
    Vhp2 = ifelse(LTtemp(t %% 365) > 273.15, 0, Vhp)
    
    # Modeled microbial dynamics MIMICS
    tempC = (LTtemp(t %% 365)-273.15)
    Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
    Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
    Klm = exp(Kslope*tempC + Kint)*Klm_mod
    Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
    
    dL = alpha_A2*PA*PA + alpha_R2*P*P + (1-SUEh)*Vhp2*H*PA + th2*H+ A_W*tw*W*W -Vlm*L*M/(Klm + M) - A_W*Vlw*L*W
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - A_W*tw*W*W
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vpf_R*N*P/(Kpf+N)
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = Vpf_R*N*P/(Kpf+N) - alpha_R2*P*P - A_P*P
    
    dPA = A_P*P - Vhp2*H*PA - alpha_A2*PA*PA
    
    dH = SUEh*Vhp2*H*PA - th2*H
    
    dy = c(dPA, dP, dL, dM, dW, dN, dS, dH)
    
    return(list(c(dy), alpha_R_pre, Vpf_pre))
    
  }
  )
}

# initialrun = ode(y=yint,times = 1:(365*2), func=singlemodel, parms=params)
# 
# plot(initialrun)

(ystable = stode(y=yint, func=singlemodel, parms=params)$y)

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2_expt = 1200 + 365

tmax2 = tmax2_expt + 365*20

timelist2 = c(timelist2_expt, rep(seq(4,23, by=1), each=2)*365 + rep(c(121,305), 20))

eshock_WE <- data.frame(var = rep("W", length(timelist2)),
                        time =  timelist2,
                        value = rep(0.65, length(timelist2)),
                        method = rep("mult", length(timelist2)))


eshock <- data.frame(var = c("P", "PA", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 4), timelist2),
                     value = c(0.1, 0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+4))

eadd <- data.frame(var = c("P", "PA", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 4), timelist2),
                   value = c(0.1, 0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, 40)),
                   method = c(rep("mult",4),rep("add", length(timelist2))))


varparams = data.frame(alpha_A = params["alpha_A"]*c(1),
                       alpha_R = params["alpha_R"]*c(1), 
                       Vpf = params["Vpf"]*        c(1))

reps= dim(varparams)[1]

outputall = matrix(NA, ncol=11, nrow=tmax2*6*reps)

params_save = params

for(i in 1: reps){
  
  params = params_save
  params["alpha_A"] = varparams[i,"alpha_A"]
  params["alpha_R"] = varparams[i,"alpha_R"]
  params["Vpf"] = varparams[i,"Vpf"]
  
  
  initialrun = ode(y=yint,times = 1:(365*100), func=singlemodel, parms=params)
  
  print(paste("Done Ecosystem ",i," of ", reps, "spin-up"))
  
  yint3= yint2 = initialrun[365*100,2:9]
  
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
  
  
  outputall[(((i-1)*tmax2*6+1):(i*tmax2*6)),] = rbind(output2_WE_Return,output2_WE_Remove,output2_HW, output2_H, output2_W, output2_0)
  
  print(paste("Done Ecosystem ",i," of ", reps))
}

outputall = as.data.frame(outputall)
names(outputall) = c("time", names(ystable))
outputall["Treatment"] = as.factor(rep(c(5,4,3, 2, 1, 0), each=tmax2, length=tmax2*6*reps))
outputall["Run"] = rep(seq(1,reps), each=tmax2*6)

colnames(outputall)[10:11] = c("alpha_R_pre", "Vpf_pre")

outputall = as_tibble(outputall)

outputPARAMS = as_tibble(varparams) %>% mutate(Run = seq(1, reps))

outputall %>% left_join(outputPARAMS) %>% write_csv(paste ("interaction_ALL_PARAMS_",format(Sys.time(), "%Y-%m-%d_%H%M"), ".csv", sep=""))

outputALL = as_tibble(outputall) %>% rename(R = P) %>% 
  rename(P = PA) %>%
  gather(-time, -Treatment, -Run, key=StateVar, value=Biomass)

jpeg(paste ("model_results_interaction",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

outputALL %>% filter(time < 1500) %>% ggplot() + geom_line(aes(x=time, y=Biomass,group=Run), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y")  + scale_color_manual(values=c("purple", "brown", "green", "blue", "pink", "red"),labels=c("None (Expt)", "Worm (Expt)", "Hopper (Expt)", "Both (Expt)", "Removal", "Return")) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)")+ geom_errorbar(data = datatomodel2, aes(x=time, ymin = lower, ymax = upper, color=Treatment)) + geom_point(data = datatomodel2, aes(x=time, y= Biomass, color=Treatment), size=2)

dev.off()

outputALL2 = outputALL

outputALL2["Treatment"] = as.factor(rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2, length=tmax2*6*reps))

outputPARAMS

outputALL2_0 = outputALL2 %>% spread(key=Treatment, value=Biomass) %>% mutate(WE = T5-T4, Both = T3-T0, H = T2-T0, W = T1-T0) %>% mutate(Interaction = Both-(H+W)) %>% select(time, Run, StateVar, WE, Both, H, W,Interaction) %>% gather(-time,-Run, -StateVar, key=Treatment, value=Biomass) %>% left_join(outputPARAMS) 

jpeg(paste ("model_longterm_interaction",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

ggplot(outputALL2_0) + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .))

dev.off()

facet_names2 <- c(
  `H` = "Add Grashoppers (H)",
  `M` = "Soil Microbes (M)",
  `N` = "Inorganic N (N)",
  `W` = "Add Earthworms (W)",
  `R` = "Plant Roots (R)",
  `S` = "Soil N (S)",
  `Interaction` = "Interaction effect",
  `P` = "Plant (P)"
)

jpeg(paste ("longterm_main_interaction",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("H", "W", "Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .))

dev.off()


# ..............Plot varparm and internal param together --------------------------------

outputall <- read_csv("output_ALL_PARAMS_2019-02-28_0132.csv") %>% select(time:Run) %>% mutate(Treatment = as.factor(Treatment))

outputall_anvar <- read_csv("interaction_ALL_PARAMS_2019-02-28_0203.csv") %>% select(time:Run) %>% mutate(Treatment = as.factor(Treatment)) %>% select(-alpha_R_pre, -Vpf_pre)

outputall_anvar$Run = rep(8, dim(outputall_anvar)[1])

outputALL = outputall %>% bind_rows(outputall_anvar) %>% rename(R = P) %>% 
  rename(P = PA) %>%
  gather(-time, -Treatment, -Run, key=StateVar, value=Biomass)

outputALL2 = outputALL

outputALL2["Treatment"] = as.factor(rep(c("T5","T4","T3","T2","T1","T0"), each=tmax2, length=tmax2*6*8))

outputALL2_0 = outputALL2 %>% spread(key=Treatment, value=Biomass) %>% mutate(WE = T5-T4, Both = T3-T0, H = T2-T0, W = T1-T0) %>% mutate(Interaction = Both-(H+W)) %>% select(time, Run, StateVar, WE, Both, H, W,Interaction) %>% gather(-time,-Run, -StateVar, key=Treatment, value=Biomass)

jpeg(paste ("interactioneffect_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

ggpubr::ggarrange(
  outputALL2_0 %>% filter(time > 8500 & Treatment == "Interaction"& Run !=8) %>% ggplot(aes(x=as.character(Run), y=Biomass, color=StateVar)) + geom_boxplot() + theme_classic() + ylab("Interaction Effect") + xlab("") + scale_x_discrete(labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .)),
  outputALL2_0 %>% filter(time > 8500 & Treatment == "Interaction" & Run ==8) %>% ggplot(aes(x=as.character(Run), y=Biomass, color=StateVar)) + geom_boxplot() + theme_classic() + ylab("") + xlab("") + scale_x_discrete(labels=expression(f(G,H))),
  ncol=2, nrow=1, common.legend = T, legend="top", labels = "auto", widths = c(6,1))

dev.off()

jpeg(paste ("model_longterm_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=10, res=600, units="in")

dev.off()

ggplot(outputALL2_0) + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y") + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

dev.off()

facet_names2 <- c(
  `H` = "Add Grashoppers (H)",
  `M` = "Soil Microbes (M)",
  `N` = "Inorganic N (N)",
  `W` = "Add Earthworms (W)",
  `R` = "Plant Roots (R)",
  `S` = "Soil N (S)",
  `Interaction` = "Interaction effect",
  `P` = "Plant (P)"
)

jpeg(paste ("longterm_main_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("H", "W", "Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

dev.off()

outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

p1 = outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("H")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Biomass~(g[N]~m^-2))) + xlab("") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

p2 = outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("W")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab("") + xlab("Time (Days since start)") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

p3 = outputALL2_0 %>% filter(StateVar %in% c("M", "N", "P", "R","S") & Treatment %in% c("Interaction")) %>% ggplot() + geom_line(aes(x=time, y=Biomass, color=as.factor(Run)), alpha=0.5) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab("") + xlab("") + geom_hline(yintercept = 0, lty=2) + scale_colour_discrete(name=c("Parameter Set"),labels=expression(Base, V[PF] %down% .,V[PF] %up% .,alpha[R] %down% .,alpha[R] %up% .,alpha[A] %down% .,alpha[A] %up% .,f(G,H)))

jpeg(paste ("longterm_main_",format(Sys.time(), "%Y-%m-%d_%H"), "split.jpeg", sep=""), height=7, width=7, res=600, units="in")

ggpubr::ggarrange(p1,p2,p3, nrow=1, ncol=3, common.legend = T, legend="top", labels=c("e", "f", "g"))

dev.off()

rm(p1,p2,p3)

# Figure for dissertation presentation

jpeg(paste ("diss1_",format(Sys.time(), "%Y-%m-%d_%H"), "split.jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(StateVar %in% c("P","S") & Treatment %in% c("H", "W") & Run!=9) %>% ggplot() + geom_hline(yintercept = 0, lty=2) + geom_line(aes(x=time, y=Biomass, color=as.factor(Run))) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Change~"in"~Biomass~(g[N]~m^-2))) + xlab("Days") + scale_colour_manual(guide=F, values=c("black", "grey", "grey", "grey", "grey", "grey", "grey", "red"))

dev.off()

jpeg(paste ("diss2_",format(Sys.time(), "%Y-%m-%d_%H"), "split.jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(StateVar %in% c("P","S") & Treatment %in% c("Interaction") & Run!=9) %>% ggplot() + geom_hline(yintercept = 0, lty=2) + geom_line(aes(x=time, y=Biomass, color=as.factor(Run))) + theme_classic() + facet_grid(StateVar~Treatment, scales="free_y", labeller= as_labeller(facet_names2)) + ylab(expression(Change~"in"~Biomass~(g[N]~m^-2))) + xlab("Days") + scale_colour_manual(guide=F, values=c("black", "grey", "grey", "grey", "grey", "grey", "grey", "red"))

dev.off()

# Another attempt
ggpubr::ggarrange(
  outputALL2_0 %>% filter(time > 8500 & StateVar %in% c("P") & Treatment %in% c("H", "W", "Interaction") & Run %in% c(1,8)) %>% ggplot(aes(x=as.character(Run), y=Biomass, color=time)) + geom_jitter() + theme_classic() + facet_grid(Treatment~StateVar, scales="free"),
  outputALL2_0 %>% filter(time > 8500 &StateVar %in% c("S") & Treatment %in% c("H", "W", "Interaction") & Run %in% c(1,8)) %>% ggplot(aes(x=as.character(Run), y=Biomass, color=time)) + geom_jitter() + theme_classic() + facet_grid(Treatment~StateVar, scales="free"),
  common.legend = T)


outputALL2_0 %>% filter(time > 8500 & StateVar %in% c("P") & Treatment %in% c("H", "W", "Interaction") & Run %in% c(1,8)) %>% ggplot(aes(x=as.character(Run), y=Biomass, color=time)) + geom_jitter() + geom_line() + theme_classic() + facet_grid(Treatment~StateVar, scales="free", labeller= as_labeller(facet_names2)) + scale_x_discrete(labels=c("No plant change", "Plant change")) + xlab("")


facet_names3 <- c(
  `H` = "Add Grashoppers (H)",
  `W` = "Add Earthworms (W)",
  `1` = "No plant change",
  `Interaction` = "Interaction effect",
  `8` = "Plant change"
)

jpeg(paste ("diss1_plant_",format(Sys.time(), "%Y-%m-%d_%H"), ".jpeg", sep=""), height=7, width=7, res=600, units="in")

outputALL2_0 %>% filter(time > 8500 & StateVar %in% c("P") & Treatment %in% c("H", "W", "Interaction") & Run %in% c(1,8)) %>% mutate(time = time-8500) %>% ggplot(aes(x=time, y=Biomass, color=as.character(Run)))+ geom_hline(yintercept = 0, lty=2)+ geom_line(size=2) + theme_classic() + facet_grid(Treatment~Run, scales="free", labeller= as_labeller(facet_names3)) + xlab("Day of the Year")+ ylab(expression(Change~"in"~Plant~Biomass~(g[N]~m^-2))) + scale_color_manual(values=c("black", "green"), guide=F)

dev.off()

outputALL2_0 %>% filter(time > 8500 & StateVar %in% c("S") & Treatment %in% c("H", "W", "Interaction") & Run %in% c(1,8)) %>% mutate(time = time-8500) %>% ggplot(aes(x=time, y=Biomass, color=as.character(Run)))+ geom_hline(yintercept = 0, lty=2)+ geom_line(size=2) + theme_classic() + facet_grid(Treatment~Run, scales="free", labeller= as_labeller(facet_names3)) + xlab("Day of the Year")+ ylab(expression(Change~"in"~Soil~Biomass~(g[N]~m^-2))) + scale_color_manual(values=c("black", "green"), guide=F)