# Code to loop the community invasion over multiple parameter sets -------------

require("FME")
require(tidyverse)
needtoload = T # Load best parameter sets and produce graphics from cluster data
if(needtoload){
  # This analysis runs a code that selects the best model runs
  # and produces graphics in the MS. See the source code.
  source("analysis_of_cluster_data_v2.R")
}

variable_names <- c(
  "H" = "Herbivore (H)" ,
  "M" = "Microbial (M)",
  "N" = "Inorganic (N)",
  "P" = "Plant (P)",
  "P1" = "Plant (P)",
  "W" = "Earthworm (W)",
  "L" = "Litter (L)",
  "R" = "Root (R)",
  "R1" = "Root (R)",
  "S" = "Stable Soil (S)",
  "U" = "Unstable Soil (U)",
  "percentStable" = "Stable SOM (%)",
  "GPP" = "Gross Primary Production",
  "Nmin" = "N Mineralization"
)

## Get and clean the parameter sets and the stable values for each simulation

partoplant = par2 %>% 
              as_tibble() %>%
              filter(Best == "Yes") %>% # select the best runs
              select(PARS, Run, NPAR) %>%
              spread(key = NPAR, value = PARS) %>%
              as.data.frame()

stabletoplant = stable1 %>% 
                  as_tibble() %>%
                  spread(key = NVAR, value = YSTABLE) %>%
                  filter(Run %in% partoplant[,"Run"])%>%
                  as.data.frame()

# Add the model and default starting conditions -----

setplantfunction = function(ID, PARAMS = partoplant, YSTABLE = stabletoplant, SUBSET = T){
  require("FME")
  # Special version of the model that sums up plant effects
  plantmodel_effect <-function(t, y,pars){
    
    with(as.list(c(pars,y)),{
      
      # Model of earthworm growth. From ASA Johnston
      A_W = exp(-Ea/Kappa*(1/LTtemp(t %% 365)-1/Tref_W)) 
      
      # Model of temperature sensitive plant growth from FENG 1990
      A_P = exp(-B/LTtemp(t %% 365))/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/LTtemp(t %% 365))))
      
      # Modeled microbial dynamics MIMICS
      tempC = (LTtemp(t %% 365)-273.15)
      Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
      Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
      Klm = exp(Kslope*tempC + Kint)*Klm_mod
      Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
      
      # Plant root nitrogen uptake
      VR = exp(Vslope*tempC + Vint)*c(VN_mod1, VN_mod2, VN_mod3)
      
      KN = c(KN_mod1, KN_mod2, KN_mod3)
      
      # Cut off at freezing (want night, so add 5C to freezing cut off)
      if(LTtemp(t %% 365) > (273.15 + 5) ){
        AP = A_P*c(A_P_mod1, A_P_mod2, A_P_mod3)
        VH = c(VH_mod1, VH_mod2, VH_mod3)
        alphaR = c(alphaR_mod1, alphaR_mod2, alphaR_mod3)
        alphaA = c(alphaA_mod1, alphaA_mod2, alphaA_mod3)
        th2 = th
      }else{
        AP = rep(0, Nplant) # no AG growth when freezing
        VH = rep(0, Nplant) # no herbivory when freezing
        alphaR = rep(0, Nplant) # no root death when freezing
        alphaA = rep(alphaA_winter, Nplant) # rapid aboveground senescence
        th2 = ifelse(H < Hmin, 0, th_winter) # rapid herbivore death to min pop
      }
      
      Rvec = c(R1, R2, R3)
      Pvec = c(P1, P2, P3)
      
      #Equations
      dL = sum(alphaR*Rvec*Rvec) + sum(alphaA*Pvec) + 
        (1-SUEh)*sum(VH*H*Pvec) + th2*H*H + 
        A_W*tw*W*W - Vlm*L*M/(Klm + M) - 
        A_W*Vlw*L*W - l*L
      
      dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) - 
        tm*M - SUEwm*A_W*Vsw*W*M
      
      dW = AEw*SUEwl*A_W*Vlw*L*W + 
        AEw*SUEws*A_W*Vsw*S*W +
        AEw*SUEws*A_W*Vsw*U*W + 
        AEw*SUEwm*A_W*Vsw*W*M - 
        A_W*tw*W*W
      
      dN = IN - q*N - fi*N + fo*S + 
        (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) - 
        sum(VR*N*Rvec/(KN+N)) + 
        (1-AEh)*SUEh*sum(VH*H*Pvec) +
        (1-AEw)*SUEwl*A_W*Vlw*L*W + 
        (1-AEw)*SUEws*A_W*Vsw*S*W +
        (1-AEw)*SUEws*A_W*Vsw*U*W +
        (1-AEw)*SUEwm*A_W*Vsw*W*M
      
      dS = fstable*tm*M - 
        StU*Vsm*S*M/(Ksm + M) + 
        fstableworm*(1-SUEwl)*A_W*Vlw*L*W +
        fstableworm*(1-SUEws)*A_W*Vsw*S*W -
        A_W*Vsw*S*W +
        fi*N - fo*S +
        fus*U
      
      dU = (1-fstable)*tm*M + 
        StU*Vsm*S*M/(Ksm + M) + 
        (1-fstableworm)*(1-SUEwl)*A_W*Vlw*L*W + 
        (1-fstableworm)*(1-SUEws)*A_W*Vsw*S*W - 
        SUEws*A_W*Vsw*U*W - 
        Vsm*U*M/(Ksm + M) -
        fus*U
      
      dR = VR*N*Rvec/(KN+N) - alphaR*Rvec*Rvec - AP*Rvec
      
      dP = AP*Rvec - VH*H*Pvec - alphaA*Pvec
      
      dH = AEh*SUEh*sum(VH*H*Pvec) - th2*H*H
      
      Nmin = IN + fo*S + 
        (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*U*M/(Ksm + M)) + 
        (1-AEh)*SUEh*sum(VH*H*Pvec) +
        (1-AEw)*SUEwl*A_W*Vlw*L*W + 
        (1-AEw)*SUEws*A_W*Vsw*S*W +
        (1-AEw)*SUEws*A_W*Vsw*U*W +
        (1-AEw)*SUEwm*A_W*Vsw*W*M
      
      GPP = sum(VR*N*Rvec/(KN+N))
      
      SU = S/(S+U)
      
      dy = c(dR, dP, dL, dM, dW, dN, dS, dU, dH)
      
      return(list(dy, Nmin=Nmin, GPP=GPP, percentStable=SU))
      
    }
    )
  }
  
  LTtemp = function(doy){
    
    -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
    
  }
  
  # Select the right run
  RUNcur = PARAMS[ID,1]
  
  # Clean up parameters and add plants
  PARAMScur = subset(PARAMS, Run == RUNcur)[,-1]
  
  paramscur = as.numeric(PARAMScur)
  
  names(paramscur) = colnames(PARAMScur)
  
  plantparams = c("VN_mod", "KN_mod", "A_P_mod", 
                  "VH_mod", "alphaR_mod", "alphaA_mod")
  
  pp0 = pp1 = paramscur[plantparams]
  
  names(pp1) = paste0(names(pp0), "1")
  pp2 = pp1*c(0.5, 0.5,1,0.1,1,1)
  names(pp2) = paste0(names(pp0), "2")
  pp3 = pp1*c(1.5, 1.5,1,10,1,1)
  names(pp3) = paste0(names(pp0), "3")
  
  paramscur = c(paramscur, pp1, pp2, pp3)
  
  # Clean up initial vector and add extra plants
  YSTABLEcur = subset(YSTABLE, Run == RUNcur)[,-1]

  YSTABLEcur = cbind(YSTABLEcur, R2=0, R3=0, P2=0,P3=0)
  
  YSTABLEcur = YSTABLEcur[c("R", "R2", "R3", "P","P2","P3", "L", "M", "W", "N", "S", "U", "H")]
  
  colnames(YSTABLEcur)[c(1,4)] = c("R1", "P1")
  
  ystable = as.numeric(YSTABLEcur)
  
  names(ystable) = colnames(YSTABLEcur)
  
  ystable["H"] = paramscur["Hmin"] # fix error where rounding of H in output makes it zero
  
  # Run model for 5 years just to be sure
  stablerun = ode(y=ystable,times = 1:(365*5),
                  func=plantmodel_effect,
                  parms=paramscur)
    
  yint = stablerun[(365*4),c(2:14)]
  
  years = 100
  
  yint2 = yint
  
  stablerunHW = ode(y=yint2,times = 1:(365*years),
                    func=plantmodel_effect,
                    parms=paramscur)
  
  yint2[1:3] = yint["R1"]/3
  yint2[4:6] = yint["P1"]/3
  
  stablerunHWi = ode(y=yint2,times = 1:(365*years), 
                     func=plantmodel_effect, 
                     parms=paramscur)
  
  yint2 = yint
  yint2["H"] = 0
  
  stablerunW = ode(y=yint2,times = 1:(365*years),
                   func=plantmodel_effect,
                   parms=paramscur)
  
  yint2[1:3] = yint["R1"]/3
  yint2[4:6] = yint["P1"]/3
  
  
  stablerunWi = ode(y=yint2,times = 1:(365*years), 
                    func=plantmodel_effect, 
                    parms=paramscur)
  
  yint2 = yint
  yint2["W"] = 0
  
  stablerunH = ode(y=yint2,times = 1:(365*years),
                   func=plantmodel_effect,
                   parms=paramscur)
  
  yint2[1:3] = yint["R1"]/3
  yint2[4:6] = yint["P1"]/3
  
  
  stablerunHi = ode(y=yint2,times = 1:(365*years), 
                    func=plantmodel_effect, 
                    parms=paramscur)
  
  yint2 = yint
  yint2["W"] = 0
  yint2["H"] = 0
  
  stablerun0 = ode(y=yint2,times = 1:(365*years),
                   func=plantmodel_effect,
                   parms=paramscur)
  
  
  yint2[1:3] = yint["R1"]/3
  yint2[4:6] = yint["P1"]/3
  
  
  stablerun0i = ode(y=yint2,times = 1:(365*years), 
                    func=plantmodel_effect, 
                    parms=paramscur)
  
  rowsel = seq(180, by=365, length=years)
  
  output_base <- data.frame(rbind(stablerun0i,
                                  stablerun0,
                                  stablerunHi,
                                  stablerunH,
                                  stablerunWi,
                                  stablerunW,
                                  stablerunHWi,
                                  stablerunHW))
  if(SUBSET){
    output = output_base[output_base$time %in% rowsel,] 
    output$Treatment = rep(c("N", "H", "W", "HW"), each=years*2)
    output$Species = rep(c("Multiple", "Single"), each=years, times=4)
    output$Run = rep(RUNcur, 4*years*2)
  }else{
    output = output_base
    output$Treatment = rep(c("N", "H", "W", "HW"), each=365*years*2)
    output$Species = rep(c("Multiple", "Single"), each=365*years, times=4)
    output$Run = rep(RUNcur, 4*365*years*2)
  }
  
  return(output)
}

# Simulate and example trajectory of the model to show dynamics -----

simexample = setplantfunction(ID = 10, PARAMS = partoplant, YSTABLE = stabletoplant, SUBSET = F)

png(paste0("simexample_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
simexample %>%
  as_tibble() %>%
  gather(-time, -Treatment, -Species, -Run, key=StateVar, value = Size) %>%
  filter(Species == "Single") %>%
  filter(!(StateVar %in% c("R2", "R3", "P2", "P3"))) %>%
  mutate(time = time/365) %>%
  filter(time < 10) %>%
  ggplot(aes(x=time, y=Size, color= Treatment)) + 
    facet_wrap(.~StateVar, scales="free",
               labeller=labeller(StateVar = variable_names)) +
  geom_line(alpha=0.5) + theme_classic() +
  scale_x_continuous(breaks = c(0,5,10)) +
  theme(legend.position = "top") +
  scale_color_discrete(limits = c("N", "H", "W", "HW"),
                       labels = c("Neither", "Grasshopper",
                                  "Earthworm", "Grasshopper \n and Earthworm"))
dev.off()

rm(simexample)

# Simulate full model ------

require(parallel)
cl <- makeCluster(getOption("cl.cores", 3))

repseq = seq(1, 50, 1)

t1 <- Sys.time()
out1 = parLapply(cl = cl, X=repseq, fun=setplantfunction, PARAMS = partoplant, YSTABLE = stabletoplant, SUBSET = T)
t2 <- Sys.time()

stopCluster(cl)

t2-t1

outdata = do.call("rbind", out1)

saveRDS(outdata, file = "Data/plantdata_31Oct2019.rds")

dir.create(paste0("effect_plots_from_",Sys.Date()))

# Hypothesis #2: Within variation is large relative to across variation

pcadata = outdata %>% as_tibble() %>%
  mutate(R = R1+R2+R3,
         P = P1+P2+P3,
         time = time/365) %>%
  gather(-time, -Treatment,-Run, -Species, key=StateVar, value=Size) %>%
  filter(StateVar %in% 
           c("L", "M", "N", "P", "R", "S", "U")) %>%
  filter(Species == "Single") %>%
  filter(time >99) %>%
  spread(key=StateVar, value=Size) %>%
  select(-time, -Species) %>%
  as.data.frame()

library(vegan)
pca1 = rda(pcadata[,c(3:9)])
summary(pca1)

png(paste0("effect_plots_from_",Sys.Date(),"/","PCA_whyRandS_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
biplot(pca1, type = c("text", "points"),
       xlab = "Component 1 (70%)",
       ylab = "Component 2 (25%)")
dev.off()

png(paste0("effect_plots_from_",Sys.Date(),"/","RunvsTreatment_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
pcadata %>%
  ggplot(aes(x=S, y=R)) + 
  geom_line(aes(group=as.character(Run)), col="grey", alpha=0.8) +
  geom_point(aes(color= as.character(Run), shape=Treatment), size=4) + 
  theme_classic() +
  scale_color_discrete(guide=F) +
  ylab("Root Nitrogen (R)") +
  xlab("Stable Soil Organic Nitrogen (S)") +
  scale_shape_discrete(labels=c("Grashoppers (H)", "Grasshoppers (H) + \n Earthworms (W)",
                                "No animals", "Earthworms (W)")) +
  theme(legend.position = "top")
dev.off()  

cleandata_v1 <- outdata %>% as_tibble() %>%
  mutate(R = R1+R2+R3,
         P = P1+P2+P3,
         time = time/365) %>%
  gather(-time, -Treatment,-Run, -Species, key=StateVar, value=Size) %>%
  spread(key=Treatment, value=Size) %>%
  filter(StateVar %in% 
           c("L", "M", "N", "P", "R", "S", "U", "GPP", "Nmin"))

cleandata_v2 <- cleandata_v1 %>% mutate(EOH1 = (H-N)/N,
                                        EOW1 = (W-N)/N) %>%
  mutate(IE = (HW-N - (H-N + W-N))/N) %>% select(-H:-W) %>%
  gather(-time, -StateVar, -Run, -Species, key=Effect, value=Size)

part1 <- cleandata_v2 %>% 
  group_by(StateVar, Species, Effect, Run) %>%
  summarize(MS = sum(Size)) # OK because a change high-low should cancel

# Find quantiles (based on position) and then select them into data frame
plotdata_v1 <- part1 %>% arrange(desc(MS)) %>% 
  do(slice(., 12)) %>% rename(upper = Run) %>% select(-MS) %>%
  full_join(
    part1 %>% arrange(desc(MS)) %>% 
      do(slice(., 38)) %>% rename(lower = Run) %>% select(-MS)
  ) %>%
  full_join(
    part1 %>% arrange(desc(MS)) %>% 
      do(slice(., 24)) %>% rename(median = Run) %>% select(-MS)
  ) %>%
  gather(-StateVar:-Effect, key=QT, value=Run) %>%
  left_join(cleandata_v2) %>% # add back in data, 54 groups, 100 times, 3 levels
  select(-Run)

plotdata_v2 <- plotdata_v1 %>%
  ungroup() %>%
  mutate(EFQT = paste0(Effect, QT)) %>%
  select(-Effect, -QT) %>%
  spread(key=EFQT, value=Size)

# part 2 is the data with the actual Runs to select

#Plot the result over years
png(paste0("effect_plots_from_",Sys.Date(),"/","Effectsovertime_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
plotdata_v1 %>% 
  spread(key=QT, value=Size) %>%
  ggplot(aes(x=time, y=median,color=Effect)) +
  geom_line(aes(linetype = Species)) +
  geom_linerange(aes(ymin=lower, ymax=upper), alpha=0.1)+
  facet_wrap(.~StateVar, scales="free") + theme_classic() +
  ylab("Effect Size (%)") + xlab("Year") + 
  scale_color_discrete(labels=c("Effect of grasshoppers", "Effect of earthworms",
                                "Interaction Effect")) +
  theme(legend.position = "top")
dev.off()

plotdata4 = plotdata_v2 %>% filter(time > 3.5 & time < 5) %>%
  filter(StateVar %in% c("P", "R", "S", "U"))
plotdata100 = plotdata_v2 %>% filter(time > 99) %>%
  filter(StateVar %in% c("P", "R", "S", "U"))

png(paste0("effect_plots_from_",Sys.Date(),"/","clean_IEovertime_earthworm_",Sys.Date(),".png"), width=8, height=6, units="in", res=300)
plotdata_v2 %>%
  filter(StateVar %in% c("P", "R", "S", "U")) %>%
  ggplot(aes(x=EOW1median, y=IEmedian, color=Species)) +
  geom_vline(xintercept=0,linetype = "dashed", color="grey") +
  geom_hline(yintercept=0,linetype = "dashed", color="grey") +
  geom_path(aes(linetype=Species)) +
  geom_point(data = plotdata4 ,shape=1, size=4) +
  geom_point(data = plotdata100,shape=19, size=4) +
  geom_errorbar(data = plotdata100, aes(ymin=IElower, ymax=IEupper), width=0) +
  geom_errorbarh(data = plotdata100, aes(xmin=EOW1lower, xmax=EOW1upper), height=0) +
  facet_wrap(.~StateVar, scales="free",
             labeller=labeller(StateVar = variable_names)) + 
  theme_classic() + theme(legend.position= "top") +
  xlab("Earthworm Effect (%)") + ylab("Interaction Effect (%)") +
  scale_color_manual(values=c("red","black"))
dev.off()

plotdata4 = plotdata_v2 %>% filter(time > 3.5 & time < 5)
plotdata100 = plotdata_v2 %>% filter(time > 99)

png(paste0("effect_plots_from_",Sys.Date(),"/","IEovertime_earthworm_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
plotdata_v2 %>%
  ggplot(aes(x=EOW1median, y=IEmedian, color=Species)) +
  geom_vline(xintercept=0,linetype = "dashed", color="grey") +
  geom_hline(yintercept=0,linetype = "dashed", color="grey") +
  geom_path(aes(linetype=Species)) +
  geom_point(data = plotdata4 ,shape=1, size=4) +
  geom_point(data = plotdata100,shape=19, size=4) +
  geom_errorbar(data = plotdata100, aes(ymin=IElower, ymax=IEupper), width=0) +
  geom_errorbarh(data = plotdata100, aes(xmin=EOW1lower, xmax=EOW1upper), height=0) +
  facet_wrap(.~StateVar, scales="free",
             labeller=labeller(StateVar = variable_names)) + 
  theme_classic() + theme(legend.position= "top") +
  xlab("Earthworm Effect (%)") + ylab("Interaction Effect (%)") +
  scale_color_manual(values=c("red","black"))
dev.off()

png(paste0("effect_plots_from_",Sys.Date(),"/","IEovertime_grasshopper_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
plotdata_v2 %>%
  ggplot(aes(x=EOH1median, y=IEmedian, color=Species)) +
  geom_vline(xintercept=0,linetype = "dashed", color="grey") +
  geom_hline(yintercept=0,linetype = "dashed", color="grey") +
  geom_path(aes(linetype=Species)) +
  geom_point(data = plotdata4 ,shape=1, size=4) +
  geom_point(data = plotdata100,shape=19, size=4) +
  geom_errorbar(data = plotdata100, aes(ymin=IElower, ymax=IEupper), color="black",width = 0) +
  geom_errorbarh(data = plotdata100, aes(xmin=EOH1lower, xmax=EOH1upper), color="black", height=0) +
  facet_wrap(.~StateVar, scales="free",
             labeller=labeller(StateVar = variable_names)) + 
  theme_classic() + theme(legend.position= "top") +
  xlab("Grasshopper Effect (%)") + ylab("Interaction Effect (%)") +
  scale_color_manual(values=c("red","black"))
dev.off()

png(paste0("effect_plots_from_",Sys.Date(),"/","PlantCommunity_n100_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
outdata %>% as_tibble() %>% filter(time > (100-1)*365) %>%
  select(R1:P3, Treatment, Species, Run) %>%
  gather(-Treatment,-Species,-Run, key= Pool, value=Size) %>%
  separate(Pool, into=c("Pool", "ID"), sep=1) %>%
  filter(Pool == "P", Species=="Multiple") %>%
  ggplot(aes(x=Treatment, y=Size, fill=ID)) +
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(labels= scales::percent_format()) +
  theme_classic() +
  scale_x_discrete(limits=c("N","H","W", "HW")) +
  facet_wrap(.~Run) +
  scale_fill_discrete(labels=c("Base species",
                                "Slow-growing species",
                                "Fast-growing species")) +
  theme(legend.position = "top")
dev.off()

png(paste0("effect_plots_from_",Sys.Date(),"/","PlantCommunity_n4_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
outdata %>% as_tibble() %>% filter(time > 3.5*365 & time < 5*365) %>%
  select(R1:P3, Treatment, Species, Run) %>%
  gather(-Treatment,-Species,-Run, key= Pool, value=Size) %>%
  separate(Pool, into=c("Pool", "ID"), sep=1) %>%
  filter(Pool == "P", Species=="Multiple") %>%
  ggplot(aes(x=Treatment, y=Size, fill=ID)) +
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(labels= scales::percent_format()) +
  theme_classic() +
  scale_x_discrete(limits=c("N","H","W", "HW")) +
  facet_wrap(.~Run)
dev.off()


CWM_per <- outdata %>% as_tibble() %>% filter(time > (100-1)*365) %>%
  select(R1:P3, Treatment, Species, Run) %>%
  gather(-Treatment,-Species,-Run, key= Pool, value=Size) %>%
  separate(Pool, into=c("Pool", "ID"), sep=1) %>%
  filter(Pool == "P", Species=="Multiple") %>%
  group_by(Treatment, Run) %>%
  summarise (sm = sum(Size))

variable_names_trends = c("1" = "Base species",
                          "2" = "Slow-growing species",
                          "3" = "Fast-growing species",
                          "N" = "Neither",
                          "H" = "Grasshopper",
                          "W" = "Earthworm",
                          "HW"= "Grasshopper \n and Earthworm")

png(paste0("effect_plots_from_",Sys.Date(),"/","PlantCommunity_trends_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
outdata %>% as_tibble() %>%
  select(time, R1:P3, Treatment, Species, Run) %>%
  gather(-time, -Treatment,-Species,-Run, key= Pool, value=Size) %>%
  separate(Pool, into=c("Pool", "ID"), sep=1) %>%
  filter(Pool == "P", Species=="Multiple") %>%
  filter(time %in% seq(180, by=365, length=100)[seq(1, 100, 10)]) %>%
  mutate(time = time/365) %>%
  left_join(
    CWM_per
  ) %>%
  mutate(Prop = Size/sm) %>%
  ggplot(aes(x=time, y=Prop, group= Run)) +
  geom_line(alpha = 0.5, color = "grey") + 
  theme_classic() +
  facet_grid(ID~Treatment,
             labeller = as_labeller(variable_names_trends)) +
  theme(legend.position = "top")
dev.off()


CWM =  outdata %>% as_tibble() %>% filter(time > (100-1)*365) %>%
  select(R1:P3, Treatment, Species, Run) %>%
  gather(-Treatment,-Species,-Run, key= Pool, value=Size) %>%
  separate(Pool, into=c("Pool", "ID"), sep=1) %>%
  filter(Pool == "P", Species=="Multiple") %>%
  mutate(Pool = paste0(Pool,ID)) %>%
  left_join(
    CWM_per
  ) %>%
  mutate(freq = Size / sm) %>%
  select(Treatment, Run, Pool, freq) %>%
  spread(key=Pool, value=freq) %>%
  right_join(
    partoplant %>% as_tibble() %>% 
      gather(-Run, key = NPAR, value=PAR ) %>%
      filter((NPAR %in% c("VN_mod", "KN_mod", "A_P_mod","VH_mod", "alphaR_mod", "alphaA_mod")))
  ) %>% 
  left_join(
    data.frame(NPAR= c("VN_mod", "KN_mod", "A_P_mod", 
                      "VH_mod", "alphaR_mod", "alphaA_mod"),
               mod2 = c(0.5, 0.5,1,0.1,1,1),
               mod3 = c(1.5, 1.5,1,10,1,1))
  ) %>%
  mutate(CW = P1*PAR +P2*PAR*mod2 + P3*PAR*mod3) %>%
  select(Treatment, Run, NPAR, CW)

png(paste0("effect_plots_from_", Sys.Date(),"/","PlantCommunityWeightedMean_",Sys.Date(),".png"), width=8, height=8, units="in", res=300)
CWM %>%
  ggplot(aes(x=Treatment, y=CW)) + 
  geom_jitter(aes(color=Run), width = 0.2, height=0) +
  facet_wrap(.~NPAR, scales="free") + theme_classic() +
  stat_summary(fun.y = mean, color="black", size=4, shape=17, geom="point") +
  theme(legend.position = "none") +
  ylab("Community Weighted Mean Parameter")  +
  scale_x_discrete(labels = c("Neither", "Grasshopper",
                              "Earthworm", "Grasshopper \n and Earthworm"),
                   limits = c("N","H","W", "HW"))
dev.off()

variable_names3 <- c("KN_mod" = "Plant~Half~Saturation~(K[NP])",
                     "VN_mod" = "Plant~Max~Uptake~(V[NP]^0)",
                     "VH_mod" = "Herbivory~Rate~(V[HP])")

png(paste0("effect_plots_from_",Sys.Date(),"/","clean_PlantCommunityWeightedMean_",Sys.Date(),".png"), 
    width=10, height=4, units="in", res=300)
CWM %>% filter(NPAR %in% c("VH_mod", "VN_mod", "KN_mod")) %>%
  ggplot(aes(x=Treatment, y=CW)) + 
  geom_line(aes(group=Run), alpha=0.2) +
  geom_jitter(aes(color = Run), width = 0.2, height=0) +
  facet_wrap(.~NPAR, scales="free", 
             labeller = as_labeller(variable_names3,label_parsed)) + 
  theme_classic() +
  stat_summary(fun.y = mean, color="black", size=4, shape=17, geom="point") +
  theme(legend.position = "none") +
  ylab("Community Weighted Mean Parameter") +
  scale_x_discrete(labels = c("Neither", "\n Grasshopper",
                              "Earthworm", "\n Grasshopper \n and Earthworm"),
                   limits = c("N","H","W", "HW"))
dev.off()


# Plot the baseline params versus effect size at 100 years ----------------

variable_names4 <- c("KN_mod" = "K[NP]",
                     "VN_mod" = "V[NP]^0",
                     "VH_mod" = "V[HP]", 
                     "A_P_mod" = "A[P]^0",
                     "alphaA_mod" = "alpha[A]",
                     "alphaR_mod" = "alpha[R]",
                     "H" = "Herbivore~(H)" ,
                     "M" = "Microbial~(M)",
                     "N" = "Inorganic~(N)",
                     "P" = "Plant~(P)",
                     "W" = "Earthworm~(W)",
                     "L" = "Litter~(L)",
                     "R" = "Root~(R)",
                     "S" = "Stable~Soil~(S)",
                     "U" = "Unstable~Soil~(U)",
                     "GPP" = "Gross~Primary~Production",
                     "Nmin" = "Mineralization")

png(paste0("effect_plots_from_",Sys.Date(),"/","effect_of_plant_PAR_",Sys.Date(),".png"), 
    width=8, height=11, units="in", res=300)
part1 %>% 
  filter(Species=="Single" & Effect == "IE") %>%
  left_join(
    partoplant %>% as_tibble() %>% 
      gather(-Run, key = NPAR, value=PAR ) %>%
      filter((NPAR %in% c("VN_mod", "KN_mod", "A_P_mod",
                          "VH_mod", "alphaR_mod", "alphaA_mod"))) %>%
      spread(key = NPAR, value = PAR)
  ) %>%
  ungroup() %>%
  select(-Species, -Effect) %>%
  gather(-StateVar, -Run, - MS, key=NPAR, value=PAR) %>%
  filter(MS > -1e3 & MS < 1e3) %>%
  ggplot(aes(x=PAR, y=MS)) + 
  geom_point() + 
  facet_grid(StateVar~NPAR, scale="free",
             labeller = as_labeller(variable_names4, 
                                    default = label_parsed)) + 
  theme_bw() +
  stat_smooth(method="lm", col="red") +
  ylab("Interactive Effect (%)") +
  xlab("Parameter Value")
dev.off()

# Compare single vs. multiple effects within Runs over time: 

compdata <- cleandata_v2 %>%
  spread(key=Species, value=Size) %>%
  mutate(Diff = abs(Multiple)-abs(Single)) %>%
  filter(Effect =="IE")

compdata_a <- compdata %>% 
  group_by(StateVar, Run) %>%
  summarize(MS = sum(Diff)) # OK because a change high-low should cancel

# Find quantiles (based on position) and then select them into data frame
compdata_b <- compdata_a %>% arrange(desc(MS)) %>% 
  do(slice(., 12)) %>% rename(upper = Run) %>% select(-MS) %>%
  full_join(
    compdata_a %>% arrange(desc(MS)) %>% 
      do(slice(., 38)) %>% rename(lower = Run) %>% select(-MS)
  ) %>%
  full_join(
    compdata_a %>% arrange(desc(MS)) %>% 
      do(slice(., 24)) %>% rename(median = Run) %>% select(-MS)
  ) %>%
  gather(-StateVar, key=QT, value=Run) %>%
  left_join(compdata) %>% # add back in data, 54 groups, 100 times, 3 levels
  select(-Run, -Multiple, -Single) %>%
  mutate(EFQT = paste0(Effect,QT)) %>%
  select(-Effect, -QT) %>%
  spread(key=EFQT, value=Diff)

png(paste0("effect_plots_from_",Sys.Date(),"/","MultSingSP_",Sys.Date(),".png"), 
    width=8, height=8, units="in", res=300)
compdata_b %>% mutate(time2 = round(time)) %>%
  filter(time2 %in% seq(0, 100, by=5))%>% 
  ggplot(aes(x=time, y=IEmedian)) + 
  geom_pointrange(aes(ymin=IElower, ymax=IEupper)) + 
  facet_wrap(.~StateVar, scales="free",
             labeller = as_labeller(variable_names4,
                                    default = label_parsed)) +
  theme_bw() +
  ylab("Difference in Interactive Effect (Multiple - Single)") + 
  xlab("Year")
dev.off()
