# Model with baseline parameters
# Feb 28 2020

# The baseline version of this model doesn't do as good of a job replicating the data
# The simple model is better...

source("treatment_apply.R")
require(deSolve)
singlerun <- function(idx){
  ID = round(runif(1)*1e8,0) # set round ID
  ID2 = round(runif(1)*1e8,0) # set backup round ID
  paramscur = params # import params
  spinyears = c(500, 1500, 3000) # spin length
  ERRRRR = T # set loop
  
  while(ERRRRR){
    
    ystable = yint
    
    for(i in 1:3){
      # Try first quick spin
      stablerun = ode(y=ystable,times = 1:(365*spinyears[i]), 
                      func=plantmodel, 
                      parms=paramscur)
      
      stableannual <- as.data.frame(stablerun)
      stableannual<- stableannual[stableannual$time %in% seq(180, 365*spinyears[i], by=365),]
      
      stabletest <- stableannual[dim(stableannual)[1],c(-1, -11, -12, -13)] -
        stableannual[(dim(stableannual)[1]-1),c(-1, -11, -12, -13)]
      
      ystable = stablerun[365*(spinyears[i]-1),c(-1, -11, -12, -13)]
      
      if(!all(ystable > 1e-4)) break
      if(max(abs(stabletest)) < 1e-3) break
    }
    
    ERRRRR = !all(c(max(abs(stabletest)) < 1e-3, ystable > 1e-4))
  }
  
  ystable2 = ystable
  
  output2_WE_Return = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur)
  
  output2_WE_Remove = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                          events = list(data=eshock_WE))
  
  output2_HW = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                   events = list(data=eadd))
  
  output2_H = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eshock))
  
  ystable["H"] = 0
  
  output2_W = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eadd))
  
  output2_0 = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eshock))
  
  output2_WE_Hopper = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur)
  
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
                YSTABLE = c(unname(ystable), rep(-2, lg - length(ystable2))))
  
  return(out3)
  
}

baseline <-singlerun(1)

write.csv(baseline, "Data/complexmodel_baseline.csv", row.names = F)

View(baseline)

library(tidyverse)

data3 <- baseline %>% as_tibble() %>% select(-PARS, -YSTABLE, -Run2) %>% 
  filter(time!=-1) %>% 
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>%
  select(time, Treatment, Expt,Run, W,P,H,M,N) %>%
  gather(- Expt, -Treatment, -time, -Run, key=StateVar, value=Model)


data4 = data3 %>% left_join(
  data3 %>% group_by(Expt,StateVar) %>% 
    summarize(Msd = mad(Model))) 


datatomodel3 = datatomodel2 %>% 
  group_by(time, Treatment, StateVar) %>% 
  summarize(mBio = median(Biomass),
            lowBio = quantile(Biomass, 0.25),
            upBio = quantile(Biomass, 0.75)
  ) %>% ungroup() %>%
  filter(!(Treatment %in% c("RmH", "RtH") & StateVar=="M")
  ) %>%
  # Get rid of experimental time before model "treatments" are imposed (i.e before 290 days)
  filter(time != 114)

errord = datatomodel3 %>% left_join(data4) %>%
  filter(!is.na(Model)) %>% # remove any cases without simulation
  filter(Msd !=0) %>% # remove hopper treatments where standard deviation is zero by definition
  group_by(Run) %>% # create group by Run so standardize function works
  summarize(fit = sqrt(sum(((Model-mBio)/Msd)^2))) %>% # function for error ==> minimum is better
  mutate(Best = ifelse(fit < quantile(fit, cut), "Yes", "No"))


(databest = data4 %>% 
  select(-Msd) %>%
  group_by(Treatment, time, StateVar) %>%
  summarize(Model = median(Model)) %>%
  ungroup() %>%
  left_join(datatomodel3) %>% 
  filter(!is.na(mBio)) %>% #add data and get ride of model draws without comparison
  filter(Treatment %in% c("N", "W", "H", "HW")) %>%
  group_by(StateVar) %>%
  summarize(rsq= cor(Model, mBio)^2)
)

databest2 = databest %>% 
  mutate(RSQ = paste0("italic(R^2)==",round(rsq,2)))

parplot1 = data4 %>% 
  select(-Msd) %>%
  group_by(Treatment, time, StateVar) %>%
  summarize(
    lowModel = quantile(Model, 0.25),
    upModel = quantile(Model, 0.75),
    Model = median(Model)
  ) %>%
  ungroup() %>%
  left_join(datatomodel3) %>% filter(!is.na(mBio)) %>%
  mutate(Year = time/365)

# Create a dummy data frame so the axes are the same size
parplot2 = parplot1 %>% select(StateVar, upModel, upBio) %>%
  gather(-StateVar, key=Type, value=V1) %>%
  group_by(StateVar) %>%
  summarize(Model_1 = max(V1)) %>%
  mutate(mBio_1 = Model_1,
         Model_2 = 0,
         mBio_2 = 0) %>%
  gather(-StateVar, key=Type, value=V1) %>%
  separate(Type, into=c("Type", "Level"), sep="_") %>%
  spread(key=Type, value=V1)

variable_names <- c(
  "H" = "Herbivore" ,
  "M" = "Microbial",
  "N" = "Inorganic",
  "P" = "Plant",
  "W" = "Earthworm",
  "L" = "Litter",
  "R" = "Root",
  "S" = "Stable Soil",
  "U" = "Unstable Soil",
  "GPP" = "GPP",
  "Nmin" = "Mineralization"
)

parplot1 %>%
  ggplot(aes(x=Model, y=mBio)) +
  geom_blank(data=parplot2) + 
  geom_abline(intercept=0, slope=1) + 
  geom_errorbar(aes(ymin=lowBio, ymax=upBio, color=Treatment))+
  geom_errorbarh(aes(xmin=lowModel, xmax=upModel, color=Treatment))+
  geom_point(alpha = 0.6, aes(size=Year, color=Treatment)) +
  facet_wrap(.~StateVar, scale="free", 
             labeller=labeller(StateVar = variable_names)) + 
  theme_classic() +
  ylab("Field Biomass") + xlab("Model Biomass") + 
  geom_text(
    data    = databest2,
    mapping = aes(x = -Inf, y = Inf, label = RSQ),
    hjust   = -1,
    vjust   = 1, parse=T
  ) +
  scale_color_discrete(breaks=c("N", "H", "W", "HW",
                                "Rt", "RmW", "RtH", "RmH"),
                       labels=c("None (Expt)", "Herbivore (Expt)",
                                "Earthworm (Expt)",
                                "Both (Expt)", 
                                "Worm control",
                                "Worm removal", 
                                "Herbivore control",
                                "Herbivore removal"))
