# Analysis of cluster data
# Mar 28 2020

require(tidyverse)
verbose = F
loadcsv = T

# NEED TO CLEAN THIS UP. STILL LOTS OF OLDER CODE

# Load parameter vector and data -----

datatomodel2 = read_csv("Data/datatomodel2.csv")

# Load (Cluster) data -----

if(loadcsv){
  dirtoload = "Model_parameter_reps/"
  
  ftoload = list.files(dirtoload)
  
  listOfDataFrames <- vector(mode = "list", length = length(ftoload))
  
  for(ii in 1:length(ftoload)){
    listOfDataFrames[[ii]] <- read.csv(paste0(dirtoload,ftoload[ii]))
  }
  
  data2 <- do.call("rbind", listOfDataFrames)
  
  rm(listOfDataFrames)
  
  
  data2$Run = paste0(data2$Run, data2$Run2)
  
  # check for duplicate Run names and replace them
  table(data2$Run)[table(data2$Run) > 160]

  saveRDS(data2, file = paste0("Data/fullmodeloutput_",Sys.Date(),".rds"))
}


# The following .rds file was created by using the function "rbind" to join individual model outputs
# data2 <- readRDS("Data/fullmodeloutput_12Aug2019.rds")

nparams = names(params)

# get parameter values
par1 <- data2 %>% select(PARS, Run) %>% 
  mutate(NPAR = rep(c(nparams, rep("REMOVE", 127)), length(unique(data2$Run)))) %>%
  filter(NPAR != "REMOVE")

# get stable state variable values

stable1 <- data2 %>% select(YSTABLE, Run) %>% filter(YSTABLE!=-2) %>% 
  mutate(NVAR = rep(colnames(data2)[2:8], length(unique(data2$Run))))

# get model data to compare with empirical data
data3 <- data2 %>% select(-PARS, -YSTABLE, -Run2) %>% filter(time!=-1) %>% rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>%
  # mutate(P = P1 + P2) %>%
  select(time, Treatment, Expt,Run, W,P,H,M,N) %>%
  gather(- Expt, -Treatment, -time, -Run, key=StateVar, value=Model)

# Match and select best runs ----

# calculate median absolute deviation (instead of standard deviation) for unique model runs
# based on van der Vaart et al. 2015

useexptdata = T
useextradata = T

if(useexptdata){
  if(useextradata){
    data4 = data3 %>% left_join(
      data3 %>% group_by(Expt,StateVar) %>% 
        summarize(Msd = mad(Model))) 
  }else{
    data4 = data3 %>% left_join(
      data3 %>% group_by(Expt,StateVar) %>% 
        summarize(Msd = mad(Model)))%>%
      filter(Treatment %in% c("N", "H", "W",  "HW"))
    }
  
}else{
  data4 = data3 %>% left_join(
    data3 %>% group_by(Expt,StateVar) %>% 
      summarize(Msd = mad(Model)))%>%
    filter(Treatment %in% c("RmH", "RmW", "Rt",  "RtH"))
}

rm(data3)

# calculate mean of each treatment by StateVar by time in the model and mean of each variable

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

  
(Nobs = length(unique(data4$Run)))
(cut = (50/Nobs))

errord = datatomodel3 %>% left_join(data4) %>%
  filter(!is.na(Model)) %>% # remove any cases without simulation
  filter(Msd !=0) %>% # remove hopper treatments where standard deviation is zero by definition
  group_by(Run) %>% # create group by Run so standardize function works
  summarize(fit = sqrt(sum(((Model-mBio)/Msd)^2))) %>% # function for error ==> minimum is better
  mutate(Best = ifelse(fit < quantile(fit, cut), "Yes", "No"))

table(errord %>% select(Best))

# Save ID of selected Runs 
# errord %>% filter(Best == "Yes") %>% write_csv(paste0("Data/selected_runs_", Sys.Date(),".csv"))


if(verbose){
  # check for extinctions (None observed, unless expected (i.e. no herbivores in N or W treatments))
  data4 %>% as_tibble() %>%
    left_join(errord %>% select(Run, Best)) %>% filter(Best=="Yes") %>% # pick best runs
    select(-Msd, -Best) %>%
    filter(StateVar %in% c("H", "W")) %>%
    group_by(Treatment, Run, StateVar) %>%
    summarize(Model = prod(Model)) %>% # use product to search for zeros in each vector
    spread(key = Treatment, value=Model) %>% View()
}


# calculate R2 by treatment using the average model runs of best models

(databest = data4 %>% 
  left_join(errord %>% select(Run, Best)) %>% filter(Best=="Yes") %>% # pick best runs
  select(-Msd, -Best) %>%
  group_by(Treatment, time, StateVar) %>%
  summarize(Model = median(Model)) %>%
  ungroup() %>%
  left_join(
    datatomodel3) %>% 
    filter(!is.na(mBio)
           ) %>% #add data and get ride of model draws without comparison
  group_by(StateVar) %>%
  summarize(rsq= cor(Model, mBio)^2))

databest2 = databest %>% 
  mutate(RSQ = paste0("italic(R^2)==",round(rsq,2)))

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

variable_names_Expt <- c(
  "1" = "Other field data",
  "2" = "Cage experiment"
)

# Plot results of selection -----------------------------------------------

#create directory for results
dir.create(paste0("plots_from_",Sys.Date()))

parplot1 = data4 %>% 
  left_join(
    errord %>% select(Run, Best)
  ) %>% 
  filter(Best=="Yes") %>% # pick best runs
  select(-Msd, -Best) %>%
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

png(paste0("plots_from_",Sys.Date(),"/","model_expt_", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
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
dev.off()

if(verbose){
  data4 %>% 
    left_join(
      errord %>% select(Run, Best)
    ) %>% 
    filter(Best=="Yes") %>% # pick best runs
    select(-Msd, -Best) %>%
    group_by(Treatment, time, StateVar, Expt) %>%
    summarize(
      lowModel = quantile(Model, 0.25),
      upModel = quantile(Model, 0.75),
      Model = median(Model)
    ) %>%
    ungroup() %>%
    left_join(datatomodel3) %>% filter(!is.na(mBio)) %>%
    select(time, Treatment, StateVar, mBio, Model)
}

# Compare effect sizes and directions!
png(paste0("plots_from_",Sys.Date(),"/","EffectSize_compare", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)

datatomodel3 %>% select(time, Treatment, StateVar, mBio) %>%
  spread(key = Treatment, value=mBio) %>%
  mutate(EOH1 = (H-N),
         EOW1 = (W-N),
         EOW2 = (Rt - RmW),
         EOH2 = (RtH - RmH)) %>%
  select(time, StateVar, EOH1:EOH2) %>%
  gather(-time, -StateVar, key=Effect, value=ExptSize) %>%
  filter(!is.na(ExptSize)) %>%
  left_join(
    #*********************************************#
    data4 %>% as_tibble() %>%
      left_join(
        errord %>% select(Run, Best)
      ) %>% 
      filter(Best=="Yes") %>% # pick best runs
      select(-Msd, -Best, -Expt) %>%
      spread(key = Treatment, value=Model) %>%
      group_by(time, StateVar) %>%
      summarize(EOH1 = median(H-N),
                EOW1 = median(W-N),
                EOW2 = median(Rt - RmW),
                EOH2 = median(RtH - RmH)) %>%
      gather(-time, -StateVar, key=Effect, value=ModelSize)
    #*********************************************#
  ) %>%
  ggplot(aes(x=ModelSize, y=ExptSize, color=Effect)) + 
  geom_point(aes(size=time)) + 
  facet_wrap(.~StateVar, scale="free") + theme_classic() +
  geom_vline(xintercept=0, color="grey", linetype="dashed")+
  geom_hline(yintercept=0, color="grey", linetype="dashed")

dev.off()

prepselect <- data4 %>% as_tibble() %>%
  left_join(
    errord %>% select(Run, Best)
  ) %>% 
  filter(Best=="Yes") %>%
  select(-Msd, -Expt, -Best) %>%
  group_by(Treatment, StateVar, Run) %>%
  summarize(MS = sum(Model)) %>%
  arrange(desc(MS))

# 50 runs, so take Run 25, 50, 75
ppdata <- prepselect %>% 
  do(slice(., 12)) %>% rename(upper = Run) %>% select(-MS) %>%
  full_join(
    prepselect %>% 
      do(slice(., 25)) %>% rename(median = Run) %>% select(-MS)
  ) %>%
  full_join(
    prepselect %>% 
      do(slice(.,38)) %>% rename(lower = Run) %>% select(-MS)
  ) %>%
  gather(-Treatment, -StateVar, key=QT, value=Run) %>%
  left_join(
    data4 %>% as_tibble() %>%
      left_join(
        errord %>% select(Run, Best)
      ) %>% 
      filter(Best=="Yes") %>%
      select(-Msd, -Best)
  ) %>% select(-Run)

ppdata2 <- datatomodel3 %>% 
  filter(Treatment %in% c("N", "H", "W", "HW")) %>%
  mutate(Year = time/365)

png(paste0("plots_from_",Sys.Date(),"/","model_expt_time", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
ppdata %>% filter(Expt == 2) %>%
  mutate(Year = time/365) %>%
  ggplot(aes(x=Year, color=Treatment)) +
  geom_line(aes(y=Model,linetype=QT)) + theme_classic() + 
  facet_wrap(.~StateVar, scales="free") +
  scale_linetype_manual(values=c("dashed", "solid", "dashed")) +
  geom_point(data = ppdata2, aes(y=mBio)) + 
  geom_errorbar(data = ppdata2, aes(ymin=lowBio, ymax = upBio))
dev.off()


# Look at state variable distribution for selected runs

if(verbose){
  png(paste0("plots_from_",Sys.Date(),"/","state_var_dist_", Sys.Date(),".png"), width=8,
      height=5, units="in", res=300)
  data4 %>% 
    left_join(errord %>% select(Run, Best)) %>% 
    filter(Best =="Yes") %>%
    ggplot(aes(x=Model, fill=as.factor(time))) + 
    geom_histogram() + facet_wrap(.~StateVar, scales="free")
  dev.off()
  
  # Look at selected Runs
  errord %>% ggplot(aes(x=fit, fill = Best)) + geom_histogram(binwidth = 10) + theme_classic()
}

# Plot root versus shoot results for the experimental runs

png(paste0("plots_from_",Sys.Date(),"/","shootTOroot", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
data2 %>% 
  as_tibble() %>%
  select(time, R, P, Treatment, Run, Expt) %>%
  left_join(errord %>% select(Run, Best)) %>% 
  filter(Best =="Yes") %>%
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN, -Best) %>%
  left_join(
    data.frame(time = unique(data2$time),
               Season = c("Spring", "Summer", "Fall", "Fall", "Spring",
                          "Spring", "Summer", "Fall", "Fall", "Fall",
                          "Spring", "Spring", "Spring", "Summer", "Summer",
                          "Fall", "Fall", "Spring", "Summer", "Fall"))
  ) %>%
  mutate(Expt = ifelse(Expt == 2, "Cage experiment", "Field plots")) %>%
  filter(time > 365) %>%
  ggplot(aes(x = R, y=P, color=Treatment)) + 
  geom_point(shape=1) + theme_classic() +
  facet_wrap(Expt~Season, scales="free") + 
  ylab("Shoot nitrogen (P)") + xlab("Root nitrogen (R)") +
  scale_color_discrete(breaks=c("N", "H", "W", "HW",
                                "Rt", "RmW", "RtH", "RmH"),
                       labels=c("None (Expt)", "Herbivore (Expt)",
                                "Earthworm (Expt)",
                                "Both (Expt)", 
                                "Worm control",
                                "Worm removal", 
                                "Herbivore control",
                                "Herbivore removal"))
dev.off()


png(paste0("plots_from_",Sys.Date(),"/","unstableTOstbale", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
data2 %>% 
  as_tibble() %>%
  select(time, U, S, Treatment, Run, Expt) %>%
  left_join(errord %>% select(Run, Best)) %>% 
  filter(Best =="Yes") %>%
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN, -Best) %>%
  filter(time > 365) %>%
  ggplot(aes(x = U, y=S, color=Treatment)) + 
  geom_point(shape=1) + theme_classic() +
  ylab("Unstable soil nitrogen (U)") + xlab("Stable soil nitrogen (S)") +
  scale_color_discrete(breaks=c("N", "H", "W", "HW",
                                "Rt", "RmW", "RtH", "RmH"),
                       labels=c("None (Expt)", "Herbivore (Expt)",
                                "Earthworm (Expt)",
                                "Both (Expt)", 
                                "Worm control",
                                "Worm removal", 
                                "Herbivore control",
                                "Herbivore removal"))
dev.off()

# look at posterior

par2 = par1 %>% left_join(errord) %>% 
  mutate(Manip = ifelse(NPAR %in% names(params[c(1:6,8:31)]),
                        "Yes", "No"))

parVT = par2 %>% filter(Manip=="Yes") %>% filter(!is.na(fit))
pvec = unique(parVT$NPAR)
sigvec = rep(-1, length(pvec))
for(i in 1:length(pvec)){
  vec1 = subset(parVT,NPAR == pvec[i] & Best=="No")$PARS
  vec2 = subset(parVT,NPAR == pvec[i] & Best=="Yes")$PARS
  sigvec[i] = var.test(vec1,vec2, ratio=1, alternative = "greater")$p.value
  
}

rm(vec1, vec2)

sigframe = data.frame(NPAR = pvec, 
                      pval = ifelse(sigvec < 0.05, "*", ""))

par_plot = par2 %>% filter(Manip=="Yes") %>% 
  left_join(data.frame(NPAR = names(params),litvalue = unname(params))) %>% # add literature values
  mutate(StdValue = PARS/litvalue) %>% # standardize estimates for easy display
  select(NPAR, Best, StdValue) %>% group_by(NPAR, Best) %>%
  summarize(lower = quantile(StdValue, 0.055), 
            upper = quantile(StdValue, 0.945),
            median = quantile(StdValue, 0.5))

parameters = c("A[P]^0", "epsilon[H]", "epsilon[W]",
               "alpha[A]","alpha[R]", 
               "f[i]", "f[o]", "f[MS]", "f[WS]", "f[US]", 
               "I[N]", "K[LM]^0",
               "K[NP]", "K[SM]^0", "q", "sigma[US]", 
               "a[M]", "a[H]", "a[WL]", "a[WM]", "a[WS]",
               "tau[H]", "tau[M]", "tau[W]", "V[HP]",
               "V[LM]^0", "V[LW]", "V[NP]^0", "V[SM]^0", "V[SW]")

png(paste0("plots_from_",Sys.Date(),"/","parameter_range_",Sys.Date(),".png"), width=12, height=4,
    units="in", res=300)
ggplot(par_plot, aes(x=NPAR, y=median), size=2) + 
  geom_pointrange(aes(ymin=lower, ymax=upper, color=Best, shape=Best, linetype=Best)) + theme_classic() +
  scale_color_manual(values=c("grey", "black")) +
  scale_shape_manual(values=c(19, 1)) + 
  geom_text(data = sigframe, aes(x=NPAR, y = 2, label=pval),size=10) +
  ylab("Scaled value") + xlab("Parameter") +
  scale_x_discrete(labels = parse(text = parameters))
dev.off()

rm(par_plot)

# check for correlations amoung the selected parameters
# low correlation amoung values selected so far... but not a huge selection gradient
png(paste0("plots_from_",Sys.Date(),"/","parameter_correlation_",Sys.Date(),".png"), width=12, height=12,
    units="in", res=600)
par2 %>% filter(Manip=="Yes" & Best=='Yes') %>% 
  select(NPAR, PARS, Run) %>% 
  spread(key=NPAR, value=PARS) %>% select(-Run) %>%
  pairs(., lower.panel = panel.cor, diag.panel = panel.hist)
dev.off()

# ....Impact of treatments on plant abundance -------------------------------

#Confirm this is working

if(verbose){
  data2 %>% as_tibble() %>% select(-PARS) %>% 
    rename(TreatmentN = Treatment) %>%
    left_join(
      data.frame(TreatmentN = seq(0,7,1),
                 Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                               "RtH"))
    ) %>% select(-TreatmentN) %>% 
    left_join(errord %>% select(Run, Best)) %>%
    filter(Best=="Yes") %>% select(-Best) %>%
    gather(-time, -Expt, -Treatment, -Run, key=StateVar, value=Size) %>%
    filter(time == 1053) %>%
    mutate(Run = as.character(Run)) %>%
    ggplot(aes(x=Treatment, y=Size, group=Run)) + 
    geom_point(color="red") + geom_line() + 
    facet_wrap(Expt~StateVar, scales="free") + theme_classic() 
}

plotdata = data2 %>% as_tibble() %>% select(-PARS, -YSTABLE) %>% 
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>% 
  left_join(errord %>% select(Run, Best)) %>%
  filter(Best=="Yes") %>% 
  select(-Best, -Expt) %>%
  gather(-time, -Treatment, -Run, key=StateVar, value=Size) %>%
  spread(key=Treatment, value=Size) %>%
  mutate(EOH1 = (H-N)/N, EOH2 = (RmH-RtH)/RtH,
         EOW1 = (W-N)/N, EOW2 = (RmW-Rt)/Rt) %>%
  mutate(IE = (HW-N - (H-N + W-N))/N) %>%
  select(-H:-W) %>%
  gather(-time, -Run, -StateVar, key=Effect, value=Size) %>%
  filter(!(StateVar %in% c("H","W", "percentStable"))) %>%
  mutate(Run = as.character(Run))

meandata =plotdata %>% group_by(time, StateVar, Effect) %>%
  summarize(lower = quantile(Size, probs = 0.945, na.rm=T),
            upper = quantile(Size, probs = 0.055, na.rm=T),
            Size = median(Size, na.rm=T))

variable_names_effect <- c(
  "EOH1" = "Effect of \n herbivores (Expt)" ,
  "EOH2" = "Effect of \n herbivores",
  "EOW1" = "Effect of \n earthworms (Expt)",
  "EOW2" = "Effect of \n earthworms",
  "IE" = "Interaction \n effect"
)

p1 = plotdata %>% filter(Effect != "IE") %>%
  ggplot(aes(x=time, y=Size)) +
  geom_line(aes(group=Run), color="grey", alpha=0.5) + 
  geom_pointrange(data= meandata%>% filter(Effect != "IE"), color="red",
                  aes(ymin=lower, ymax=upper))+
  geom_line(data= meandata %>% filter(Effect != "IE"), color="black") +
  facet_grid(StateVar~Effect, scales="free", 
             labeller=labeller(Effect = variable_names_effect,
                               StateVar = variable_names)) + 
  theme_classic()

p2 = plotdata %>% filter(Effect == "IE") %>%
  ggplot(aes(x=time, y=Size)) + ylab("") + 
    geom_line(aes(group=Run), color="grey", alpha=0.5) + 
    geom_pointrange(data= meandata%>% filter(Effect == "IE"), color="red",
                    aes(ymin=lower, ymax=upper))+
    geom_line(data= meandata%>% filter(Effect == "IE"), color="black") +
  facet_grid(StateVar~Effect, scales="free", 
             labeller=labeller(Effect = variable_names_effect,
                               StateVar = variable_names)) + 
  theme_classic()
  
png(paste0("plots_from_",Sys.Date(),"/","interaction_effect_",Sys.Date(),".png"),
    width=8, height=8, units="in", res=300)
ggpubr::ggarrange(p1, p2, widths=c(2.5,1))
dev.off()


plotdata2 = plotdata %>% 
  filter(Effect %in% c("EOH1", "EOW1", "IE") &
           StateVar %in% c("L", "M", "N", "P", "R", "S", "U")) %>%
  mutate(Year = time/365)

meandata2 =plotdata2 %>% group_by(time, StateVar, Effect) %>%
  summarize(lower = quantile(Size, probs = 0.945, na.rm=T),
            upper = quantile(Size, probs = 0.055, na.rm=T),
            Size = median(Size, na.rm=T)) %>%
  mutate(Year = time/365)

plotexpdata2 = datatomodel3 %>% 
  select(time, Treatment, StateVar, mBio) %>%
  spread(key = Treatment, value=mBio) %>%
  mutate(EOH1 = (H-N),
         EOW1 = (W-N),
         IE = (HW-N - (H-N + W-N))/N,
         EOW2 = (Rt - RmW),
         EOH2 = (RtH - RmH)) %>%
  select(time, StateVar, EOH1:EOH2) %>%
  gather(-time, -StateVar, key=Effect, value=Size) %>%
  filter(!is.na(Size)) %>% 
  filter(Effect %in% c("EOH1", "EOW1", "IE") &
           StateVar %in% c("L", "M", "N", "P", "R", "S", "U")) %>%
  mutate(Year = time/365)


p1 = plotdata2 %>% filter(Effect == "EOH1") %>%
  ggplot(aes(x=Year, y=Size)) +
  geom_line(aes(group=Run), color="grey", alpha=0.5) + 
  geom_line(data= meandata2 %>% filter(Effect == "EOH1"), color="black") +
  geom_point(data= plotexpdata2 %>% filter(Effect == "EOH1"), color="red") +
  facet_grid(StateVar~Effect, scales="free", 
             labeller=labeller(Effect = variable_names_effect,
                               StateVar = variable_names)) + 
  theme_classic()

p2 = plotdata2 %>% filter(Effect == "EOW1") %>%
  ggplot(aes(x=Year, y=Size)) +
  geom_line(aes(group=Run), color="grey", alpha=0.5) + 
  geom_line(data= meandata2 %>% filter(Effect == "EOW1"), color="black") +
  geom_point(data= plotexpdata2 %>% filter(Effect == "EOW1"), color="red") +
  facet_grid(StateVar~Effect, scales="free", 
             labeller=labeller(Effect = variable_names_effect,
                               StateVar = variable_names)) + 
  theme_classic()

p3 = plotdata2 %>% filter(Effect == "IE") %>%
  ggplot(aes(x=Year, y=Size)) + ylab("") + 
  geom_line(aes(group=Run), color="grey", alpha=0.5) + 
  geom_line(data= meandata2 %>% filter(Effect == "IE"), color="black") +
  geom_point(data= plotexpdata2 %>% filter(Effect == "IE"), color="red") +
  facet_grid(StateVar~Effect, scales="free", 
             labeller=labeller(Effect = variable_names_effect,
                               StateVar = variable_names)) + 
  theme_classic()



png(paste0("plots_from_",Sys.Date(),"/","clean_interaction_effect_",Sys.Date(),".png"),
    width=8, height=8, units="in", res=300)
ggpubr::ggarrange(p1, p2,p3, ncol=3, labels="auto", font.label = list(face="plain"))
dev.off()


# New plot of raw data over time 

rawplot1 = data2 %>% as_tibble() %>% select(-PARS, -YSTABLE) %>% 
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>% 
  left_join(errord %>% select(Run, Best)) %>%
  filter(Best=="Yes") %>% 
  select(-Best, -Expt) %>%
  gather(-time, -Treatment, -Run, key=StateVar, value=Size)

treatment_names <- c(
  "N" = "None" ,
  "H" = "Grasshoppers",
  "W" = "Earthworms",
  "HW" = "Grasshoppers and Earthworms"
)

png(paste0("plots_from_",Sys.Date(),"/","raw_temporal_",Sys.Date(),".png"),
    width=8, height=8, units="in", res=300)

rawplot1 %>% 
  mutate(Year = (time/365)-1) %>%
  filter(Year >=1) %>%
  filter(Treatment %in% c("N", "W", "H", "HW")) %>%
  filter(StateVar %in% c("H","W", "M", "N", "P")) %>%
  ggplot(aes(x=Year, y=Size)) +
  geom_line(aes(group=Run), color="grey", alpha=0.5) +
  geom_pointrange(data = ppdata2 %>%
                    mutate(Year = Year - 1) %>%
                    filter(Year >=1) %>%
                    filter(Treatment %in% c("N", "W", "H", "HW")) %>%
                    filter(StateVar %in% c("H","W", "M", "N", "P")),
                  aes(x=Year, y=mBio, ymin= lowBio, ymax = upBio), color = "red") +
  facet_grid(StateVar~Treatment, scales="free", 
             labeller=labeller(Treatment = treatment_names,
                               StateVar = variable_names)) +
  theme_classic()

dev.off() 

png(paste0("plots_from_",Sys.Date(),"/","raw_temporal_entire",Sys.Date(),".png"),
    width=8, height=8, units="in", res=300)

rawplot1 %>% 
  mutate(Year = (time/365)-1) %>%
  filter(Treatment %in% c("N", "W", "H", "HW")) %>%
  filter(StateVar %in% c("H","W", "M", "N", "P")) %>%
  ggplot(aes(x=Year, y=Size)) +
  geom_line(aes(group=Run), color="grey", alpha=0.5) +
  geom_pointrange(data = ppdata2 %>%
                    mutate(Year = Year - 1) %>%
                    filter(Treatment %in% c("N", "W", "H", "HW")) %>%
                    filter(StateVar %in% c("H","W", "M", "N", "P")),
                  aes(x=Year, y=mBio, ymin= lowBio, ymax = upBio), color = "red") +
  facet_grid(StateVar~Treatment, scales="free", 
             labeller=labeller(Treatment = treatment_names,
                               StateVar = variable_names)) +
  theme_classic()

dev.off() 
  

# Then show the measurement error associated with each pool,
# which is the residual variance in models after all variables considered

# then project forward...

# Correlate effect with parameters
png(paste0("plots_from_",Sys.Date(),"/","parameter_effect_",Sys.Date(),".png"),
    width=20, height=10, units="in", res=300)
plotdata %>% filter(Effect == "EOW1" &
                      StateVar %in% c("L", "M", "N", "P", "R", "S", "U")) %>% left_join(
  par2 %>% filter(Manip=="Yes") %>% 
    mutate(Run = as.character(Run)) %>%
    select(-fit, -Best, -Manip)
) %>%
  ggplot(aes(x=PARS, y=Size)) +
  geom_point(aes(color=time)) +
  stat_smooth(method="lm", color="red") +
  facet_grid(StateVar~NPAR, scales="free") + theme_minimal()
dev.off()
