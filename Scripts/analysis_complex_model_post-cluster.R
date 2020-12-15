# Analysis of Complex model data from the cluster
# Robert W. Buchkowski
# June 23 2020

# The main purpose of this analysis is to fit the model data to the empirical data and summarize the model results. Several graphs are produced that are not part of the final paper.

# This produces the fit results presented in Table 1.

require(tidyverse) # version 1.2.1
verbose = F # Do you want extra graphs and analyses?
loadcsv = F # Do you need to load the raw inputs from a cluster, if true you will need to run new model simulations, otherwise you can just load the ones provided.

# Load in the data from field experiments -----

datatomodel2 = read_csv("Data/datatomodel2.csv")

# Load (Cluster) data -----

if(loadcsv){ # Only important if you are loading in data from the raw outputs of "complex_model_non-eqm_to_cluster.R" It would need to be modified for the specific directory on your machine
  
  dirtoload = "Model_noneqm_Jun2020/"
  
  ftoload = list.files(dirtoload)
  
  listOfDataFrames <- vector(mode = "list", length = length(ftoload))
  
  for(ii in 1:length(ftoload)){
    listOfDataFrames[[ii]] <- read.csv(paste0(dirtoload,ftoload[ii]))
  }
  
  data2 <- do.call("rbind", listOfDataFrames)
  
  rm(listOfDataFrames,ftoload)
  
  
  data2$Run = paste0(data2$Run, data2$Run2)
  
  # check for duplicate Run names and replace them
  table(data2$Run)[table(data2$Run) != 160]
  
  # Get rid of RUNS without full list
  data2 = data2 %>% filter(!(Run %in% c("138385645710656", "2071216147931349", "5029548420248414", "5777566158453009", "7718466888195712")))

  saveRDS(data2, file = paste0("Data/fullmodeloutput_",Sys.Date(),".rds"))
}

# The following .rds file was created by using the function "rbind" to join individual model outputs
data2 <- readRDS("Data/fullmodeloutput_2020-06-22.rds")
# 277,499 unique runs with full data

# Baseline parameter values, loaded for use in comparisons
source("Scripts/parameters.R"); rm(yint)

# get state variable values when they are at a stable equilibirum before the simulations start
if(verbose){
  stable1 <- data2 %>% select(YSTABLE, Run) %>% filter(YSTABLE!=-2) %>% 
    mutate(NVAR = rep(colnames(data2)[2:8], length(unique(data2$Run))))
}

# get model data to compare with empirical data
data3 <- data2 %>% select(-PARS, -YSTABLE, -Run2) %>% filter(time!=-1) %>% rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>%
  select(time, Treatment, Expt,Run, W,P,H,M,N)

PARMS = data2 %>% select(PARS, Run) %>% 
  mutate(NPAR = rep(c(names(params), rep("REMOVE", 127)), 
                    length(unique(data2$Run)))) %>%
  filter(NPAR != "REMOVE")

saveRDS(PARMS, file = paste0("Data/PARMS_",Sys.Date(),".rds"))

# Clean out memory
rm(data2, PARMS)
gc()

# Match and select best runs ----

# calculate median absolute deviation (instead of standard deviation) for unique model runs
# based on van der Vaart et al. 2015

# The full dataset is too big to run, so I need to do it in parts. I will work with it by state variable
write_rds(data3)

MSD = data3 %>% group_by(Expt) %>%
  summarise(W = mad(W),
            P = mad(P),
            H = mad(H),
            M = mad(M),
            N = mad(N))

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
  
(Nobs = length(unique(data3$Run))) # Number of runs
(cut = (50/Nobs)) # Cut off to get 50 best models

# Create dataframe to hold the error data

errord = tibble(Run = unique(data3$Run))

# List of state variables
SV = c("W", "P", "H", "M", "N")

# Get the individual fit for each state variable in a loop (because the full dataset is too big) and assign it the name errord + state variable letter
for(i in 1:5){
  assign(paste0("errord", SV[i]), 
         data3[,c("Run","time", "Treatment","Expt",SV[i])] %>%
           rename(Model = SV[i]) %>%
           right_join(
             datatomodel3 %>% filter(StateVar == SV[i]) %>% select(time, Treatment, mBio)
           ) %>%
           left_join(
             MSD %>% select(Expt, SV[i]) %>% rename(Msd = SV[i])
           ) %>%
           filter(!is.na(Model)) %>% # remove any cases without simulation
           filter(Msd !=0) %>% # remove hopper treatments where standard deviation is zero by definition
           group_by(Run) %>% # create group by Run so standardize function works
           summarize(fit = sqrt(sum(((Model-mBio)/Msd)^2)))
         )
}

# Combine the different errors into one file and clean out the separate ones
errord = errordH %>% rename(H = fit) %>%
  full_join(
    errordW %>% rename(W = fit)
  ) %>%
  full_join(
    errordP %>% rename(P = fit)
  ) %>%
  full_join(
    errordN %>% rename(N = fit)
  ) %>%
  full_join(
    errordM %>% rename(M = fit)
  )

rm(errordH, errordM, errordN, errordP, errordW)

# Find the best fit by summing across all state variables: remember this is OK becuase they are standardized above using Msd. The quantile of the fit was calculated above as approximately the best fitting (i.e. lowest fit) 0.02%
errord = errord %>% mutate(fit = H + W + P + N + M) %>%
  select(Run, fit) %>%
  mutate(Best = ifelse(fit < quantile(fit, cut), "Yes", "No"))

# Confirm that there are only 50 Best runs
table(errord %>% select(Best)) 

# Save ID of selected Runs if you want
if(verbose){
  errord %>% filter(Best == "Yes") %>% write_csv(paste0("Data/selected_runs_", Sys.Date(),".csv"))
}

# calculate R2 by treatment using the average model runs of best models and save it in data5
bestruns = errord %>% filter(Best=="Yes") %>% select(Run) %>% pull()

data5 = data3[data3$Run %in% bestruns,]
data5 = as_tibble(data5) %>%
  select(-Expt) %>%
  gather(-time, -Treatment, -Run, key = StateVar, value = Model)

(databest = data5  %>%
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

# Create a conversion list to turn state variable letters into full names to create nicer plots
variable_names <- c(
  "H" = "Herbivore" ,
  "M" = "Microbial",
  "N" = "Inorganic N",
  "P" = "Plant",
  "W" = "Earthworm",
  "L" = "Litter",
  "S" = "Soil"
)

variable_names_Expt <- c(
  "1" = "Other field data",
  "2" = "Cage experiment"
)

# Plot results of the parameter selection -----------------------------------------------

#create directory for plots if necessary that sets a specific date
if(!dir.exists(paste0("modelresults_",Sys.Date()))){dir.create(paste0("modelresults_",Sys.Date()))}

# Plot the comparison between model and the empirical data: These plots have model on the x-axis and field data on the y-axis. Much of this code is summarizing and reorganizing the two data sources.
parplot1 = data5 %>% 
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

modelexptplot = parplot1 %>%
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

png(paste0("modelresults_",Sys.Date(),"/","model_expt_", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
modelexptplot
dev.off()

rm(parplot1, parplot2)

# Look at state variable distribution for selected runs
statevardist = data5 %>% filter(time ==1395) %>%
  ggplot(aes(x=Model, fill = Treatment)) + 
  geom_histogram() + facet_wrap(.~StateVar, scales="free",labeller=labeller(StateVar = variable_names)) + theme_classic() +
  scale_fill_discrete(breaks=c("N", "H", "W", "HW",
                               "Rt", "RmW", "RtH", "RmH"),
                      labels=c("None (Expt)", "Herbivore (Expt)",
                               "Earthworm (Expt)",
                               "Both (Expt)", 
                               "Worm control",
                               "Worm removal", 
                               "Herbivore control",
                               "Herbivore removal"))

png(paste0("modelresults_",Sys.Date(),"/","state_var_dist_", Sys.Date(),".png"), width=8,
    height=5, units="in", res=300)
statevardist
dev.off()

# Look at parameter matches
parVT = read_rds("Data/PARMS_2020-06-22.rds") %>% left_join(errord) %>% 
  filter(!NPAR %in% c("Kappa","Tref_W", "Tref_P") & !is.na(fit))

# Compare the variance for the best parameters and the rest using the 'var.test' function. This analysis is used to determine whether the 50 best fitting parameter sets actually fit better than the remaining ones.

pvec = unique(parVT$NPAR) # vector for each variable
sigvec = rep(-1, length(pvec)) # vector for significance
for(i in 1:length(pvec)){
  vec1 = subset(parVT,NPAR == pvec[i] & Best=="No")$PARS
  vec2 = subset(parVT,NPAR == pvec[i] & Best=="Yes")$PARS
  sigvec[i] = var.test(vec1,vec2, ratio=1, alternative = "greater")$p.value # significance tested with var.test function
  
}

rm(vec1, vec2)

# Create a plot of the results:
sigframe = data.frame(NPAR = pvec, 
                      pval = ifelse(sigvec < 0.05, "*", ""))

# Nice dataframe for plot
par_plot = parVT %>% 
  left_join(data.frame(NPAR = names(params),litvalue = unname(params))) %>% # add literature values
  mutate(StdValue = PARS/litvalue) %>% # standardize estimates for easy display
  select(NPAR, Best, StdValue) %>% group_by(NPAR, Best) %>%
  summarize(lower = quantile(StdValue, 0.055), 
            upper = quantile(StdValue, 0.945),
            median = quantile(StdValue, 0.5))

# Parameter names
parameters = c("B", "D", "E[a]",
               "f[i]", "f[o]", 
               "I[N]",
               "K[int]", "K[LM]^0", 
               "K[NP]","K[slope]", "K[SM]^0",
               "l", "q",
               "a[M]", "a[H]", "a[WL]", "a[WM]", "a[WS]",
               "tau[H]", "tau[M]", "alpha","tau[W]", 
               "V[HP]", "V[int]",
               "V[LM]^0", "V[LW]", "V[NP]",
               "V[slope]","V[SM]^0", "V[SW]")

# Make the plot
parameterrange = ggplot(par_plot, aes(x=NPAR, y=median), size=2) + 
  geom_pointrange(aes(ymin=lower, ymax=upper, color=Best, shape=Best, linetype=Best)) + theme_classic() +
  scale_color_manual(values=c("grey", "black")) +
  scale_shape_manual(values=c(19, 1)) + 
  geom_text(data = sigframe, aes(x=NPAR, y = 2, label=pval),size=10) +
  ylab("Scaled value") + xlab("Parameter") +
  scale_x_discrete(labels = parse(text = parameters))

png(paste0("modelresults_",Sys.Date(),"/","parameter_range_",Sys.Date(),".png"), width=12, height=4,
    units="in", res=300)
parameterrange
dev.off()
