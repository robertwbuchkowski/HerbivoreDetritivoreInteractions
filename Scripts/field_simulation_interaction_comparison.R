# Calculate distribution of feedback size ----

require(tidyverse) # version 1.2.1

# The following analysis calculates the interaction effect size relative to the herbivore effect size and detritvore effect size used in Figure 6.

#create directory for plots if necessary that sets a specific date
if(!dir.exists(paste0("modelresults_",Sys.Date()))){dir.create(paste0("modelresults_",Sys.Date()))}

# A function for displaying scientific notation
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

scientific3 <- function(x){
  x = signif(x, digits = 1)
  ifelse(x==0, "0", gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
}

# Create a dictionary to rename state variables
variable_names <- c(
  "H" = "Grasshopper" ,
  "M" = "Microbial biomass",
  "N" = "Inorganic N",
  "P" = "Plant biomass",
  "W" = "Earthworm",
  "L" = "Litter",
  "S" = "Soil organic matter"
)

# Load in the simulated data and the selected runs for best fit
data2 <- readRDS("Data/fullmodeloutput_2020-06-22.rds")
BEST <- read_csv("Data/selected_runs_2020-06-22.csv", 
                 col_types = cols(Run = col_character()))

# Clean up the data
data3 = data2 %>% 
  select(-PARS, -YSTABLE, -Run2) %>% # clean out extra data: parameters, stability
  filter(time!=-1) %>% # get ride of extra QCQA returns
  rename(TreatmentN = Treatment) %>%
  left_join(
    data.frame(TreatmentN = seq(0,7,1),
               Treatment = c("N", "W", "H", "HW", "RmW", "Rt", "RmH",
                             "RtH"))
  ) %>% select(-TreatmentN) %>%
  filter(Treatment %in% c("N", "W", "H", "HW")) %>% 
  select(-Expt)

runID = data3 %>% 
  select(-time, -Treatment) %>% 
  gather(-Run, key = StateVar, value = Biomass) %>%
  mutate(Biomass = ifelse(Biomass < 0, 0,1)) %>%
  group_by(Run) %>%
  summarize(Total = sum(Biomass)) %>% 
  filter(Total == 560) %>% select(Run) %>%
  filter(Run != "9481143743947965") # Get rid of problem Run!

# Calcaulte the effects (IE = interaction; WE = detritivore; HE = herbivore). This only includes the "best" runs:
(MAXTIME = max(data3$time)) # Get the maximum time (i.e. 4-years)

# Load in the data and identify the "best" data
out1 = runID %>%
  left_join(
    data3 %>% select(-H, -W) %>% 
      filter(time == MAXTIME) %>% 
      select(-time) 
  ) %>%
  gather(-Treatment, -Run, key = StateVar, value = Biomass) %>%
  spread(key = Treatment, value = Biomass) %>% 
  mutate(IE = (HW - H - W + N)) %>%
  mutate(WE = (W - N), HE = (H - N)) %>% 
  mutate(IEpred = WE + HE, IEacc = HW - N) %>%
  select(Run, StateVar, IE, WE, HE, IEacc, IEpred) %>%
  left_join(BEST) %>%
  select(-fit) %>%
  mutate(Best = ifelse(is.na(Best),"No", Best))

# Plot for the best fitting parameter sets only
png(paste0("modelresults_",Sys.Date(),"/interactioneffect.png"), width = 8, height = 5, units = "in", res = 600)
out1 %>% filter(Best == "Yes") %>% select(-Best, -IEacc, -IEpred) %>% gather(-Run, -StateVar, key = Effect, value = value) %>%
  mutate(value = 100*abs(value) + 1e-6) %>%
  ggplot(aes(x = Effect, y = value, fill = Effect)) + geom_boxplot(alpha = 0.7) + theme_classic() + 
  facet_wrap(.~StateVar,labeller=labeller(StateVar = variable_names)) +
  scale_y_log10(name = "Effect magnitude", labels = scientific) + 
  scale_x_discrete(limits = c("HE", "WE", "IE")) + 
  scale_fill_manual(values = c("blue", "orange", "grey"), limits = c("HE", "WE", "IE"), labels = c("Grasshopper (HE)", "Earthworm (WE)", "Interaction (IE)"), name = "Effect") +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.box = "horizontal")
dev.off()

png(paste0("modelresults_",Sys.Date(),"/interactioneffect_all.png"), width = 8, height = 5, units = "in", res = 600)
out1 %>% select(-Best, -IEacc, -IEpred) %>% gather(-Run, -StateVar, key = Effect, value = value) %>%
  mutate(value = 100*abs(value) + 1e-6) %>%
  ggplot(aes(x = Effect, y = value, fill = Effect)) + geom_boxplot(alpha = 0.7) + theme_classic() + 
  facet_wrap(.~StateVar,labeller=labeller(StateVar = variable_names)) +
  scale_y_log10(name = "Effect magnitude", labels = scientific) + 
  scale_x_discrete(limits = c("HE", "WE", "IE")) + 
  scale_fill_manual(values = c("blue", "orange", "grey"), limits = c("HE", "WE", "IE"), labels = c("Grasshopper (HE)", "Earthworm (WE)", "Interaction (IE)"), name = "Effect") +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.box = "horizontal")
dev.off()

# Test whether interaction effect is always smaller than the individual effects of herbivores and detritivores at the end of the simulation
test1 = out1 %>% filter(Best == "Yes") %>% mutate(IW = -abs(IE) + abs(WE), IH = -abs(IE) + abs(HE)) %>%
  mutate(IW = IW > 0, IH = IH >0)

table(test1$IW,test1$IH)

test1 = out1 %>% mutate(IW = -abs(IE) + abs(WE), IH = -abs(IE) + abs(HE)) %>% filter(abs(HE) + abs(WE) != 0) %>%
  mutate(IW = IW > 0, IH = IH >0)

table(test1$IW,test1$IH)
