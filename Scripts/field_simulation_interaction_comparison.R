# Calculate distribution of feedback size ----

require(tidyverse) # version 1.2.1

# The following analysis produces the complex model data used in the simple_model.R script. It also calculates the interaction effect size relative to the herbivore effect size and detritvore effect size used in Figure 6.

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

out2 = out1 %>%
  left_join(
    data.frame(StateVar = c("P", "L", "N", "S", "M"),
               ncol = as.character(c("#009E73","#E69F00","#56B4E9","#F0E442", "#0072B2")),
               npch = c(1,2,3,4,5))
  ) %>%
  left_join(
    data.frame(Best = c("Yes", "No"),
               ncex = c(1, 0.5))
  )

out3 = out2[1:(5*10000),] # Take the first 10,000 model runs to compare with the simple model

# This is the data that is used in the simple_model.R script to plot the comparison
write_rds(out3, "Data/TrueLinearComplex.rds")

# Go back and get the state variable sizes for the selected 10,000 model runs (line 89)
out3 %>% select(Run) %>% 
  left_join(
    data3 %>% 
      filter(time == MAXTIME)
  ) %>% 
  select(-time) %>%
  write_rds("Data/TrueLinearComplex_StateVar.rds")

rm(data3)
gc()

# Plot the full interaction effects available the grasshopper and earthworm effects

# Plot for the best fitting parameter sets only
png(paste0("modelresults_",Sys.Date(),"/interactioneffect.png"), width = 8, height = 5, units = "in", res = 600)
out1 %>% filter(Best == "Yes") %>% select(-Best, -IEacc, -IEpred) %>% gather(-Run, -StateVar, key = Effect, value = value) %>%
  mutate(value = 100*abs(value) + 1e-6) %>%
  ggplot(aes(x = value, fill = Effect)) + geom_density(alpha = 0.7) + theme_classic() + 
  facet_wrap(.~StateVar,labeller=labeller(StateVar = variable_names)) +
  scale_x_log10(name = "Effect", labels = scientific) + 
  scale_fill_manual(values = c("blue", "orange", "grey"), limits = c("HE", "WE", "IE"), labels = c("Grasshopper", "Earthworm", "Interaction"), name = "Effect") +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.box = "horizontal")
dev.off()

# Build a data frame for the plant plot: Figure 6
plantplot = out1 %>% 
  filter(Best == "Yes") %>%
  select(-Best, -IEacc, -IEpred) %>%
  gather(-Run, -StateVar, key = Effect, value = value) %>%
  mutate(value = abs(value) + 1e-6) %>%
  filter(StateVar == "P")

# Plot the interaction effect for on plants
png(paste0("modelresults_",Sys.Date(),"/interactioneffect_Figure6.png"), width = 5, height = 5, units = "in", res = 600)
gt =  plantplot %>%
  group_by(Effect) %>%
  summarise(X = median(value)) %>%
  mutate(Y = c(0.65, 0.65, 0.65)) %>%
  mutate(t = signif(X, 2)) %>%
  mutate(X = X*c(0.1, 1,1)) %>%
  mutate(t = scientific3(t))

plantplot %>%
  ggplot(aes(x = value)) +
  geom_text(aes(x = X, y = Y, label = t, col = Effect),data = gt, parse = T) + 
  geom_density(aes(fill = Effect),alpha = 0.7) + theme_classic() + 
  scale_x_log10(name = "Interaction effect onto plants", labels = scientific) +
  scale_fill_manual(values = c("blue", "orange", "grey"), limits = c("HE", "WE", "IE"), labels = c("Herbivore", "Detritivore", "Interaction"), name = "Effect") +
  scale_color_manual(values = c("blue", "orange", "grey"), limits = c("HE", "WE", "IE"), labels = c("Herbivore", "Detritivore", "Interaction"), name = "Effect") +
  theme(legend.position = c(0.25, 0.75),
        legend.justification = c(1, 0),
        legend.box = "horizontal") + ylab("Density")
dev.off()


# Test whether interaction effect is always smaller than the individual effects of herbivores and detritivores at the end of the simulation
test1 = out1 %>% filter(Best == "Yes") %>% mutate(IW = -abs(IE) + abs(WE), IH = -abs(IE) + abs(HE)) %>%
  mutate(IW = IW > 0, IH = IH >0)

table(test1$IW,test1$IH)

test1 = out1 %>% mutate(IW = -abs(IE) + abs(WE), IH = -abs(IE) + abs(HE)) %>% filter(abs(HE) + abs(WE) != 0) %>%
  mutate(IW = IW > 0, IH = IH >0)

table(test1$IW,test1$IH)