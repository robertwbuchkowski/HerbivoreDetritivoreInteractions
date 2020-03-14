# Predict the next time step approach

# Run the model with annual time steps (fall to fall)

# Consider each "cage" as a separate predictor of the gradient change given each state variable starting value

# then fit the model to the correct gradient

library(tidyverse)

datatomodel2 = readRDS("Data/datatomodel2.Rds")

sort(unique(datatomodel2$time))

datatomodel2 %>%
  left_join(
    unique(rbind(wormavgsample, soilavgsample,wormWEavgsample))
  ) %>%
  separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>%
  filter(Season == "Fall") %>%
  filter(Treatment %in% c("N", "H", "W", "HW")) %>%
  ggplot(aes(x=Year, y = Biomass, group = Plot, color= Treatment)) +
  geom_line() + facet_wrap(.~StateVar, scales = "free") +
  theme_classic()

# Model to fit
# dP = apx*P*X - tp*P - ahp*H*P

dataplant = datatomodel2 %>%
  left_join(
    unique(rbind(wormavgsample, soilavgsample,wormWEavgsample))
  ) %>%
  separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>%
  filter(Season == "Fall") %>%
  filter(Treatment %in% c("N", "H", "W", "HW")) %>%
  filter(Year == "17") %>%
  select(Plot, StateVar, Biomass) %>%
  spread(key = StateVar, value = Biomass) %>%
  left_join(
    datatomodel2 %>%
      left_join(
        unique(rbind(wormavgsample, soilavgsample,wormWEavgsample))
      ) %>%
      separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>%
      filter(Season == "Fall") %>%
      filter(Treatment %in% c("N", "H", "W", "HW")) %>%
      filter(Year == "18") %>%
      select(Plot, StateVar, Biomass) %>%
      spread(key = StateVar, value = Biomass) %>%
      select(Plot, P) %>%
      rename(P18 = P)
  ) %>% rename(X = N) %>%
  mutate(dP = P18 - P)

m1 = nls(dP ~ apx*P*X - tp*P - ahp*H*P, data = dataplant,
         start = c(apx = 0.1, tp = 0.1, ahp = 0.1))

summary(m1)

# bad fit, maybe because worms directly effect the population size of plants, but N effects don't necessarily cover that part
# bad fit, because low N levels means more biomass

plot(P~X, dataplant)

dataworm = datatomodel2 %>%
  left_join(
    unique(rbind(wormavgsample, soilavgsample,wormWEavgsample))
  ) %>%
  separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>%
  filter(Season == "Fall") %>%
  filter(Treatment %in% c("N", "H", "W", "HW")) %>%
  filter(Year == "16") %>%
  select(Plot, StateVar, Biomass) %>%
  spread(key = StateVar, value = Biomass) %>%
  left_join(
    datatomodel2 %>%
      left_join(
        unique(rbind(wormavgsample, soilavgsample,wormWEavgsample))
      ) %>%
      separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>%
      filter(Season == "Fall") %>%
      filter(Treatment %in% c("N", "H", "W", "HW")) %>%
      filter(Year == "17") %>%
      select(Plot, StateVar, Biomass) %>%
      spread(key = StateVar, value = Biomass) %>%
      select(Plot, W) %>%
      rename(W18 = W)
  ) %>%
  mutate(dW = W18 - W)

m1 = nls(dW ~ rw*W - rK*W*W, data = dataworm,
         start = c(rw = 0.1, rK = 14))

summary(m1)

wormfucn <- function(t,y,pars){
  
  with(as.list(c(pars,y)),{
    dW = rw*W - rK*W*W
    
    return(list(c(dW)))
  })
}

PAR = c(
  rw = 2.1056,
  rK = 0.8405
)

output = data.frame(W0 = dataworm$W,
                    Wreal = dataworm$W18,
                    Wfit = NA)

for(i in 1:60){
  output[i,"Wfit"] = output[i,"W0"] + wormfucn(t = 0, y = c(W = output[i,"W0"]), pars = PAR)[[1]]
}


plot(Wfit~Wreal, data=output)
summary(lm(Wfit~Wreal, data=output))

