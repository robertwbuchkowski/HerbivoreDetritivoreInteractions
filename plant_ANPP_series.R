# Try to create a long time series of plant ANPP

library(tidyverse)


read_csv("Data/schmitz_data_repo.csv") %>%
  mutate(gm2 = `Mean PB (g)`/`Area (m2)`) %>%
  select(`Experiment Year`, Treatment, Type, gm2) %>%
  rename(Year = `Experiment Year`) %>%
  bind_rows(
    read_csv("Data/adamhouston_data.csv") %>%
      filter(N.Treatment=="None" & Trophic.Treatment < 3 & !is.na(totmass) & Field == "Xmas") %>%
      select(Cage.ID, smass, gmass, fmass, Trophic.Treatment) %>%
      mutate(Grass = gmass, Herb = smass + fmass) %>%
      select(Cage.ID, Trophic.Treatment, Grass, Herb) %>%
      gather(-Cage.ID, -Trophic.Treatment, key = Type, value = Total) %>%
      group_by(Trophic.Treatment, Type) %>%
      summarize(Total = mean(Total)) %>%
      mutate(Treatment = ifelse(Trophic.Treatment ==1, "Ctrl", "MEFE")) %>%
      mutate(Year = 2017) %>%
      mutate(Total2 = Total/0.75) %>%
      ungroup() %>%
      select(Year, Type, Treatment, Total2) %>%
      bind_rows(
        read_csv("Data/zackmiller_data.csv") %>%
          filter(Env == "C" & Trophic <3) %>%
          mutate( Grass = Gr_Biomass, 
                  Herb = Sol_Biomass + Dsol_Biomass + Etc_Biomass) %>%
          select(ID,Trophic, Grass, Herb) %>%
          gather(-ID, -Trophic, key = Type, value = Total) %>%
          group_by(Trophic, Type) %>%
          summarize(Total = mean(Total)) %>%
          mutate(Year = 2016) %>% # CHECK!! 
          mutate(Treatment = ifelse(Trophic ==1, "Ctrl", "MEFE")) %>%
          mutate(Total2 = Total/0.75) %>%
          ungroup() %>%
          select(Year, Type, Treatment, Total2)
      ) %>%
      rename(gm2 = Total2)
  ) %>%
  spread(key = Type, value = gm2) %>%
  mutate(Total2 = Grass + Herb) %>%
  mutate(Total2 = ifelse(is.na(Total2), Total, Total2)) %>%
  select(Year, Treatment, Total2) %>%
ggplot(aes(x = Year, y = Total2, color = Treatment)) + geom_line() + 
  geom_point() + theme_classic()






