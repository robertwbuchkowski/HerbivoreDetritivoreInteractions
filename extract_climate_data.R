# Extract climate data from NOAA product

library("tabulizer")

# Get the page area I need to extract from
locate_areas("Data/Climate_data_West_Thompson_Lake.pdf", pages = 1)

Climdata <- extract_tables("Data/Climate_data_West_Thompson_Lake.pdf",
                           area = list(c(217,3.4,680,603.5))); 

# Fix an error in reading
Climdata[[91]] = Climdata[[91]][c(12:44),c(1:15)]
Climdata[[45]] = Climdata[[45]][c(13:43),c(1:15)]

Climdata2 = do.call("rbind", Climdata)

View(Climdata2)

Climdata3 = Climdata2[,c(1:11)]

colnames(Climdata3) = c("Year", "Month", "Day", "Tmin", "Tmax", "Temp", "Rain", "FlagRain", "Snow", "FlagSnow", "GroundSnow")


Climdata3 = as.data.frame(Climdata3)

Climdata3o1 = Climdata3[!(Climdata3$Year == "Summary" | Climdata3$Tmin == "Summary"),]

Climdata3o1[Climdata3o1 == "T"] = NA
Climdata3o1[Climdata3o1 == ""] = NA

Climdata3o1 %>% as_tibble()

Climdata3o1$Year

fixdata <- function(vvv){
  as.numeric(levels(vvv))[vvv]
}

Climdata3o2 = sapply(Climdata3o1, fixdata) %>% as.data.frame()

Climdata3o2 = Climdata3o2[!(is.na(Climdata3o2$Year) | is.na(Climdata3o2$Day)),1:6]

Climdata3o3 = data.frame(Date = ymd(paste0(Climdata3o2$Year, "-",str_pad(as.character(Climdata3o2$Month), 2, "left", "0"), "-",str_pad(as.character(Climdata3o2$Day), 2, "left", "0"))),
           Temp = (Climdata3o2$Temp-32)*5/9  + 273.15)

Climdata3o4 = subset(Climdata3o3, !is.na(Temp))

Climdata3o4[,"Year"] = year(Climdata3o4$Date)
class(Climdata3o4$Date)

Climdata3o4[,"Time"] = as.numeric(difftime(Climdata3o4$Date,as.Date("1996-01-01"), units = "days"))


m1 = nls(Temp ~ a*cos(2*3.14/365*Time + c) + b, start=c(a=-12.8244, b=281.15, c= -0.37), data=Climdata3o4)

coefficients(m1)

Climdata3o4[,"zPredicted"] = coefficients(m1)["a"]*cos(2*3.14/365*Climdata3o4$Time + coefficients(m1)["c"]) + coefficients(m1)["b"]

Climdata3o4 %>% select(Time, Temp, zPredicted) %>%
  gather(-Time, key = Type, value = Temp) %>% ggplot(aes(x=Time, y = Temp, color = Type)) + geom_line(lwd = 0.5) + theme_classic() + geom_hline(yintercept = 273.15, lty = 2)


range(Climdata3o4$Time)

Climdata3o4 %>% filter(Time %in% seq(1,8815, by = 30)) %>% ggplot(aes(x=Time, y = Temp)) + geom_line() + theme_classic() + geom_hline(yintercept = 273.15, lty = 2)

yearavg = Climdata3o4 %>% as_tibble() %>%
  filter(month(Date) > 3 | month(Date) < 10) %>%
  group_by(Year) %>%
  summarize(Temp = mean(Temp))

# Try to create a long time series of plant ANPP ----

library(tidyverse)


ANPP = read_csv("Data/schmitz_data_repo.csv") %>%
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
  select(Year, Treatment, Total2)

ANPP %>% ggplot(aes(x = Year, y = Total2, color = Treatment)) + geom_line() + 
  geom_point() + theme_classic()

ANPPtemp = ANPP %>% filter(Treatment == "Ctrl") %>%
  select(-Treatment) %>%
  left_join(yearavg) %>%
  mutate(Temp2 = Temp*Temp) %>%
  mutate(Tempexp = exp(1/Temp))

ANPPtemp %>%
  ggplot(aes(x = Temp, y = Total2)) + geom_point() + theme_classic()

summary(lm(Total2 ~ Temp, data = ANPPtemp))

# No good relationship over time...sadness :(
