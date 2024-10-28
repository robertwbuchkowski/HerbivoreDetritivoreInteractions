# New model to fit the climate data for ESW:

#Library
library(tidyverse)

#Load in the climate data:

cd <- read_csv("Data/climate-daily.csv") %>%
  select(MEAN_TEMPERATURE, LOCAL_DATE, TOTAL_PRECIPITATION)%>%
  bind_rows(
    read_csv("Data/climate-daily (1).csv") %>%
      select(MEAN_TEMPERATURE, LOCAL_DATE,TOTAL_PRECIPITATION)
  ) %>%
  bind_rows(
    read_csv("Data/climate-daily (2).csv") %>%
      select(MEAN_TEMPERATURE, LOCAL_DATE,TOTAL_PRECIPITATION)
  )

cd %>% ggplot(aes(x = LOCAL_DATE, y = MEAN_TEMPERATURE)) + geom_point()

cd %>%
  filter(LOCAL_DATE > "1980-01-01") %>% ggplot(aes(x = LOCAL_DATE, y = MEAN_TEMPERATURE)) + geom_point()

# Filter to only 1980 forwards:

cd = cd %>%
  filter(LOCAL_DATE > "1997-04-14") %>%
  mutate(doy = lubridate::yday(LOCAL_DATE)) %>%
  reframe(MEAN_TEMPERATURE = mean(MEAN_TEMPERATURE, na.rm=  T), TOTAL_PRECIPITATION = mean(TOTAL_PRECIPITATION, na.rm=  T), .by = doy)

cd %>% ggplot(aes(x = doy, y = MEAN_TEMPERATURE)) + geom_point()

cd = cd %>% filter(doy != 366) %>%
  mutate(MEAN_TEMPERATURE = MEAN_TEMPERATURE + 273.15)

# Function for Connecticut:
# -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846

m1 = nls(MEAN_TEMPERATURE~ -a*cos(2*3.14/365*doy - b) + c, data = cd, start = c(a = 12.8244, b = 0.3666, c = 281.9846))

m1

LTtemp = function(doy){
  
  -12.7084*cos(2*3.14/365*doy-0.3758)+281.8291
  
}

cd %>%
  mutate(tf = LTtemp(doy)) %>% ggplot(aes(x = doy, y = MEAN_TEMPERATURE)) + geom_point() + geom_line(aes(y = tf), color = "red")
