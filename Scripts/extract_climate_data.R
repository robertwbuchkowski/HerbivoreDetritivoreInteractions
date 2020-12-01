# Extract climate data from NOAA product:

library("tabulizer") # Version 0.2.2

# This script loads in the climate data from a PDF provided by NOAA and looks at relationships

# Get the page area I need to extract from
locate_areas("Data/Climate_data_West_Thompson_Lake.pdf", pages = 1)

# Read the tables in the Viewer: this requires manual work
Climdata <- extract_tables("Data/Climate_data_West_Thompson_Lake.pdf",
                           area = list(c(217,3.4,680,603.5))); 

# Fix an error in reading
Climdata[[91]] = Climdata[[91]][c(12:44),c(1:15)]
Climdata[[45]] = Climdata[[45]][c(13:43),c(1:15)]

# Bind data
Climdata2 = do.call("rbind", Climdata)

# Clean up data frame
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


# Run the model to predict temperature based on the time of year
m1 = nls(Temp ~ a*cos(2*3.14/365*Time + c) + b, start=c(a=-12.8244, b=281.15, c= -0.37), data=Climdata3o4)

coefficients(m1)

Climdata3o4[,"zPredicted"] = coefficients(m1)["a"]*cos(2*3.14/365*Climdata3o4$Time + coefficients(m1)["c"]) + coefficients(m1)["b"]

Climdata3o4 %>% select(Time, Temp, zPredicted) %>%
  gather(-Time, key = Type, value = Temp) %>% ggplot(aes(x=Time, y = Temp, color = Type)) + geom_line(lwd = 0.5) + theme_classic() + geom_hline(yintercept = 273.15, lty = 2)


Climdata3o4 %>% filter(Time %in% seq(1,8815, by = 30)) %>% ggplot(aes(x=Time, y = Temp)) + geom_line() + theme_classic() + geom_hline(yintercept = 273.15, lty = 2)

yearavg = Climdata3o4 %>% as_tibble() %>%
  filter(month(Date) > 3 | month(Date) < 10) %>%
  group_by(Year) %>%
  summarize(Temp = mean(Temp))