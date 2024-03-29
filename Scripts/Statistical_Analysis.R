# Statistical Analysis

library(vegan) # version 2.5-6
library(nlme) # version 3.1-137
library(lme4) # version 1.1-19
library(effects) # version 4.1-0
library(piecewiseSEM) # version 2.0.2
library(tidyverse) # version 1.2.1

verbose = F

# Load in the various data sets
wormdata2 <- read_rds("Data/wormdata2.rds")
pbdata <- read_rds("Data/pbdata.rds")
RDAdata <- read_rds("Data/RDAdata.rds")
wormWE <- read_rds("Data/wormWE.rds")
RDAWdata <- read_rds("Data/RDAWdata.rds")
SOAL2017 <- read_rds("Data/SOAL2017.rds")
pbdataWE <- read_rds("Data/pbdataWE.rds")

# 1.0 Experiment Analysis ----
# .... 1.0.1 Explore CH4 worm data ----

m1 = lmer(WORM_N~ Addition + Season + Year + (1|Plot), data=wormdata2 %>% mutate(Year= as.factor(Year))); summary(m1)
summary(multcomp::glht(m1, linfct = multcomp::mcp(Addition="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Year="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Season="Tukey")))
plot(m1)
plot(resid(m1)~Addition, data=wormdata2)
boxplot(resid(m1)~Season, data=wormdata2)
boxplot(resid(m1)~Year, data=wormdata2)

m1 = lmer(AP_N~ Addition + Season + Year + (1|Plot), data=wormdata2 %>% mutate(Year= as.factor(Year))); summary(m1)
summary(multcomp::glht(m1, linfct = multcomp::mcp(Addition="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Year="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Season="Tukey")))

m1 = lmer(LUM_N~ Addition + Season + Year + (1|Plot), data=wormdata2 %>% mutate(Year= as.factor(Year))); summary(m1)
summary(multcomp::glht(m1, linfct = multcomp::mcp(Addition="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Year="Tukey")))
summary(multcomp::glht(m1, linfct = multcomp::mcp(Season="Tukey")))

# SUMMARY: LUM not well manipulated by treatments, but AP was well manipulated. LUM actually higher in removal treatments, which is consistent with their ecological role as a primary colonizer to worm-free areas. Soil temperautre, season, and treatment determine the number of worms removed. Can use SoilTavg as a covariate to models considering the full worm numbers

# Plot the earthworm data
wormdata2 %>% ggplot(aes(x=SeasonYear, y=WORM_N, color=Addition)) + geom_boxplot() + geom_jitter() + theme_classic() + scale_color_discrete(labels = c("Add", "Control", "Remove")) + scale_x_discrete(limits=c("Spring16", "Fall16","Spring17", "Fall17","Spring18", "Fall18"), labels = c("Spring '16", "Fall '16","Spring '17", "Fall '17","Spring '18", "Fall '18")) + ylab("Earthworms (#)") + xlab("Season and Year")

wormdata2 %>% select(WORM_N, AP_N, LUM_N, SoilTavg,VWCavg,CloudCover, Air_Temp) %>% gather(-WORM_N, -AP_N, -LUM_N, key=abiotic_key, value=abiotic) %>% gather(-abiotic_key, -abiotic, key=wormtype, value=Number) %>% filter(!is.na(abiotic)) %>% ggplot(aes(x=abiotic, y=Number)) + geom_point() + facet_wrap(abiotic_key~wormtype, scales="free_x") + theme_classic()


# ... 1.0.2 Correct worm data for variable extraction efficiency ----

# Create a dataframe with only fall 2017 and 2018 data
wormdata3 = wormdata2 %>% filter(SeasonYear == "Fall17"|SeasonYear == "Fall18")

ch4_correctworm0 <- lm(WORM_N ~ VWCavg + SoilTavg + CloudCover + Air_Temp, data= wormdata3); summary(ch4_correctworm0) # least significant CloudCover and Air_Temp
ch4_correctworm <- lm(WORM_N ~ VWCavg + SoilTavg, data= wormdata3); summary(ch4_correctworm) # both significant
1/(1-summary(ch4_correctworm)$r.squared)
car::vif(ch4_correctworm)

# Correct endogeic earthworm numbers
ch4_correctAP0 <- lm(AP_N ~ VWCavg + SoilTavg + CloudCover + Air_Temp, data= wormdata3); summary(ch4_correctAP0) # least significant CloudCover and Air_Temp

ch4_correctAP <- lm(log(AP_N+0.1) ~ SoilTavg, data= wormdata3); summary(ch4_correctAP) # least significant CloudCover and Air_Temp

# Create table for supplemental material
sjPlot::tab_model(ch4_correctworm0,ch4_correctworm, collapse.ci = T,pred.labels = c("Intercept", "Volumetric Water Content", "Soil Temperature","Cloud Cover (%)","Air Temperature"), dv.labels = c("Worm (#; all variables)", "Worm (#; selected variables)"),p.style="stars")

# Replace the earthworm data for the earthworm abundance from the model with extraction conditions
wormdata3 = wormdata3 %>% 
  select(Plot, Season, Year) %>%
  mutate(WORM_N_STD = as.vector(resid(ch4_correctworm)),
         AP_N_STD = as.vector(resid(ch4_correctAP))) %>%
  select(-Season)

pbdata = pbdata %>% 
  full_join(wormdata3) %>%
  rename(WORM_N_old = WORM_N, AP_N_old = AP_N) %>%
  rename(WORM_N = WORM_N_STD, AP_N = AP_N_STD)

# Replace Earthworm in the RDA datafile for running the redundancy analysis
RDAdata = RDAdata %>% rename(Earthworm_old = Earthworm, AP_N_old = AP_N) %>% 
  left_join(pbdata %>% select(Plot, Year, WORM_N, AP_N)) %>%
  rename(Earthworm = WORM_N) %>%
  mutate(Year = as.factor(Year))

pbdata %>% ggplot(aes(x = WORM_N_old, y = WORM_N, shape=Year, color = Addition)) + geom_point(size = 2) + theme_classic() + xlab("Earthworm (#)") + ylab("Earthworm (residuals)") + scale_color_manual(name = "Earthworm Treatment", values = c("#E69F00", "#56B4E9", "#009E73")) + scale_shape_discrete(labels= c("2017", "2018"))

# .... Export pbdata to RDS file for using in other code sections
write_rds(pbdata, "Data/pbdata_modified.rds")

# ... 1.0.3 Plot the outputs by treatment ----

pbdata_temp = tibble(
  Treatment = factor(c("Control", "Grasshopper", "Earthworm", "G & W", "Biomass Ctrl","Biomass Ctrl","Biomass Ctrl"), levels = c("Control", "Biomass Ctrl", "Grasshopper", "Earthworm", "G & W")),
  HopperAdd = c("Remove", "Add", "Remove", "Add", "Frozen", "Remove", "Frozen"),
  Addition = as.factor(c("Remove", "Remove", "Add", "Add", "Frozen", "Frozen", "Remove"))
) %>%
  right_join(
    pbdata, by = c("HopperAdd", "Addition")
  ) %>%
  mutate(Year = ifelse(Year == "17", "2017", "2018"))

jpeg(paste0("Figure2_",Sys.Date(), ".jpeg"), units="in", width=9, height=6, res=600)
ggpubr::ggarrange(
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = PlantBiomass, color = Year)) + geom_boxplot()  + theme_classic() + ylab("Aboveground plant biomass") + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = NTs, color = Year)) + geom_boxplot() + theme_classic() + ylab("Field N mineralization") + 
    theme(legend.position = c(0.1, 0.8),
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = NTm, color = Year)) + geom_boxplot() + theme_classic() + ylab("Lab N mineralization")  + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = SIR, color = Year)) + geom_boxplot() + theme_classic()  + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm"))
)
dev.off()

pdf("figure_2.pdf", width=9, height=6)
ggpubr::ggarrange(
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = PlantBiomass, color = Year)) + geom_boxplot()  + theme_classic() + ylab("Aboveground plant biomass") + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = NTs, color = Year)) + geom_boxplot() + theme_classic() + ylab("Field N mineralization") + 
    theme(legend.position = c(0.1, 0.8),
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = NTm, color = Year)) + geom_boxplot() + theme_classic() + ylab("Lab N mineralization")  + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm")),
  
  pbdata_temp %>%
    ggplot(aes(x = Treatment, y = SIR, color = Year)) + geom_boxplot() + theme_classic()  + 
    theme(legend.position = "none",
          axis.ticks.length=unit(-0.2, "cm"))
)
dev.off()

# .. 1.1 Linear models --------------------

# A function that runs the different linear models with and without interactions
# This function keeps the following code clean and allows me to change variable specifications without copying the code multiple times
linearmodels<- function(inputdata, inter = T){
  # Plant model
  inputdata[,"Base16"] = inputdata$PlantBiomass16
  PBinter = lmer(PlantBiomass ~ WORM_N*HopperN + Year + Base16 + (1|Plot), data= inputdata)
  PBfinal = lmer(PlantBiomass ~ WORM_N + HopperN + Year + Base16 + (1|Plot), data= inputdata)
  
  # Microbial biomass (SIR) model
  inputdata[,"Base16"] = inputdata$SIR16
  SIRinter = lmer(SIR ~ WORM_N*HopperN + Base16 + Year + (1|Plot), data= inputdata)
  SIRfinal = lmer(SIR ~ WORM_N + HopperN + Base16 + Year + (1|Plot), data= inputdata)
  
  # Field nitrogen mineralization model
  inputdata[,"Base16"] = inputdata$NTs16
  NTsinter = lmer(log(NTs) ~ WORM_N*HopperN + Base16 + Year + (1|Plot), data= inputdata, na.action = na.omit)
  NTsfinal = lmer(log(NTs) ~ WORM_N+ HopperN + Base16 + Year + (1|Plot), data= inputdata, na.action = na.omit)
  
  # Laboratory nitrogen mineralization model
  inputdata[,"Base16"] = inputdata$NTm16
  NTminter = lmer(NTm ~ WORM_N*HopperN + Base16 + Year + (1|Plot), data= inputdata)
  NTmfinal = lmer(NTm ~ WORM_N + HopperN + Base16 + Year + (1|Plot), data= inputdata)
  
  # Create table to output results
  if(inter){
    sjPlot::tab_model(PBinter,SIRinter, NTsinter, NTminter, collapse.ci = T, 
                      pred.labels = c("Intercept", "Earthworm (#)", "Grasshopper (#)","Year [2018]","Baseline '16 (of dependent variable)", "Earthworm × Grasshopper"),
                      dv.labels = c("Plant Biomass", "Substrate Induced Respiration","Field N mineralization (log)", "Lab N mineralization"), p.style="stars")
  }else{
    sjPlot::tab_model(PBfinal,SIRfinal, NTsfinal, NTmfinal, collapse.ci = T, 
                      pred.labels = c("Intercept", "Earthworm (#)", "Grasshopper (#)","Year [2018]","Baseline '16 (of dependent variable)"),
                      dv.labels = c("Plant Biomass", "Substrate Induced Respiration","Field N mineralization (log)", "Lab N mineralization"), p.style="stars")
  }

}

# ... 1.1.1 Original model analysis ----
linearmodels(pbdata) # With interaction

linearmodels(pbdata, inter= F) # Without interaction

# ... 1.1.2 Model analysis with raw earthworm number ----
linearmodels(pbdata %>% mutate(WORM_N = WORM_N_old))

# ... 1.1.3 Model analysis with high quality earthworm data ----

# Find the quantiles graphically to subset only high quality data
pbdata %>% select(Year, WORM_N_old) %>% group_by(Year) %>% summarize(lq = quantile(WORM_N_old, c(0.25)),med = quantile(WORM_N_old, c(0.5)), uq = quantile(WORM_N_old, c(0.75)))

pbHQ = pbdata %>% filter(Decimated ==0) %>% 
  filter(HopperAdd=="Add" & HopperN >0 | 
           HopperAdd!="Add" & HopperN ==0) %>% 
  filter(Addition=="Add" & Year =="17" & WORM_N_old >9 | 
           Addition!="Add" & Year =="17" & WORM_N_old <9 |
           Addition=="Add" & Year =="18" & WORM_N_old >1.5 | 
           Addition!="Add" & Year =="18" & WORM_N_old <1.5 )

linearmodels(pbHQ)

linearmodels(pbHQ %>% mutate(WORM_N = WORM_N_old))

# ... 1.1.4 Model analysis with endogeic earthworms ----
linearmodels(pbdata %>% mutate(WORM_N = AP_N))
# ... 1.1.5 Model analysis with uncorrected endogeic earthworms ----
linearmodels(pbdata %>% mutate(WORM_N = AP_N_old))
# ... 1.1.6 Model analysis of two years separately ----

# A new version of the 'linearmodels' function that runs the years separately: The variable names are the same
YearSeparateModels<- function(inputdata, inter = T){
  
  inputdata[,"Base16"] = inputdata$PlantBiomass16
  PBinter = lm(PlantBiomass ~ WORM_N*HopperN + Base16 , data= inputdata)
  PBfinal = lm(PlantBiomass ~ WORM_N + HopperN + Base16 , data= inputdata)
  
  inputdata[,"Base16"] = inputdata$SIR16
  SIRinter = lm(SIR ~ WORM_N*HopperN + Base16 , data= inputdata)
  SIRfinal = lm(SIR ~ WORM_N + HopperN + Base16 , data= inputdata)
  
  inputdata[,"Base16"] = inputdata$NTs16
  NTsinter = lm(log(NTs) ~ WORM_N*HopperN + Base16 , data= inputdata, na.action = na.omit)
  NTsfinal = lm(log(NTs) ~ WORM_N+ HopperN + Base16 , data= inputdata, na.action = na.omit)
  
  inputdata[,"Base16"] = inputdata$NTm16
  NTminter = lm(NTm ~ WORM_N*HopperN + Base16 , data= inputdata)
  NTmfinal = lm(NTm ~ WORM_N + HopperN + Base16 , data= inputdata)
  
  # Create table to output results
  if(inter){
    sjPlot::tab_model(PBinter,SIRinter, NTsinter, NTminter, collapse.ci = T, 
                      pred.labels = c("Intercept", "Earthworm (#)", "Grasshopper (#)", "Baseline '16 (of dependent variable)","Earthworm × Grasshopper"),
                      dv.labels = c("Plant Biomass", "Substrate Induced Respiration","Field N mineralization (log)", "Lab N mineralization"), p.style="stars")
  }else{
    sjPlot::tab_model(PBfinal,SIRfinal, NTsfinal, NTmfinal, collapse.ci = T, 
                      pred.labels = c("Intercept", "Earthworm (#)", "Grasshopper (#)", "Baseline '16 (of dependent variable)","Earthworm × Grasshopper"),
                      dv.labels = c("Plant Biomass", "Substrate Induced Respiration","Field N mineralization (log)", "Lab N mineralization"), p.style="stars")
  }
}

# No interactions in 2017
YearSeparateModels(pbdata %>% filter(Year == "17"))

# No interactions in 2018
YearSeparateModels(pbdata %>% filter(Year == "18"))

# ... 1.1.7 Plot the results of the univariate linear models ----

# This plot is found in the Appendix S1: Figure S9
pbdata = pbdata %>% mutate(Year = as.factor(Year))


par(mfrow=c(4,2), omi=c(0.5,0,0,0), mar=c(2,4.5,2,1))
# ---- Panel A ----- #
pbdata[,"Base16"] = pbdata$PlantBiomass16
PBfinal = lmer(PlantBiomass ~ WORM_N + HopperN + Year + Base16 + (1|Plot), data= pbdata)

plot(PlantBiomass~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17),
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="Plant Biomass (g)")
ee2 = as.data.frame(Effect("WORM_N",PBfinal))
points(fit~WORM_N, data=ee2, type="l", col="white", lwd=5, lty=1)
points(fit~WORM_N, data=ee2, type="l", col="black", lwd=3, lty=1)

legend("topleft", legend="a", bty="n")
legend("topright", legend=c("2017", "2018"), col=c("blue", "orange"), bty="n", pch=c(19,17), title = "Year")

# ---- Panel B ----- #
plot(PlantBiomass~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17),
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="")

legend("topleft", legend="b", bty="n")
legend("topright", legend = c("Year ***",
                              "Worms *"), bty = "n")
# ---- Panel C ----- #
pbdata[,"Base16"] = pbdata$SIR16
SIRfinal = lmer(SIR ~ WORM_N + HopperN + Base16 + Year + (1|Plot), data= pbdata)

plot(SIR~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab=expression(SIR~(g[C]~g[DMES]^-1~hr^-1)))
legend("topleft", legend="c", bty="n")

# ---- Panel D ----- #
plot(SIR~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="")
legend("topright", legend = c("Year ***"), bty = "n")
legend("topleft", legend="d", bty="n")

# ---- Panel E ----- #
pbdata[,"Base16"] = pbdata$NTs16
NTsfinal = lmer(log(NTs) ~ WORM_N+ HopperN + Base16 + Year + (1|Plot), data= pbdata, na.action = na.omit)

plot(NTs~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab=expression(Field~Mineralization~(mu*g[N]~cm^-2~day^-1)))
ee2 = as.data.frame(Effect("WORM_N",NTsfinal))
ee2$fit = exp(ee2$fit)
points(fit~WORM_N, data=ee2, type="l", col="white", lwd=5, lty=2)
points(fit~WORM_N, data=ee2, type="l", col="black", lwd=3, lty=2)
legend("topleft", legend="e", bty="n")

# ---- Panel F ----- #
plot(NTs~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab="")
ee2 = as.data.frame(Effect("HopperN",NTsfinal))
ee2$fit = exp(ee2$fit)
points(fit~HopperN, data=ee2, type="l", col="white", lwd=5, lty=1)
points(fit~HopperN, data=ee2, type="l", col="black", lwd=3, lty=1)

legend("topleft", legend="f", bty="n")
legend("topright", legend = c("Year ***",
                              "Worm ***",
                              "Hopper **"), bty = "n")

# ---- Panel G ----- #
pbdata[,"Base16"] = pbdata$NTm16
NTmfinal = lmer(NTm ~ WORM_N + HopperN + Base16 + Year + (1|Plot), data= pbdata)

plot(NTm~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab=expression(Lab~Mineralization~(mu*g[N]~g[DMES]^-1~day^-1)))
legend("topleft", legend="g", bty="n")

# ---- Panel H ----- #
plot(NTm~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab="")
legend("topleft", legend="h", bty="n")
legend("topright", legend = c("Year ***"), bty = "n")

mtext(text="Earthworm (residuals)", side=1, outer=T, adj=0.25, line=1)
mtext(text="Grasshopper (#)", side=1, outer=T, adj=0.85, line=1)

pbdata = pbdata %>% select(-Base16)

rm(PBfinal,SIRfinal, NTsfinal, NTmfinal)

#.... 1.1.8 Model the effects of biomass control and year ----

# Subset the data to only include the two controls
bmdata = pbdata %>% filter(HopperAdd!="Add" & Addition!="Add")

bm1 = lmer(PlantBiomass~BiomassCtrl + Year + (1|Plot), data=bmdata); summary(bm1)
bm2 = lmer(SIR~BiomassCtrl+Year + (1|Plot), data=bmdata); summary(bm2)
bm3 = lmer(NTs~BiomassCtrl+ Year + (1|Plot), data=bmdata); summary(bm3)
bm4 = lmer(NTm~BiomassCtrl+ Year + (1|Plot), data=bmdata); summary(bm4)

plot(bm1)
plot(bm2)
plot(bm3)
plot(bm4)

key_label = c(NTm="Lab N Mineralization",
              NTs="Field N Mineralization",
              PlantBiomass= "Plant Biomass",
              SIR="Microbial Biomass")

bmdata2 = bmdata %>% select(BiomassCtrl, PlantBiomass, SIR, NTs, NTm, Year) %>% gather(-BiomassCtrl, -Year, key=Key, value=Measure) %>%
  mutate(Year = ifelse(Year == "17", "2017", "2018"))

ann_text <- data.frame(Year = c(0.75, 0.75, 0.75, 2.1, 2.1), Measure = c(0.6,1.5, 60, 1.5, 1.45),
                       Key = as.factor(c("NTm", "NTs", "PlantBiomass", "SIR", "SIR")),
                       BiomassCtrl= rep("No", 5))

# This is Figure 3
ggplot(bmdata2, aes(x=Year, y=Measure, color=BiomassCtrl)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() + 
  theme_classic() + 
  facet_wrap("Key", scale="free_y", labeller = labeller(Key=key_label)) + 
  geom_text(data = ann_text,label = c("Year*", "Year***", "Year*", "Year***", "Biomass Ctrl*"), color="black") +
  scale_color_discrete(name = "Biomass Control")

# .. 1.2 Plant Community Composition Analysis: RDA -------------------------------

# Scale the data with a hellinger transformation
planthel = decostand(RDAdata %>% select(ACMI:VICR), "hellinger")

# Interaction (RDA with an interaction term)
rda1_1 = rda(planthel~ Earthworm*Grasshopper + Year + PlantBiomass + BioCtrl + Condition(SoilV) + Condition(PlantBiomass16) + Condition(Grass16) + Condition(TRPR16)  + Condition(SoilAdded_Mar16), data=RDAdata)

# No interaction (RDA without an interaction term)
rda1 = rda(planthel~ Earthworm + Grasshopper + Year + PlantBiomass + BioCtrl + Condition(SoilV) + Condition(PlantBiomass16) + Condition(Grass16) + Condition(TRPR16)  + Condition(SoilAdded_Mar16), data=RDAdata)

anova(rda1, rda1_1) # interaction not significant

summary(rda1)
coef(rda1)

(R2 <- RsquareAdj(rda1)$r.squared)

# Global test of the RDA result
anova.cca(rda1, step=1000)
# Tests of all canonical axes
anova.cca(rda1, by="axis", step=1000)
# Tests of all terms
anova.cca(rda1, by="terms", step=1000)
# Tests of all terms
anova.cca(rda1, by="margin", step=1000)

# Get the proportion of variation explained
(prop_explained = round(rda1$CCA$eig/rda1$tot.chi*100,1))


pdf("figure_4.pdf", width=7, height=7)
# Plot the RDA: Figure 4 in the main text
plot(rda1, type="n", xlab=paste0("RDA1 (", prop_explained[1],"%)"), ylab=paste0("RDA2 (", prop_explained[2],"%)"), tck = 0.02)
points(rda1, display= "sites", choices=c(1,2), scaling=2, cex=0.5, col="grey",pch=ifelse(RDAdata$Year=="17", 19, 17))
ordiellipse(rda1, groups=RDAdata$Year,col= alpha(c("black", "grey40"), 0.7), lwd=3, label=F)

text(rda1,labels = c("                              Earthworm", 
                     "     Grasshopper", 
                     "Year",
                     "Plant Biomass", 
                     "                       Biomass Control"), 
     display="bp", choices=c(1,2), scaling=2, col=c("blue","blue","orange","blue"), select=c(T,T, F,T, T), 
     cex = 1)

text(rda1, display= "sp", choices=c(1,2), scaling=2, cex=0.5, col="red")
legend("bottomleft", legend=bquote(italic(R)^2 == .(format(R2, digits = 2))), bty="n")

legend("bottomright", legend=c("2017", "2018"), col=alpha(c("black", "grey40"), 0.7), bty="n", pch=c(19,17))
text(-0.8,1.2, "Significant", col = "orange")
text(-0.8,1.1, "Not significant", col = "blue")
text(-0.3,0.55, "Year", col = "orange")

par(new=TRUE) # overlay existing plot
par(mar=c(0,0,0,0)) # strip out the margins for the inset plot
par(fig=c(0.47,0.94,0.6,0.882)) # fig shrinks and places relative to figure region
plot(0,0, pch=3,type="p", col="grey", xlab="", ylab="", xlim=c(-0.11, 0.09), ylim=c(-0.03, 0.07),xaxt='n',yaxt='n')
abline(h=0, lty = 3); abline(v=0, lty=3)
text(rda1, display= "sp", choices=c(1,2), scaling=2, cex=0.45, col="red")
dev.off()

# .. 1.3 Run Path Analysis (Structural equation models [SEM]) ----

# Modify and standardize the data for the SEM analysis

# Total plant cover calculation
Totals = RDAdata %>% select(ACMI:VICR) %>% rowSums()

# Add in the goldenrod and clover to the data frame standardized for percent of total community
pcdata = pbdata %>% 
  left_join(RDAdata %>% select(Plot, Year, SOAL, TRPR, TRPR16)) %>% 
  mutate(SOAL = SOAL/Totals, TRPR = TRPR/Totals) %>% 
  left_join(SOAL2017 %>% filter(DOE == 612) %>% select(Plot, HtTOT) %>% 
              rename(SOALbase = HtTOT)) %>% mutate(NTs = replace_na(NTs,1.05)) 

rm(Totals)

pbdata_SEM = pcdata %>% 
  
  # Select the right data
  select(c("WORM_N","AP_N", "HopperN", "SIR", "NTs", "NTm", 
           "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
           "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
           "TRPR16", "SOALbase", "PlantBiomass", "Plot", "Year")) %>%
  
  # Rescale the data for an SEM
  mutate_at(c("WORM_N","AP_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
              "TRPR16", "SOALbase"),list(scale = scale)) %>%
  mutate(Year = as.character(Year))

# Reorganize the data for an SEM
pbdata_SEM = pbdata_SEM %>% select(PlantBiomass, Plot, Year, WORM_N_scale:SOALbase_scale) %>%
  gather(-PlantBiomass, -Plot, -Year,key = Variable, value = value) %>%
  separate(Variable, into=c("Variable", "sclae"), sep="_s") %>%
  select(-sclae) %>%
  spread(key = Variable, value = value)

# Define and run the SEMs with (SEM1), with one (SEM2), and without (SEM3) interactions
SEM1 <- psem(lme(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SOAL ~WORM_N*HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
             lme(TRPR ~WORM_N*HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTs ~ WORM_N*HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTm ~ WORM_N*HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SIR ~ WORM_N*HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
             SOAL %~~% PlantBiomass16,
             PlantBiomass %~~% TRPR16,
             PlantBiomass %~~% SIR16,
             SOAL %~~% TRPR16,
             TRPR%~~%SIR,
             TRPR%~~%NTs,
             TRPR%~~%SOAL
)

(SEM1a = summary(SEM1, .progressBar = T))


SEM2 <- psem(lme(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N + HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               PlantBiomass %~~% SIR16,
               SOAL %~~% TRPR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL
)

(SEM2a = summary(SEM2, .progressBar = T))

SEM3 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SOAL ~WORM_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
             lme(TRPR ~WORM_N + HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTs ~ WORM_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTm ~ WORM_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SIR ~ WORM_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
             SOAL %~~% PlantBiomass16,
             PlantBiomass %~~% TRPR16,
             PlantBiomass %~~% SIR16,
             SOAL %~~% TRPR16,
             TRPR%~~%SIR,
             TRPR%~~%NTs,
             TRPR%~~%SOAL
)

(SEM3a = summary(SEM3, .progressBar = T))

# Check the significant clover interaction
SEM4 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SOAL ~WORM_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
             lme(TRPR ~WORM_N*HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTs ~ WORM_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTm ~ WORM_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SIR ~ WORM_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
             SOAL %~~% PlantBiomass16,
             PlantBiomass %~~% TRPR16,
             PlantBiomass %~~% SIR16,
             SOAL %~~% TRPR16,
             TRPR%~~%SIR,
             TRPR%~~%NTs,
             TRPR%~~%SOAL
)

(SEM4a = summary(SEM4, .progressBar = T))

# Try the final SEM with Endogeic earthworms instead of total earthworms
SEM5 <- psem(lme(PlantBiomass ~ AP_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SOAL ~AP_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
             lme(TRPR ~AP_N + HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTs ~ AP_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(NTm ~ AP_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
             lme(SIR ~ AP_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
             SOAL %~~% PlantBiomass16,
             PlantBiomass %~~% TRPR16,
             PlantBiomass %~~% SIR16,
             SOAL %~~% TRPR16,
             TRPR%~~%SIR,
             TRPR%~~%NTs,
             TRPR%~~%SOAL
)

(SEM5a = summary(SEM4, .progressBar = T))

# Generate AIC values for each SEM model
data.frame(Model = c("All H x W interations", "Only plant interaction", "No interaction", "Only clover interaction"),
           AIC = c(AIC(SEM1),AIC(SEM2),AIC(SEM3), AIC(SEM4))
)


# Create an output table for the SEM results
SEM3a$coefficients[,1:8] %>%
  left_join(
    data.frame(Response = unique(SEM3a$coefficients$Response),
               Response2 = c("PB",
                             "Goldenrod",
                             "Clover",
                             "Field N Min.",
                             "Lab N Min.",
                             "SIR",
                             "Goldenrod",
                             "PB",
                             "Clover"))
  ) %>%
  left_join(
    data.frame(Predictor = unique(SEM3a$coefficients$Predictor),
               Predictor2 = c("Earthworm",
                              "Grasshopper",
                              "SIR",
                              "Field N Min.",
                              "Lab N Min.",
                              "PB 2016",
                              "Goldenrod",
                              "Clover",
                              "Year [2018]",
                              "Goldenrod 2016",
                              "Clover 2016",
                              "Field N Min. 2016",
                              "Lab N Min. 2016",
                              "SIR 2016",
                              "PB 2016",
                              "Clover 2016",
                              "SIR 2016",
                              "SIR",
                              "Field N Min.",
                              "Goldenrod")
    )
) %>%
  select(-Response, -Predictor) %>%
  rename(Response = Response2, Predictor = Predictor2) %>%
  select(Response, Predictor, Estimate, Std.Error, DF, P.Value, Std.Estimate) %>%
  mutate(Std.Error = ifelse(is.na(Std.Error), "Correlation", Std.Error)) %>%
  rename('Std Error' = Std.Error,
         p = P.Value,
         'Std Estimate' = Std.Estimate) %>%
    sjPlot::tab_df()

fplot = SEM3a$coefficients[(SEM3a$coefficients$P.Value < 0.05 & !is.na(SEM3a$coefficients$Std.Error)),]

fplot[,"EstimateScale"] = scales::rescale(abs(fplot$Estimate), to = c(1,6))

fplot[,c("Response","Predictor", "EstimateScale")]

# 2.0 Worm Extraction Plot analysis ----

# This is the analysis of the 1-m^2 plots that I surveyed in the field next to the main experiment. This analysis is presented in the supplemental material of this manuscript.

# It follows the same general steps as above for the experimental analysis

# The only treatments are worm removal and addition

# .... 2.0.1 Explore WE worm data ----

wormWE %>%
  filter(SeasonYear == "Fall16") %>%
  summarise(mean(WORM_N, na.rm = T),
            sd(WORM_N, na.rm = T))

ggplot(wormWE,aes(x=SeasonYear, y=WORM_N,color=Treatment)) + geom_boxplot() + theme_classic() + geom_jitter(height = 0) + scale_x_discrete(limits=c("Fall16","Spring17", "Fall17","Spring18", "Fall18"), labels = c("Fall '16","Spring '17", "Fall '17","Spring '18", "Fall '18")) + ylab("Earthworm (#)") + xlab("Season and Year") 
# Fall 2018 is the lowest year, significant relative to other falls, but not the spring numbers.

ggplot(wormWE,aes(x=SeasonYear, y=AP_N,color=Treatment)) + geom_boxplot() + theme_classic() + geom_jitter(height = 0) + scale_x_discrete(limits=c("Fall16","Spring17", "Fall17","Spring18", "Fall18"), labels = c("Fall '16","Spring '17", "Fall '17","Spring '18", "Fall '18")) + ylab("Earthworm (#)") + xlab("Season and Year") 

ggplot(wormWE,aes(x=SeasonYear, y=LUM_N,color=Treatment)) + geom_boxplot() + theme_classic() + geom_jitter(height = 0) + scale_x_discrete(limits=c("Fall16","Spring17", "Fall17","Spring18", "Fall18"), labels = c("Fall '16","Spring '17", "Fall '17","Spring '18", "Fall '18")) + ylab("Earthworm (#)") + xlab("Season and Year") 

# Linear models of earthworm abundance by treatment and abiotic variables (AP = Endogeic; LUM = Epigeic-Anecic; WORM = total)
summary(lmer(WORM_N~Treatment + DOE + SoilTavg + CloudCover + (1|Plot), data= wormWE %>% filter(DOE > 500)))
summary(lm(AP_N~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500)))
summary(lm(LUM_N~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) # not significantly reduced
summary(lm(Worm_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500)))
summary(lm(AP_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) 
summary(lm(LUM_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) # significantly reduced the weight of LUM, but see above not numbers (i.e. colonization)

# .... 2.0.2 Correct worm data for variable extraction efficiency
WE_correctworm0 <- lm(WORM_N~SoilTavg + VWCavg + CloudCover + Air_Temp, data = pbdataWE); summary(WE_correctworm0)
WE_correctworm <- lm(WORM_N~ SoilTavg + VWCavg, data = pbdataWE); summary(WE_correctworm)

pbdataWE = pbdataWE %>% 
  rename(WORM_N_old = WORM_N) %>%
  mutate(WORM_N = resid(WE_correctworm))

# .. 2.1 Linear Models ----
# ... 2.1.1 Model analysis of years together ----

WEmodels = vector("list", 4)

#Plant biomass
WEmodels[[1]] = lmer(PlantBiomass~WORM_N + PlantBiomass16 + Year + (1|Plot), data= pbdataWE); summary(WEmodels[[1]])
plot(WEmodels[[1]])
plot(resid(WEmodels[[1]])~pbdataWE$WORM_N)

# Microbial biomass (SIR)
WEmodels[[2]] = lmer(SIR~WORM_N + SIR16 + Year + (1|Plot), data= pbdataWE); summary(WEmodels[[2]])
plot(WEmodels[[2]])
plot(resid(WEmodels[[2]])~pbdataWE$WORM_N)

# Field nitrogen mineralization
WEmodels[[3]] = lmer(NTs~WORM_N + NTs16 + Year + (1|Plot), data= pbdataWE); summary(WEmodels[[3]])
plot(WEmodels[[3]])
plot(resid(WEmodels[[3]])~pbdataWE$WORM_N)

# Laboratory nitrogen mineralizaiton potential
WEmodels[[4]] = lmer(NTm~WORM_N + NTm16 + Year + (1|Plot), data= pbdataWE); summary(WEmodels[[4]])
plot(WEmodels[[4]])
plot(resid(WEmodels[[4]])~pbdataWE$WORM_N)

# ... 2.1.2 Model analysis of two years separately ----
# Plant biomass
W_PBm17 = lm(PlantBiomass~WORM_N + PlantBiomass16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_PBm18 = lm(PlantBiomass~WORM_N + PlantBiomass16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

# Microbial biomass (SIR)
W_S17 = lm(SIR~WORM_N + SIR16, data= pbdataWE %>% filter(Year=="17")); summary(W_S17)

W_S18 = lm(SIR~WORM_N + SIR16, data= pbdataWE %>% filter(Year=="18")); summary(W_S18)

# Field nitrogen mineralization
W_Ns17 = lm(NTs~WORM_N + NTs16, data= pbdataWE %>% filter(Year=="17")); summary(W_Ns17)

W_Ns18 = lm(NTs~WORM_N + NTs16, data= pbdataWE %>% filter(Year=="18")); summary(W_Ns18)

# Laboratory nitrogen mineralizaiton potential
W_Nm17 = lm(NTm~WORM_N + NTm16, data= pbdataWE %>% filter(Year=="17")); summary(W_Nm17)

W_Nm18 = lm(NTm~WORM_N + NTm16, data= pbdataWE %>% filter(Year=="18")); summary(W_Nm18)

# no significant effects of previous biomass, SIR, NTs, or WORM_N on plant biomass in these plots

# .. 2.2 Plant Community Composition Analysis: RDA ----

planthel = decostand(RDAWdata %>% select(ACMI:VILA), "hellinger")

rda2 = rda(planthel~ WORM_N + Year + PlantBiomass + Treatment + Condition(PlantBiomass16), data=RDAWdata)

summary(rda2)
coef(rda2)

(R2 <- RsquareAdj(rda2)$r.squared)

# Gloabl test of the RDA result
anova.cca(rda2, step=1000)
# Tests of all canonical axes
anova.cca(rda2, by="axis", step=1000)
# Tests of all terms
anova.cca(rda2, by="terms", step=1000)
# Tests of all terms
anova.cca(rda2, by="margin", step=1000)

(prop_explained = round(rda2$CCA$eig/rda2$tot.chi*100,1))

plot(rda2, type="n", xlab=paste0("RDA1 (", prop_explained[1],"%; ns)"), ylab=paste0("RDA2 (", prop_explained[2],"%; ns)"))
points(rda2, display= "sites", choices=c(1,2), scaling=2, pch=ifelse(RDAWdata$Year=="17", 19, 17),col=ifelse(RDAWdata$Treatment=="Remove", "black", "grey"))
text(rda2, display="bp", choices=c(1,2), scaling=2, col=c("blue","orange"), select=c(T,F, T,F))
text(rda2, display= "sp", choices=c(1,2), scaling=2, cex=0.5, col="red")
legend("topright", legend=bquote(italic(R)^2 == .(format(R2, digits = 3))), bty="n")
legend("bottomright", legend=c("`17", "`18", "Control", "Remove"), col=c("grey30", "grey30", "grey", "black"), bty="n", pch=c(19,17, 15,15))