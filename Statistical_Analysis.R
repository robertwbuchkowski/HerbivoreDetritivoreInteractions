# CH4 and WE Data Analysis and Exploration

require(vegan)
require(nlme)
require(lme4)
require(multcomp)
require(car)
require(lmerTest)
require(stargazer)
require(effects)
require(piecewiseSEM)
require(tidyverse)

verbose = F

source("extract_data_model_analysis.R")

dir.create(paste0("Stats_plots_from_",Sys.Date()))

# ..Cover and worm analysis of the data-----
coverworm = percentcover %>% separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>% select(-Date, -ExperimentStart, - DOE, -Season) %>% gather(-Plot, -Year, -DOY, key=Plant, value=Cover) %>% 
  left_join(read_csv("Data/plant_groups_ch4.csv")) %>% # can add this information
  filter(DOY!=150) %>% # remove second spring survey so 2017 and 2018 have the same # of surveys
  filter(!is.na(Cover)) %>% group_by(Plot, Year,Plant, Fcn_grp) %>% 
  summarise(Cover = mean(Cover)) %>% ungroup() %>% group_by(Plot, Year, Fcn_grp) %>% 
  summarise(Cover=sum(Cover)) %>% ungroup() %>% 
  filter(Fcn_grp %in% c("legume", "forb", "grass")) %>% 
  left_join(wormdata2 %>% filter(Season=="Spring") %>% 
              select(Plot, Year, Addition, WORM_N))

coverworm %>% ggplot(aes(x=WORM_N, y=Cover, color=Year)) + geom_point() + theme_classic() + facet_grid(Fcn_grp~.) + stat_smooth(method="lm")

summary(lm(Cover~WORM_N+Year, data=coverworm %>% filter(Fcn_grp =="legume")))
summary(lm(Cover~WORM_N+Year, data=coverworm %>% filter(Fcn_grp =="grass")))
summary(lm(Cover~WORM_N+Year, data=coverworm %>% filter(Fcn_grp =="forb")))

summary(lm(Cover~Addition+Year, data=coverworm %>% filter(Fcn_grp =="legume")))
summary(lm(Cover~Addition+Year, data=coverworm %>% filter(Fcn_grp =="grass")))
summary(lm(Cover~Addition+Year, data=coverworm %>% filter(Fcn_grp =="forb")))

# worm treatment or number don't explain changes in cover...largest effect is time with forbs replacing clover.

rm(coverworm)

# ..Explore WE worm data ----------------------------------------------------
jpeg(paste0("Stats_plots_from_",Sys.Date(),"/Worm_QCQA_WE_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
ggplot(wormWE,aes(x=DOE, y=WORM_N, color=Treatment, group=Plot)) + geom_line() + theme_classic()
dev.off()

ggplot(wormWE,aes(x=DOE, y=AP_N, color=Treatment, group=Plot)) + geom_line() + theme_classic()

ggplot(wormWE,aes(x=DOE, y=LUM_N, color=Treatment, group=Plot)) + geom_line() + theme_classic()

ggplot(wormWE,aes(x=DOE, y=Worm_Wt, color=Treatment, group=Plot)) + geom_line() + theme_classic()  

ggplot(wormWE,aes(x=DOE, y=Worm_Wt/WORM_N, color=Treatment, group=Plot)) + geom_line() + theme_classic()  # no trend in weight of worms over time

summary(lmer(WORM_N~Treatment + DOE + SoilTavg + CloudCover + (1|Plot), data= wormWE %>% filter(DOE > 500)))
summary(lm(AP_N~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500)))
summary(lm(LUM_N~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) # not significantly reduced
summary(lm(Worm_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500)))
summary(lm(AP_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) 
summary(lm(LUM_Wt~Treatment + DOE+ SoilTavg + CloudCover, data= wormWE %>% filter(DOE > 500))) # significantly reduced the weight of LUM, but see above not numbers (i.e. colonization)

we_m1 = lm(WORM_N~Treatment + Season + Year + SoilTavg, data= wormWE); summary(we_m1)
summary(glht(we_m1, linfct=mcp(Year ="Tukey")))
summary(glht(we_m1, linfct=mcp(Season ="Tukey")))

we_m1 = lm(WORM_N~Treatment + SeasonYear + SoilTavg, data= wormWE); summary(we_m1)
summary(glht(we_m1, linfct=mcp(SeasonYear ="Tukey")))

ggplot(wormWE,aes(x=SeasonYear, y=WORM_N)) + geom_boxplot() + theme_classic() + geom_jitter(aes(color=Treatment), height = 0)  # Fall 2018 is the lowest year, significant relative to other falls, but not the spring numbers.

WE %>% mutate(Diff = PlantBiomass17-PlantBiomassSp18) %>% ggplot(aes(x=Treatment, y=Diff)) + geom_boxplot() + geom_jitter(color="red")

t.test(Diff~Treatment, data= WE %>% mutate(Diff = PlantBiomass17-PlantBiomassSp18))

WE %>% mutate(Diff = PlantBiomass17-PlantBiomassSp18) %>% left_join(wormWE %>% filter(SeasonYear == "Spring18") %>% select(Plot, WORM_N)) %>% ggplot(aes(x=WORM_N, y=Diff, color=Treatment)) + geom_point() + theme_classic()

summary(lm(Diff~WORM_N, data= WE %>% mutate(Diff = PlantBiomass17-PlantBiomassSp18) %>% left_join(wormWE %>% filter(SeasonYear == "Spring18") %>% select(Plot, WORM_N))))

# ..Explore CH4 worm data ---------------------------------------------------

wormdata %>% mutate(Addition = ifelse(AP_add_N>0, "Add", "Remove"))%>%ggplot(aes(x=DOE, y=WORM_N, group=Plot, color=Addition)) + geom_line() + theme_classic()

wormdata %>% mutate(Addition = ifelse(AP_add_N>0, "Add", "Remove"))%>%ggplot(aes(x=SeasonYear, y=WORM_N, color=Addition)) + geom_boxplot() + theme_classic()

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/Worm_QCQA_Expt_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
wormdata2 %>% ggplot(aes(x=SeasonYear, y=WORM_N, color=Addition)) + geom_boxplot() + geom_jitter() + theme_classic() + scale_color_discrete(labels = c("Add", "Control", "Remove")) + scale_x_discrete(limits=c("Spring16", "Fall16","Spring17", "Fall17","Spring18", "Fall18")) + ylab("Earthworms (#)") + xlab("Season and Year")
dev.off()

wormdata2 %>% ggplot(aes(x=SeasonYear, y=AP_N, color=Addition)) + geom_boxplot() + geom_jitter() + theme_classic()

wormdata2 %>% ggplot(aes(x=SeasonYear, y=LUM_N, color=Addition)) + geom_boxplot() + geom_jitter() + theme_classic()

m1 = lmer(WORM_N~ Addition + Season + Year + (1|Plot), data=wormdata2 %>% mutate(Year= as.factor(Year))); summary(m1)
summary(glht(m1, linfct = mcp(Addition="Tukey")))
summary(glht(m1, linfct = mcp(Year="Tukey")))
summary(glht(m1, linfct = mcp(Season="Tukey")))
plot(m1)
plot(resid(m1)~Addition, data=wormdata2)
boxplot(resid(m1)~Season, data=wormdata2)
boxplot(resid(m1)~Year, data=wormdata2)

# class(m1) = "lmerMod"

m1 = lmer(WORM_N~ Addition + DOE + (1|Plot), data=wormdata2); summary(m1)
summary(glht(m1, linfct = mcp(Addition="Tukey")))
plot(m1)

m1 = lm(AP_N~ Addition + DOE, data=wormdata2); summary(m1)
summary(glht(m1, linfct = mcp(Addition="Tukey")))
plot(m1) 

m1 = lm(LUM_N~ Addition + DOE, data=wormdata2); summary(m1)
summary(glht(m1, linfct = mcp(Addition="Tukey")))
plot(m1) 
# significantly more worms in the Addition plots than either the Frozen plots or the Removal plots

ch4_wm1 = lmer(WORM_N ~ VWCavg + SoilTavg + CloudCover + Air_Temp + (1|Plot), data= wormdata2); summary(ch4_wm1); plot(ch4_wm1)

ch4_wm1 = lmer(AP_N ~ VWCavg + SoilTavg + CloudCover + Air_Temp + (1|Plot), data= wormdata2); summary(ch4_wm1); plot(ch4_wm1)

ch4_wm1 = lmer(LUM_N ~ VWCavg + SoilTavg + CloudCover + Air_Temp + (1|Plot), data= wormdata2); summary(ch4_wm1); plot(ch4_wm1)

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/Worm_covar_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
wormdata2 %>% select(WORM_N, AP_N, LUM_N, SoilTavg,VWCavg,CloudCover, Air_Temp) %>% gather(-WORM_N, -AP_N, -LUM_N, key=abiotic_key, value=abiotic) %>% gather(-abiotic_key, -abiotic, key=wormtype, value=Number) %>% filter(!is.na(abiotic)) %>% ggplot(aes(x=abiotic, y=Number)) + geom_point() + facet_wrap(abiotic_key~wormtype, scales="free_x") + theme_classic()
dev.off()

ch4_wm1 = mgcv::gamm(WORM_N*10 ~ s(DOY), random=list(Plot=~1), data= wormdata2 %>% filter(Season=="Fall"), family="poisson"); summary(ch4_wm1); anova(ch4_wm1$gam); plot(ch4_wm1$gam); plot(resid(ch4_wm1$lme, type="normalized")~DOY, data= wormdata2 %>% filter(Season=="Fall"))

ch4_wm1 = mgcv::gamm(WORM_N*10 ~ s(DOY), random=list(Plot=~1), data= wormdata2 %>% filter(Season=="Spring"), family="poisson"); summary(ch4_wm1); anova(ch4_wm1$gam); plot(ch4_wm1$gam); plot(ch4_wm1$gam); plot(resid(ch4_wm1$lme, type="normalized")~DOY, data= wormdata2 %>% filter(Season=="Spring"))

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/Worm_doy_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
wormdata2 %>% select(WORM_N, AP_N, LUM_N, DOY, Season)%>% gather(-DOY, - Season, key=wormtype, value=Number) %>% ggplot(aes(x=DOY, y=Number)) + geom_jitter()+ theme_classic() + facet_wrap(Season~wormtype, scale="free_x")
dev.off()

ch4_wm2 = lm(WORM_N ~ VWCavg + SoilTavg + CloudCover + Season + Addition, data= wormdata2); summary(ch4_wm2) # keep only SoilTavg and use Season

1/(1-summary(ch4_wm2)$r.squared) # if VIF > than this, variables are more related to each other than the dependent variable

vif(ch4_wm2) # better, use the above model for worms!

wormdata2 %>% ggplot(aes(x = VWCavg, y= WORM_N, color=Year)) + geom_point() + theme_classic()

wormdata2 %>% ggplot(aes(x = SoilTavg, y= WORM_N, color=Year)) + geom_point() + theme_classic()

wormdata2 %>% ggplot(aes(x = Air_Temp, y= WORM_N, color=Year)) + geom_point() + theme_classic()

wormdata2 %>% ggplot(aes(x = DOY, y= WORM_N, color=Year)) + geom_point() + theme_classic()
# Soil temperature while removing earthworms significantly impacted the number removed

# Test AP and LUM separately

ch4_APm1 = lm(AP_N ~ VWCavg + SoilTavg + CloudCover + Season + Addition, data= wormdata2); summary(ch4_APm1)

1/(1-summary(ch4_APm1)$r.squared) # if VIF > than this, variables are more related to each other than the dependent variable

vif(ch4_APm1) # all less than 5, so OK, but VWCavg is borderline.

ch4_LUMm1 = lm(LUM_N ~ VWCavg + SoilTavg + CloudCover + Season + Addition, data= wormdata2); summary(ch4_LUMm1)

1/(1-summary(ch4_LUMm1)$r.squared) # if VIF > than this, variables are more related to each other than the dependent variable

vif(ch4_LUMm1) # all less than 5, so OK, but VWCavg is borderline.

# SUMMARY: LUM not well manipulated by treatments, but AP was well manipulated. LUM actually higher in removal treatments, which is consistent with their ecological role as a primary colonizer to worm-free areas. Soil temperautre, season, and treatment determine the number of worms removed. Can use SoilTavg as a covariate to models considering the full worm numbers

wormdata3 = wormdata2 %>% filter(SeasonYear == "Fall17"|SeasonYear == "Fall18")

ch4_correctworm <- lm(WORM_N ~ VWCavg + SoilTavg + CloudCover, data= wormdata3); summary(ch4_correctworm) # least significant CloudCover
ch4_correctworm <- lm(WORM_N ~ VWCavg + SoilTavg, data= wormdata3); summary(ch4_correctworm) # both significant
1/(1-summary(ch4_correctworm)$r.squared)
vif(ch4_correctworm)

wormdata3 = wormdata3 %>% 
  select(Plot, Season, Year, WORM_N) %>%
  mutate(WORM_N_STD = as.vector(resid(ch4_correctworm))) %>%
  select(-Season, -WORM_N)

pbdata = pbdata %>% 
  full_join(wormdata3) %>%
  rename(WORM_N_old = WORM_N) %>%
  rename(WORM_N = WORM_N_STD)

# ..Explore VWC data for weather ------------------------------------------

wormdata%>% separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F) %>% filter(Season=="Fall") %>% group_by(Year) %>% summarize(sd(VWCavg), mean(VWCavg), sd(SoilTavg), mean(SoilTavg))

m1 = lm(VWCavg~Year, data = wormdata%>% separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F) %>% filter(Season=="Fall") %>% mutate(Year = as.factor(Year)))
summary(glht(m1, linfct = mcp(Year="Tukey")))

wormWE%>% separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F) %>% filter(Season=="Fall") %>% group_by(Year) %>% summarize(sd(VWCavg), mean(VWCavg), sd(SoilTavg), mean(SoilTavg))

m1 = lm(VWCavg~Year, data = wormWE%>% separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F) %>% filter(Season=="Fall") %>% mutate(Year = as.factor(Year)))
summary(glht(m1, linfct = mcp(Year="Tukey")))

# ..Explore Soil Data -------------------------------------------------------

summary(lm(SIR~FrozenSIR, data=soildata2 %>% filter(SeasonYear=="Spring17"))) # significant negative effect of freezing on SIR data from spring 2017

# Test relationships between nitrogen data and SIR data

soildata2 %>% ggplot(aes(shape=Season)) + geom_point(aes(x=NTm, y = NTs)) + geom_point(aes(x=NO3m, y = NO3s), color="red") + geom_point(aes(x=NH4m, y = NH4s), color="blue") + theme_classic()

summary(lm(NTs~NTm, data= soildata2)) # 10% related
summary(lm(NO3s~NO3m, data= soildata2)) # 12% related
summary(lm(NH4s~NH4m, data= soildata2)) # 20% related

soildata2 %>% ggplot(aes(x=SIR, y = NTm)) + geom_point() + theme_classic()

# Linear models --------------------
# ----------------Run individual linear models---------------- #
#......Effect of biomass control --------------------------------------------

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

bmdata2 = bmdata %>% select(BiomassCtrl, PlantBiomass, SIR, NTs, NTm, Year) %>% gather(-BiomassCtrl, -Year, key=Key, value=Measure)

ann_text <- data.frame(Year = c(0.75, 0.75, 0.75, 2.1, 2.1), Measure = c(0.6,1.5, 60, 1.5, 1.45),
                       Key = as.factor(c("NTm", "NTs", "PlantBiomass", "SIR", "SIR")),
                       BiomassCtrl= rep("No", 5))

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/BiomassCtrl_Year_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
ggplot(bmdata2, aes(x=Year, y=Measure, color=BiomassCtrl)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic() + facet_wrap("Key", scale="free_y", labeller = labeller(Key=key_label)) + geom_text(data = ann_text,label = c("Year*", "Year***", "Year*", "Year***", "BiomassCtrl*"), color="black")
dev.off()
  
#......Plant Biomass Models ----------------------------------------------------

ch4_PBm17 = lm(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16, data= pbdata %>% filter(Year=="17")) ; summary(ch4_PBm17) # interaction not significant
ch4_PBm17_1 = update(ch4_PBm17, .~. - WORM_N:HopperN); summary(ch4_PBm17_1)
1/(1-summary(ch4_PBm17)$r.squared)
vif(ch4_PBm17)
par(mfrow=c(2,2)); plot(ch4_PBm17_1); par(mfrow=c(1,1))

anova(ch4_PBm17, ch4_PBm17_1)

ch4_PBm18 = lm(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16, data= pbdata %>% filter(Year=="18")) ; summary(ch4_PBm18) # interaction not significant
ch4_PBm18_1 = update(ch4_PBm18, .~. - WORM_N:HopperN); summary(ch4_PBm18_1)
1/(1-summary(ch4_PBm18)$r.squared)
vif(ch4_PBm18)
par(mfrow=c(2,2)); plot(ch4_PBm18_1); par(mfrow=c(1,1))

anova(ch4_PBm18, ch4_PBm18_1)

pbdata$Year = as.factor(pbdata$Year)

# Put the Years together
ch4_PBm = lme(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16 + Year, random = ~1|Plot, data= pbdata, method = "ML") ; summary(ch4_PBm)
vif(ch4_PBm)
plot(ch4_PBm)
ch4_PBm_1 = update(ch4_PBm, .~. - WORM_N:HopperN); summary(ch4_PBm_1)

anova(ch4_PBm, ch4_PBm_1)

ch4_PBm = lme(log(PlantBiomass) ~ WORM_N*HopperN + PlantBiomass16 + Year, random = ~1|Plot, data= pbdata) ; summary(ch4_PBm)
vif(ch4_PBm)
plot(ch4_PBm)
ch4_PBm_1 = update(ch4_PBm, .~. - WORM_N:HopperN); summary(ch4_PBm_1)
plot(ch4_PBm_1)
plot(resid(ch4_PBm_1)~pbdata$WORM_N)
boxplot(resid(ch4_PBm_1)~pbdata$HopperN)
boxplot(resid(ch4_PBm_1)~pbdata$Year)
plot(resid(ch4_PBm_1)~pbdata$PlantBiomass16)

ch4_PBm_plot = lme(PlantBiomass ~ WORM_N + HopperN + PlantBiomass16 + Year, random = ~1|Plot, data= pbdata)

ch4_PBm_1 = lme(log(PlantBiomass) ~ WORM_N + HopperN*Year + PlantBiomass16, random = ~1|Plot, data= pbdata) ; summary(ch4_PBm_1)
plot(ch4_PBm_1)
plot(resid(ch4_PBm_1)~pbdata$WORM_N)
boxplot(resid(ch4_PBm_1)~pbdata$HopperN)
boxplot(resid(ch4_PBm_1)~pbdata$Year)
plot(resid(ch4_PBm_1)~pbdata$PlantBiomass16)

#......SIR Models ----------------------------------------------------

ch4_SIRm17_1 = lm(SIR ~ WORM_N*HopperN + SIR16, data= pbdata %>% filter(Year=="17")) ; summary(ch4_SIRm17_1) # interaction not significant
ch4_SIRm17 = update(ch4_SIRm17_1, .~. - WORM_N:HopperN); summary(ch4_SIRm17)
1/(1-summary(ch4_SIRm17)$r.squared)
vif(ch4_SIRm17)

anova(ch4_SIRm17_1,ch4_SIRm17)

ch4_SIRm18_1 = lm(SIR ~ WORM_N*HopperN + SIR16, data= pbdata %>% filter(Year=="18")) ; summary(ch4_SIRm18_1) # interaction not significant
ch4_SIRm18 = update(ch4_SIRm18_1, .~. - WORM_N:HopperN); summary(ch4_SIRm18)
1/(1-summary(ch4_SIRm18)$r.squared)
vif(ch4_SIRm18)

anova(ch4_SIRm18_1,ch4_SIRm18)

# hoppers have the most convincing effect on SIR

# Put the Years together
ch4_SIRm = lme(SIR ~ WORM_N*HopperN + SIR16 + Year, random = ~1|Plot, data= pbdata, method = "ML") ; summary(ch4_SIRm)
vif(ch4_SIRm)
plot(ch4_SIRm)
ch4_SIRm_1 = update(ch4_SIRm, .~. - WORM_N:HopperN); summary(ch4_SIRm_1)

anova(ch4_SIRm, ch4_SIRm_1)

ch4_SIRm = lme(SIR ~ WORM_N*HopperN + SIR16 + Year, random = ~1|Plot, data= pbdata) ; summary(ch4_SIRm)
vif(ch4_SIRm)
plot(ch4_SIRm)
ch4_SIRm_1 = update(ch4_SIRm, .~. - WORM_N:HopperN); summary(ch4_SIRm_1)
plot(resid(ch4_SIRm_1)~pbdata$WORM_N)
boxplot(resid(ch4_SIRm_1)~pbdata$Year)
plot(resid(ch4_SIRm_1)~jitter(pbdata$HopperN))
plot(resid(ch4_SIRm_1)~pbdata$SIR16)

#......NTs Models ----------------------------------------------------

ch4_NTsm17_1 = lm(NTs ~ WORM_N*HopperN + NTs16, data= pbdata %>% filter(Year=="17")) ; summary(ch4_NTsm17_1) # interaction not significant
ch4_NTsm17 = update(ch4_NTsm17_1, .~. - WORM_N:HopperN); summary(ch4_NTsm17)
1/(1-summary(ch4_NTsm17)$r.squared)
vif(ch4_NTsm17)

anova(ch4_NTsm17_1,ch4_NTsm17)

ch4_NTsm18_1 = lm(NTs ~ WORM_N*HopperN + NTs16, data= pbdata %>% filter(Year=="18")) ; summary(ch4_NTsm18_1) # interaction not significant
ch4_NTsm18 = update(ch4_NTsm18_1, .~. - WORM_N:HopperN); summary(ch4_NTsm18)
1/(1-summary(ch4_NTsm18)$r.squared)
vif(ch4_NTsm18)

anova(ch4_NTsm18_1,ch4_NTsm18)

# both increased NTs in 2017 and 2018

# Put the Years together
ch4_NTsm = lme(log(NTs) ~ WORM_N*HopperN + NTs16 + Year, random = ~1|Plot, data= pbdata, method = "ML", na.action = na.omit) ; summary(ch4_NTsm)
vif(ch4_NTsm)
plot(ch4_NTsm)
ch4_NTsm_1 = update(ch4_NTsm, .~. - WORM_N:HopperN); summary(ch4_NTsm_1)

anova(ch4_NTsm, ch4_NTsm_1)

ch4_NTsm = lme(log(NTs) ~ WORM_N*HopperN + NTs16 + Year, random = ~1|Plot, data= pbdata, na.action = na.omit) ; summary(ch4_NTsm)
vif(ch4_NTsm)
plot(ch4_NTsm)
ch4_NTsm_1 = update(ch4_NTsm, .~. - WORM_N:HopperN); summary(ch4_NTsm_1)

ch4_NTsm_plot = lme(NTs ~ WORM_N + HopperN + NTs16 + Year, random = ~1|Plot, data= pbdata, na.action = na.omit)

#......NTm Models ----------------------------------------------------

ch4_NTmm17_1 = lm(NTm ~ WORM_N*HopperN + NTm16, data= pbdata %>% filter(Year=="17")) ; summary(ch4_NTmm17_1) # interaction not significant
ch4_NTmm17 = update(ch4_NTmm17_1, .~. - WORM_N:HopperN); summary(ch4_NTmm17)
1/(1-summary(ch4_NTmm17)$r.squared)
vif(ch4_NTmm17)

anova(ch4_NTmm17_1,ch4_NTmm17)

ch4_NTmm18_1 = lm(NTm ~ WORM_N*HopperN + NTm16, data= pbdata %>% filter(Year=="18")) ; summary(ch4_NTmm18_1) # interaction not significant
ch4_NTmm18 = update(ch4_NTmm18_1, .~. - WORM_N:HopperN); summary(ch4_NTmm18)
1/(1-summary(ch4_NTmm18)$r.squared)
vif(ch4_NTmm18)

anova(ch4_NTmm18_1,ch4_NTmm18)

# worms increased NTm in 2017, no significant effects in 2018

# Put the Years together
ch4_NTmm = lme(NTm ~ WORM_N*HopperN + NTm16 + Year, random = ~1|Plot, data= pbdata, method = "ML") ; summary(ch4_NTmm)
vif(ch4_NTmm)
plot(ch4_NTmm)
ch4_NTmm_1 = update(ch4_NTmm, .~. - WORM_N:HopperN); summary(ch4_NTmm_1)

anova(ch4_NTmm, ch4_NTmm_1)

ch4_NTmm = lme(NTm ~ WORM_N*HopperN + NTm16 + Year, random = ~1|Plot, data= pbdata) ; summary(ch4_NTmm)
vif(ch4_NTmm)
plot(ch4_NTmm)
ch4_NTmm_1 = update(ch4_NTmm, .~. - WORM_N:HopperN); summary(ch4_NTmm_1)



#......Output Models -----------------------------------------------------

stargazer(ch4_PBm17_1,ch4_PBm18_1,ch4_SIRm17,ch4_SIRm18,ch4_NTsm17,ch4_NTsm18,ch4_NTmm17,ch4_NTmm18,star.cutoffs = c(0.05, 0.01, 0.001), column.labels = rep(c("2017", "2018"),4))

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/Figure3_",Sys.Date(), ".jpeg"), units="in", width=6, height=8, res=600)

par(mfrow=c(4,2), omi=c(0.5,0,0,0), mar=c(2,4.5,2,1))
# ---- Panel A ----- #
summary(ch4_PBm_plot)
plot(PlantBiomass~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17),
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="Plant Biomass (g)")
ee2 = as.data.frame(Effect("WORM_N",ch4_PBm_plot))
points(fit~WORM_N, data=ee2, type="l", col="white", lwd=5, lty=1)
points(fit~WORM_N, data=ee2, type="l", col="black", lwd=3, lty=1)

legend("topleft", legend="a", bty="n")
legend("topright", legend=c("2017", "2018"), col=c("blue", "orange"), bty="n", pch=c(19,17), title = "Year")

# ---- Panel B ----- #
plot(PlantBiomass~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17),
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="")

legend("topleft", legend="b", bty="n")
legend("topright", legend = c("Year ***",
                              "Worms *",
                              "Starting biomass ***"), bty = "n")
# ---- Panel C ----- #
summary(ch4_SIRm_1)
plot(SIR~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab=expression(SIR~(g[C]~g[DMES]^-1~hr^-1)))
legend("topleft", legend="c", bty="n")

# ---- Panel D ----- #
plot(SIR~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col=ifelse(Year=="17", "blue","orange"), xlab="", ylab="")
legend("topright", legend = c("Year ***"), bty = "n")
legend("topleft", legend="d", bty="n")

# ---- Panel E ----- #
summary(ch4_NTsm_plot)
plot(NTs~WORM_N, data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab=expression(Field~Mineralization~(mu*g[N]~cm^-2~day^-1)))
ee2 = as.data.frame(Effect("WORM_N",ch4_NTsm_plot))
points(fit~WORM_N, data=ee2, type="l", col="white", lwd=5, lty=2)
points(fit~WORM_N, data=ee2, type="l", col="black", lwd=3, lty=2)
legend("topleft", legend="e", bty="n")

# ---- Panel F ----- #
plot(NTs~jitter(HopperN), data=pbdata, pch=ifelse(Year=="17", 19,17), 
     col= ifelse(Year=="17", "blue","orange"), xlab="", ylab="")
ee2 = as.data.frame(Effect("HopperN",ch4_NTsm_plot))
points(fit~HopperN, data=ee2, type="l", col="white", lwd=5, lty=1)
points(fit~HopperN, data=ee2, type="l", col="black", lwd=3, lty=1)

legend("topleft", legend="f", bty="n")
legend("topright", legend = c("Year ***",
                              "Worm +",
                              "Hopper *"), bty = "n")

# ---- Panel G ----- #
summary(ch4_NTmm_1)
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
dev.off()

# Test models on a subset of the data of high quality ---------------------

# Overall general effect size stays the same, but getting rid of the bunch of worm data makes the PB interaction no longer significant.

pbdata %>% select(Year, WORM_N) %>% group_by(Year) %>% summarize(lq = quantile(WORM_N, c(0.25)),med = quantile(WORM_N, c(0.5)), uq = quantile(WORM_N, c(0.75)))

pbdata %>% select(Year, WORM_N) %>% group_by(Year) %>% summarize(lq = quantile(WORM_N, c(0.33)), uq = quantile(WORM_N, c(0.66)))

pbHQ = pbdata %>% filter(Decimated ==0) %>% 
  filter(HopperAdd=="Add" & HopperN >0 | 
           HopperAdd!="Add" & HopperN ==0) %>% 
  filter(Addition=="Add" & Year =="17" & WORM_N >8 | 
           Addition!="Add" & Year =="17" & WORM_N <13 |
           Addition=="Add" & Year =="18" & WORM_N >0 | 
           Addition!="Add" & Year =="18" & WORM_N <2 )

# Full separation
pbHQ = pbdata %>% filter(Decimated ==0) %>% 
  filter(HopperAdd=="Add" & HopperN >0 | 
           HopperAdd!="Add" & HopperN ==0) %>% 
  filter(Addition=="Add" & Year =="17" & WORM_N >9 | 
           Addition!="Add" & Year =="17" & WORM_N <9 |
           Addition=="Add" & Year =="18" & WORM_N >1.5 | 
           Addition!="Add" & Year =="18" & WORM_N <1.5 )


#......Plant Biomass Models ---------------------------------------------------

ch4_PBm17 = lm(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16, data= pbHQ %>% filter(Year=="17")) ; summary(ch4_PBm17) # interaction not significant
ch4_PBm17_1 = update(ch4_PBm17, .~. - WORM_N:HopperN); summary(ch4_PBm17_1)
1/(1-summary(ch4_PBm17)$r.squared)
vif(ch4_PBm17)
par(mfrow=c(2,2)); plot(ch4_PBm17_1); par(mfrow=c(1,1))

anova(ch4_PBm17, ch4_PBm17_1)

ch4_PBm18 = lm(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16, data= pbHQ %>% filter(Year=="18")) ; summary(ch4_PBm18) # interaction not significant
ch4_PBm18_1 = update(ch4_PBm18, .~. - WORM_N:HopperN); summary(ch4_PBm18_1)
1/(1-summary(ch4_PBm18)$r.squared)
vif(ch4_PBm18)
par(mfrow=c(2,2)); plot(ch4_PBm18_1); par(mfrow=c(1,1))

anova(ch4_PBm18, ch4_PBm18_1)

pbHQ$Year = as.factor(pbHQ$Year)

# Put the Years together
ch4_PBm = lme(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16 + Year, random = ~1|Plot, data= pbHQ, method = "ML") ; summary(ch4_PBm)
vif(ch4_PBm)
plot(ch4_PBm)
ch4_PBm_1 = update(ch4_PBm, .~. - WORM_N:HopperN); summary(ch4_PBm_1)

anova(ch4_PBm, ch4_PBm_1)

ch4_PBm = lme(PlantBiomass ~ WORM_N*HopperN + PlantBiomass16 + Year, random = ~1|Plot, data= pbHQ) ; summary(ch4_PBm)
vif(ch4_PBm)
plot(ch4_PBm)
ch4_PBm_1 = update(ch4_PBm, .~. - WORM_N:HopperN); summary(ch4_PBm_1)

#......SIR Models ----------------------------------------------------

ch4_SIRm17_1 = lm(SIR ~ WORM_N*HopperN + SIR16, data= pbHQ %>% filter(Year=="17")) ; summary(ch4_SIRm17_1) # interaction not significant
ch4_SIRm17 = update(ch4_SIRm17_1, .~. - WORM_N:HopperN); summary(ch4_SIRm17)
1/(1-summary(ch4_SIRm17)$r.squared)
vif(ch4_SIRm17)

anova(ch4_SIRm17_1,ch4_SIRm17)

ch4_SIRm18_1 = lm(SIR ~ WORM_N*HopperN + SIR16, data= pbHQ %>% filter(Year=="18")) ; summary(ch4_SIRm18_1) # interaction not significant
ch4_SIRm18 = update(ch4_SIRm18_1, .~. - WORM_N:HopperN); summary(ch4_SIRm18)
1/(1-summary(ch4_SIRm18)$r.squared)
vif(ch4_SIRm18)

anova(ch4_SIRm18_1,ch4_SIRm18)

# hoppers have the most convincing effect on SIR

# Put the Years together
ch4_SIRm = lme(SIR ~ WORM_N*HopperN + SIR16 + Year, random = ~1|Plot, data= pbHQ, method = "ML") ; summary(ch4_SIRm)
vif(ch4_SIRm)
plot(ch4_SIRm)
ch4_SIRm_1 = update(ch4_SIRm, .~. - WORM_N:HopperN); summary(ch4_SIRm_1)

anova(ch4_SIRm, ch4_SIRm_1)

ch4_SIRm = lme(SIR ~ WORM_N*HopperN + SIR16 + Year, random = ~1|Plot, data= pbHQ) ; summary(ch4_SIRm)
vif(ch4_SIRm)
plot(ch4_SIRm)
ch4_SIRm_1 = update(ch4_SIRm, .~. - WORM_N:HopperN); summary(ch4_SIRm_1)

#......NTs Models ----------------------------------------------------

ch4_NTsm17_1 = lm(NTs ~ WORM_N*HopperN + NTs16, data= pbHQ %>% filter(Year=="17")) ; summary(ch4_NTsm17_1) # interaction not significant
ch4_NTsm17 = update(ch4_NTsm17_1, .~. - WORM_N:HopperN); summary(ch4_NTsm17)
1/(1-summary(ch4_NTsm17)$r.squared)
vif(ch4_NTsm17)

anova(ch4_NTsm17_1,ch4_NTsm17)

ch4_NTsm18_1 = lm(NTs ~ WORM_N*HopperN + NTs16, data= pbHQ %>% filter(Year=="18")) ; summary(ch4_NTsm18_1) # interaction not significant
ch4_NTsm18 = update(ch4_NTsm18_1, .~. - WORM_N:HopperN); summary(ch4_NTsm18)
1/(1-summary(ch4_NTsm18)$r.squared)
vif(ch4_NTsm18)

anova(ch4_NTsm18_1,ch4_NTsm18)

# both increased NTs in 2017 and 2018

# Put the Years together
ch4_NTsm = lme(log(NTs) ~ WORM_N*HopperN + NTs16 + Year, random = ~1|Plot, data= pbHQ, method = "ML", na.action = na.omit) ; summary(ch4_NTsm)
vif(ch4_NTsm)
plot(ch4_NTsm)
ch4_NTsm_1 = update(ch4_NTsm, .~. - WORM_N:HopperN); summary(ch4_NTsm_1)

anova(ch4_NTsm, ch4_NTsm_1)

ch4_NTsm = lme(log(NTs) ~ WORM_N*HopperN + NTs16 + Year, random = ~1|Plot, data= pbHQ, na.action = na.omit) ; summary(ch4_NTsm)
vif(ch4_NTsm)
plot(ch4_NTsm)
ch4_NTsm_1 = update(ch4_NTsm, .~. - WORM_N:HopperN); summary(ch4_NTsm_1)

ch4_NTsm_plot = lme(NTs ~ WORM_N + HopperN + NTs16 + Year, random = ~1|Plot, data= pbHQ, na.action = na.omit)

#......NTm Models ----------------------------------------------------

ch4_NTmm17_1 = lm(NTm ~ WORM_N*HopperN + NTm16, data= pbHQ %>% filter(Year=="17")) ; summary(ch4_NTmm17_1) # interaction not significant
ch4_NTmm17 = update(ch4_NTmm17_1, .~. - WORM_N:HopperN); summary(ch4_NTmm17)
1/(1-summary(ch4_NTmm17)$r.squared)
vif(ch4_NTmm17)

anova(ch4_NTmm17_1,ch4_NTmm17)

ch4_NTmm18_1 = lm(NTm ~ WORM_N*HopperN + NTm16, data= pbHQ %>% filter(Year=="18")) ; summary(ch4_NTmm18_1) # interaction not significant
ch4_NTmm18 = update(ch4_NTmm18_1, .~. - WORM_N:HopperN); summary(ch4_NTmm18)
1/(1-summary(ch4_NTmm18)$r.squared)
vif(ch4_NTmm18)

anova(ch4_NTmm18_1,ch4_NTmm18)

# worms increased NTm in 2017, no significant effects in 2018

# Put the Years together
ch4_NTmm = lme(NTm ~ WORM_N*HopperN + NTm16 + Year, random = ~1|Plot, data= pbHQ, method = "ML") ; summary(ch4_NTmm)
vif(ch4_NTmm)
plot(ch4_NTmm)
ch4_NTmm_1 = update(ch4_NTmm, .~. - WORM_N:HopperN); summary(ch4_NTmm_1)

anova(ch4_NTmm, ch4_NTmm_1)

ch4_NTmm = lme(NTm ~ WORM_N*HopperN + NTm16 + Year, random = ~1|Plot, data= pbHQ) ; summary(ch4_NTmm)
vif(ch4_NTmm)
plot(ch4_NTmm)
ch4_NTmm_1 = update(ch4_NTmm, .~. - WORM_N:HopperN); summary(ch4_NTmm_1)

# Explore WE data ----

W_PBm17 = lm(PlantBiomass~WORM_N + PlantBiomass16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_PBm18 = lm(PlantBiomass~WORM_N + PlantBiomass16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_S17 = lm(SIR~WORM_N + SIR16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_S18 = lm(SIR~WORM_N + SIR16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_Ns17 = lm(NTs~WORM_N + NTs16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_Ns18 = lm(NTs~WORM_N + NTs16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_Nm17 = lm(NTm~WORM_N + NTm16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm18)

W_Nm18 = lm(NTm~WORM_N + NTm16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

stargazer(W_PBm17,W_PBm18,W_S17,W_S18,W_Ns17,W_Ns18,W_Nm17,W_Nm18)

# no significant effects of previous biomass, SIR, NTs, or WORM_N on plant biomass in these plots

W_PBm17 = lm(PlantBiomass~Treatment + PlantBiomass16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_PBm18 = lm(PlantBiomass~Treatment + PlantBiomass16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_S17 = lm(SIR~Treatment + SIR16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_S18 = lm(SIR~Treatment + SIR16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_Ns17 = lm(NTs~Treatment + NTs16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm17)

W_Ns18 = lm(NTs~Treatment + NTs16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

W_Nm17 = lm(NTm~Treatment + NTm16, data= pbdataWE %>% filter(Year=="17")); summary(W_PBm18)

W_Nm18 = lm(NTm~Treatment + NTm16, data= pbdataWE %>% filter(Year=="18")); summary(W_PBm18)

stargazer(W_PBm17,W_PBm18,W_S17,W_S18,W_Ns17,W_Ns18,W_Nm17,W_Nm18)

# Run SEM Models ----
# ...Run SEM Model to test relationships between variables -----

pbdata_SEM = pbdata %>% mutate(NTs = replace_na(NTs,1.05)) %>% 
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N"),list(scale)) %>% 
  filter(Year=="17")


model1 <- psem(lm(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16, data=pbdata_SEM),
               lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
               lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
               lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM)
               )

(model1a = summary(model1, .progressBar = T))

pbdata_SEM = pbdata %>% mutate(NTs = replace_na(NTs,1.05)) %>% 
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm",
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N"),funs(scale)) %>% 
  filter(Year=="18")


model1 <- psem(lm(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16, data=pbdata_SEM),
               lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
               lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
               lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM)
)

(model1b = summary(model1, .progressBar = T))

model1ab = rbind(model1a$coefficients,model1b$coefficients)

colnames(model1ab)[9] = "Sig."

model1ab["Year"] = rep(c(2017, 2018), each=15)

model1ab = model1ab[,c("Year",colnames(model1ab)[1:9])]

stargazer(model1ab, summary = F)

# ...Run SEM Model including common plants (Presented in the paper) -----

Totals = RDAdata %>% select(ACMI:VICR) %>% rowSums()

pcdata = pbdata %>% 
  left_join(RDAdata %>% select(Plot, Year, SOAL, TRPR, TRPR16)) %>% 
  mutate(SOAL = SOAL/Totals, TRPR = TRPR/Totals) %>% 
  left_join(SOAL2017 %>% filter(DOE == 612) %>% select(Plot, HtTOT) %>% 
              rename(SOALbase = HtTOT)) %>% mutate(NTs = replace_na(NTs,1.05)) 

rm(Totals)

pbdata_SEM = pcdata %>%
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
              "TRPR16", "SOALbase"),funs(scale)) %>% filter(Year=="17") %>% 
  select(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
           "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
           "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
           "TRPR16", "SOALbase", "PlantBiomass", "Plot")) 

model1 <- psem(lm(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR, data=pbdata_SEM),
               lm(SOAL ~WORM_N + HopperN + SOALbase, data=pbdata_SEM),
               lm(TRPR ~WORM_N + HopperN + TRPR16, data=pbdata_SEM),
               lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
               lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
               lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM),
               SIR %~~% PlantBiomass16,
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               SOAL %~~% TRPR16,
               SIR %~~% SOALbase,
               NTs %~~% SOALbase,
               NTm %~~% NTs16,
               TRPR%~~%SIR16,
               SOAL%~~%NTm,
               TRPR%~~%SOAL
)

(model1a = summary(model1, .progressBar = T))

pbdata_SEM = pcdata %>%
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
              "TRPR16", "SOALbase"),funs(scale))  %>% filter(Year=="18") %>% 
  select(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
           "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
           "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
           "TRPR16", "SOALbase", "PlantBiomass", "Plot")) 

model1 <- psem(lm(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR, data=pbdata_SEM),
               lm(SOAL ~WORM_N + HopperN + SOALbase, data=pbdata_SEM),
               lm(TRPR ~WORM_N + HopperN + TRPR16, data=pbdata_SEM),
               lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
               lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
               lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM),
               SIR %~~% PlantBiomass16,
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               SOAL %~~% TRPR16,
               SIR %~~% SOALbase,
               NTs %~~% SOALbase,
               NTm %~~% NTs16,
               TRPR%~~%SIR16,
               SOAL%~~%NTm,
               TRPR%~~%SOAL
)

(model1b = summary(model1, .progressBar = T))

model1ab = rbind(model1a$coefficients,model1b$coefficients)

colnames(model1ab)[9] = "Sig."

model1ab["Year"] = rep(c(2017, 2018), each=dim(model1ab)[1]/2)

model1ab = model1ab[,c("Year",colnames(model1ab)[1:9])]

model1ab$Response = gsub("~","",model1ab$Response)
model1ab$Predictor = gsub("~","",model1ab$Predictor)

model1ab = model1ab %>% mutate(Std.Error = replace_na(Std.Error, "Correlation"))

stargazer(model1ab, summary = F, rownames = F)


par(mfrow=c(1,2))
plot(PlantBiomass~SOAL, data=pcdata, pch=ifelse(Year=="17", 19,17))
plot(PlantBiomass~TRPR, data=pcdata, pch=ifelse(Year=="17", 19,17))

# Standardize parameter estimates for visual plot 
# ... plot built in powerpoint
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

model1ab %>% filter(Sig. !="", Std.Error != "Correlation") %>% 
  mutate(Std.Estimate = 10*range01(Std.Estimate)+1)

rm(range01)

# ...Run SEM Model including common plants and years together -----

Totals = RDAdata %>% select(ACMI:VICR) %>% rowSums()

pcdata = pbdata %>% 
  left_join(RDAdata %>% select(Plot, Year, SOAL, TRPR, TRPR16)) %>% 
  mutate(SOAL = SOAL/Totals, TRPR = TRPR/Totals) %>% 
  left_join(SOAL2017 %>% filter(DOE == 612) %>% select(Plot, HtTOT) %>% 
              rename(SOALbase = HtTOT)) %>% mutate(NTs = replace_na(NTs,1.05)) 

rm(Totals)

pbdata_SEM = pcdata %>% 
  select(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
           "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
           "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
           "TRPR16", "SOALbase", "PlantBiomass", "Plot", "Year")) %>%
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
              "TRPR16", "SOALbase"),funs(scale))

model1 <- psem(lme(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N + HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               PlantBiomass %~~% SIR16,
               SOAL %~~% TRPR16,
               # TRPR%~~%SIR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL
               # TRPR%~~%NTm,
               # NTs%~~%SIR,
               # NTm%~~%NTs
)

(model1a = summary(model1, .progressBar = T))

model1 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N + HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N + HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N + HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N + HopperN + NTm16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N + HopperN + SIR16 + Year,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               PlantBiomass %~~% SIR16,
               SOAL %~~% TRPR16,
               # TRPR%~~%SIR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL
               # TRPR%~~%NTm,
               # NTs%~~%SIR,
               # NTm%~~%NTs
)

(model2a = summary(model1, .progressBar = T))





model1 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR + Year,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N+HopperN + SOALbase + Year,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N+HopperN + TRPR16 + Year,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N+HopperN + NTs16 + Year,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               SOAL %~~% TRPR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL,
               TRPR%~~%NTm,
               NTs%~~%SIR,
               NTm%~~%NTs
)

(model2b = summary(model1, .progressBar = T))


model1 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N+HopperN + NTs16,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N+HopperN + NTm16,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N+HopperN + SIR16,random=~1|Plot, data=pbdata_SEM),
               NTs%~~%SIR,
               NTm%~~%NTs
)

(model3a = summary(model1, .progressBar = T))

# .... Run with pbHQ data -----

pbdata_SEM = pcdata %>% 
  right_join(pbHQ %>% select(Plot, Year)) %>% # select HQ data
  select(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
           "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
           "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
           "TRPR16", "SOALbase", "PlantBiomass", "Plot")) %>%
  mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
              "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
              "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", "TRPR", 
              "TRPR16", "SOALbase"),funs(scale))

model1 <- psem(lme(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N + HopperN + SOALbase,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N + HopperN + TRPR16,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N + HopperN + NTs16,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N + HopperN + NTm16,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N + HopperN + SIR16,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               SIR%~~%SOALbase,
               PlantBiomass%~~%SIR16,
               SOAL %~~% TRPR16,
               TRPR%~~%SIR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL
)

(model1a = summary(model1, .progressBar = T))

model1 <- psem(lme(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR,random=~1|Plot, data=pbdata_SEM),
               lme(SOAL ~WORM_N*HopperN + SOALbase,random=~1|Plot, data=pbdata_SEM),
               lme(TRPR ~WORM_N*HopperN + TRPR16,random=~1|Plot, data=pbdata_SEM),
               lme(NTs ~ WORM_N*HopperN + NTs16,random=~1|Plot, data=pbdata_SEM),
               lme(NTm ~ WORM_N*HopperN + NTm16,random=~1|Plot, data=pbdata_SEM),
               lme(SIR ~ WORM_N*HopperN + SIR16,random=~1|Plot, data=pbdata_SEM),
               SOAL %~~% PlantBiomass16,
               PlantBiomass %~~% TRPR16,
               SIR%~~%SOALbase,
               PlantBiomass%~~%SIR16,
               SOAL %~~% TRPR16,
               TRPR%~~%SIR16,
               TRPR%~~%SIR,
               TRPR%~~%NTs,
               TRPR%~~%SOAL
)

(model2a = summary(model1, .progressBar = T))

# ....Run AIC model comparison of interactions --------------------------------

YeartoRun = c("17", "18")

for(i in 1:2){
  pbdata_SEM = pcdata %>% 
    mutate_at(c("WORM_N", "HopperN", "SIR", "NTs", "NTm", 
                "PlantBiomass16", "SoilTavg","VWCavg", "NTm16", 
                "NTs16", "SIR16", "AP_N", "LUM_N", "SOAL", 
                "TRPR", "TRPR16", "SOALbase"),funs(scale)) %>% 
    filter(Year==YeartoRun[i])
  
  model1 <- psem(lm(PlantBiomass ~ WORM_N + HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR, data=pbdata_SEM),
                 lm(SOAL ~WORM_N + HopperN + SOALbase, data=pbdata_SEM),
                 lm(TRPR ~WORM_N + HopperN + TRPR16, data=pbdata_SEM),
                 lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
                 lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
                 lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM),
                 SIR %~~% PlantBiomass16,
                 SOAL %~~% PlantBiomass16,
                 PlantBiomass %~~% TRPR16,
                 SOAL %~~% TRPR16,
                 SIR %~~% SOALbase,
                 NTs %~~% SOALbase,
                 NTm %~~% NTs16
  )
  
  model2 <- psem(lm(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR, data=pbdata_SEM),
                 lm(SOAL ~WORM_N + HopperN + SOALbase, data=pbdata_SEM),
                 lm(TRPR ~WORM_N + HopperN + TRPR16, data=pbdata_SEM),
                 lm(NTs ~ WORM_N + HopperN + NTs16, data=pbdata_SEM),
                 lm(NTm ~ WORM_N + HopperN + NTm16, data=pbdata_SEM),
                 lm(SIR ~ WORM_N + HopperN + SIR16, data=pbdata_SEM),
                 SIR %~~% PlantBiomass16,
                 SOAL %~~% PlantBiomass16,
                 PlantBiomass %~~% TRPR16,
                 SOAL %~~% TRPR16,
                 SIR %~~% SOALbase,
                 NTs %~~% SOALbase,
                 NTm %~~% NTs16
  )
  
  
  model3 <- psem(lm(PlantBiomass ~ WORM_N*HopperN + SIR + NTs + NTm + PlantBiomass16 + SOAL + TRPR, data=pbdata_SEM),
                 lm(SOAL ~WORM_N*HopperN + SOALbase, data=pbdata_SEM),
                 lm(TRPR ~WORM_N*HopperN + TRPR16, data=pbdata_SEM),
                 lm(NTs ~ WORM_N*HopperN + NTs16, data=pbdata_SEM),
                 lm(NTm ~ WORM_N*HopperN + NTm16, data=pbdata_SEM),
                 lm(SIR ~ WORM_N*HopperN + SIR16, data=pbdata_SEM),
                 SIR %~~% PlantBiomass16,
                 SOAL %~~% PlantBiomass16,
                 PlantBiomass %~~% TRPR16,
                 SOAL %~~% TRPR16,
                 SIR %~~% SOALbase,
                 NTs %~~% SOALbase,
                 NTm %~~% NTs16
  )
  print(YeartoRun[i])
  print(cbind(
    Model = c("No interaction", "Plant interaction","All interactions"),
    AIC = c(AIC(model1),AIC(model2),AIC(model3))
  ))

}

# No interactions are selected based on AIC

# Plant Community Composition Analysis: RDA -------------------------------
# ...Expt Plots ---------------------------------------------------------

planthel = decostand(RDAdata %>% select(ACMI:VICR), "hellinger")

# Interaction
rda1_1 = rda(planthel~ Earthworm*Grasshopper + Year + PlantBiomass + BioCtrl + Condition(SoilV) + Condition(PlantBiomass16) + Condition(Grass16) + Condition(TRPR16)  + Condition(SoilAdded_Mar16), data=RDAdata)

# No interaction
rda1 = rda(planthel~ Earthworm + Grasshopper + Year + PlantBiomass + BioCtrl + Condition(SoilV) + Condition(PlantBiomass16) + Condition(Grass16) + Condition(TRPR16)  + Condition(SoilAdded_Mar16), data=RDAdata)

# Decimated is significant...

anova(rda1, rda1_1) # interaction not significant

summary(rda1)
coef(rda1)

(R2 <- RsquareAdj(rda1)$r.squared)

# Gloabl test of the RDA result
anova.cca(rda1, step=1000)
# Tests of all canonical axes
anova.cca(rda1, by="axis", step=1000)
# Tests of all terms
anova.cca(rda1, by="terms", step=1000)
# Tests of all terms
anova.cca(rda1, by="margin", step=1000)

(prop_explained = round(rda1$CCA$eig/rda1$tot.chi*100,1))

dev.off()
jpeg(paste0("Stats_plots_from_",Sys.Date(),"/RDA_Expt_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
plot(rda1, type="n", xlab=paste0("RDA1 (", prop_explained[1],"%)"), ylab=paste0("RDA2 (", prop_explained[2],"%)"))
points(rda1, display= "sites", choices=c(1,2), scaling=2, cex=0.5, col="grey",pch=ifelse(RDAdata$Year=="17", 19, 17))
# ordiellipse(rda1, groups=RDAdata$wxh,col= c("black", "grey40"), lwd=3, label=T)
ordiellipse(rda1, groups=RDAdata$Year,col= c("black", "grey40"), lwd=3, label=F)
text(rda1, display="bp", choices=c(1,2), scaling=2, col=c("blue","blue","orange","blue"), select=c(T,T, F,T, T))
#ordiellipse(rda1, groups=RDAdata$BiomassCtrl,col= c("black", "grey40"), lwd=3, label=T)
text(rda1, display= "sp", choices=c(1,2), scaling=2, cex=0.5, col="red")
legend("bottomleft", legend=bquote(italic(R)^2 == .(format(R2, digits = 3))), bty="n")

legend("bottomright", legend=c("`17", "`18"), col=c("black", "grey40"), bty="n", pch=c(19,17))

par(new=TRUE) # overlay existing plot
par(mar=c(0,0,0,0)) # strip out the margins for the inset plot
par(fig=c(0.47,0.94,0.6,0.882)) # fig shrinks and places relative to figure region
plot(0,0, pch=3,type="p", col="grey", xlab="", ylab="", xlim=c(-0.11, 0.09), ylim=c(-0.03, 0.07),xaxt='n',yaxt='n')
abline(h=0, lty = 3); abline(v=0, lty=3)
text(rda1, display= "sp", choices=c(1,2), scaling=2, cex=0.45, col="red")
dev.off()

# ...Worm Extraction Plots ---------------------------------------------------

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

jpeg(paste0("Stats_plots_from_",Sys.Date(),"/RDA_WE_",Sys.Date(), ".jpeg"), units="in", width=7, height=7, res=600)
plot(rda2, type="n", xlab=paste0("RDA1 (", prop_explained[1],"%; ns)"), ylab=paste0("RDA2 (", prop_explained[2],"%; ns)"))
points(rda2, display= "sites", choices=c(1,2), scaling=2, pch=ifelse(RDAWdata$Year=="17", 19, 17),col=ifelse(RDAWdata$Treatment=="Remove", "black", "grey"))
text(rda2, display="bp", choices=c(1,2), scaling=2, col=c("blue","orange"), select=c(T,F, T,F))
#ordiellipse(rda2, groups=RDAWdata$Year,col= c("black", "grey"), lwd=3, label=T)
#ordiellipse(rda2, groups=RDAWdata$Treatment,col= c("purple", "cyan"), lwd=3, lty=2, label=T)
text(rda2, display= "sp", choices=c(1,2), scaling=2, cex=0.5, col="red")
legend("topright", legend=bquote(italic(R)^2 == .(format(R2, digits = 3))), bty="n")
legend("bottomright", legend=c("`17", "`18", "Control", "Remove"), col=c("grey30", "grey30", "grey", "black"), bty="n", pch=c(19,17, 15,15))
dev.off()
