# Extract data for model analysis

# meant to be run in-line with in the script "Statistical_Analysis.R"

verbose = F # set to true if running dynamically to see analysis at the bottom
createmodeldata = F # set true if you want to re-create the model data from the raw data frames

# ---- Set up sampling times ----
hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309)+365, c(323, 114)+365*2)
wormsampWE = c(311, c(114, 314)+365, c(99, 280)+365*2)

wormavgsample = data.frame(SeasonYear = c("Spring16","Fall16","Spring17","Fall17","Spring18","Fall18"),
                           time = wormsamp[order(wormsamp)])

soilavgsample = data.frame(SeasonYear = c("Fall15","Spring16","Fall16","Spring17","Fall17","Spring18","Fall18"),
                           time = soilsamp[order(soilsamp)])

wormWEavgsample = data.frame(SeasonYear = c("Fall16","Spring17","Fall17","Spring18","Fall18"),
                             time = wormsampWE[order(wormsampWE)])

# ---- Load the appropriate data ----

# Load WE data frame 
WE <- read_csv("Data/plotcharacteristics_WE_Dec2018.csv", col_types = cols(Plot=col_character()))

# Load worm WE data frame 
wormWE <- read_csv("Data/Worm_data_WE_Dec2018.csv", 
                   col_types = cols(Date = col_date(format = "%m/%d/%y"),
                                    Plot =col_character(),
                                    ExperimentStart = col_date(format = "%m/%d/%y"), 
                                    TimeEnd = col_time(format = "%H:%M:%S"), 
                                    TimeStart = col_time(format = "%H:%M:%S")))  %>%
  left_join(WE %>% select(Plot,Treatment)) %>% 
  separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F) %>% 
  mutate(SeasonYear = as.factor(SeasonYear), Season = as.factor(Season), Year=as.factor(Year))

wormWE$SeasonYear = as.factor(wormWE$SeasonYear)

# Load percentcoverWE data
percentcoverWE <- read_csv("Data/plant_community_summary_WE.csv", 
                           col_types = cols(Date = col_date(format = "%m/%d/%y"),
                                            Plot = col_character(),
                                            PLMA = col_integer()))%>% select(Plot:ASDU) %>% 
  replace(is.na(.), 0) %>% 
  rename(FRsp = FRVE, OXsp = 'Oxalis sp.', ASsp = ASSE) %>% mutate(VILA = FoxGrape + UK17) %>% select(-FoxGrape, -UK17) %>% # confirmed from photo
  mutate(RUHI = RUHI + UK15) %>% select(-UK15) %>% #confirmed from photo
  mutate(PLsp = PLMA + UK4 + UK8) %>% select(-UK4, - UK8) %>% #confirmed from photo
  mutate(ASLA = ASLA + UK16) %>% select(-UK16) %>% #confirmed from photo
  mutate(UKSeedling = UKSeedling + UK18) %>% select(-UK18) %>% #joined, ID uncertain
  mutate(UK = FuzzyRound + UKPointLeaf + HOCA) %>% select(-FuzzyRound, -UKPointLeaf, - HOCA) %>% separate(Date, into=c("Year", NA), sep=4, remove=F) %>% separate(Year, into=c(NA, "Year"), sep=2)

# Load experimental plot data
plotdata <- read_csv("Data/plotcharacteristicdata.csv", col_types = cols(Plot=col_character(),'Soil Group'=col_character(),SoilAdded_Mar16=col_character())) %>% 
  left_join(read_csv("Data/Hopper_14Dec2018.csv", col_types = cols(Plot = col_character()))) %>% mutate(Hopperf_N_Su18 = ifelse(is.na(Number_Hopper), 0, Number_Hopper)) %>% select(-Number_Hopper) %>% rename(SoilGroup = 'Soil Group',Submerged17=Submerged_2April17, Submerged18 = Submerged_7Nov2018) %>% mutate(Submerged17 = ifelse(Submerged17==100,1,0)) %>% left_join(read_csv("Data/hopperwt_2018.csv", col_types = cols(Plot=col_character()))) %>% mutate(HopperN18 = replace_na(N,0), HopperWt18=replace_na(Wt,0)) %>% select(-N, -Wt)

# Load experimental soil data
soildata <- read_csv("Data/soil_data_CH4_Jan2019.csv", col_types = cols(Plot = col_character())) %>% gather(-Project, -Plot, key= Type, value= Measure) %>% separate(Type, into=c("Type", "SeasonYear"), sep="_") %>% separate(SeasonYear, into=c("Season", "Year"), sep=1) %>% mutate(Season = ifelse(Season=="S", "Spring", "Fall")) %>% unite(SeasonYear, Season, Year, sep="", remove=F)

soildata2 = soildata %>% 
  filter(Project =="CH4") %>% 
  select(-Project) %>% 
  spread(key=Type, value=Measure) %>% 
  mutate(FrozenSIR = ifelse(is.na(FrozenSIR), 0, FrozenSIR))

# Load experimental earthworm data
wormdata <- read_csv("Data/Worm_data_master_Dec2018.csv", 
                     col_types = cols(AP_N = col_double(), 
                                      AP_Wt = col_double(), Date = col_date(format = "%m/%d/%y"), 
                                      ExperimentStart = col_date(format = "%m/%d/%y"), 
                                      LUM_N = col_double(), LUM_Wt = col_double(), 
                                      Plot = col_character(), SoilT1 = col_double(), 
                                      SoilT2 = col_double(), SoilT3 = col_double(), 
                                      SoilTavg = col_double(), VWC1 = col_double(), 
                                      VWC2 = col_double(), VWC3 = col_double(),
                                      TimeEnd = col_time(format = "%H:%M:%S"), 
                                      TimeStart = col_time(format = "%H:%M:%S")))

wormdata2 = wormdata %>% left_join(wormdata %>% filter(SeasonYear == "Fall17") %>% select(Plot, AP_add_N,Frozen_add) %>% mutate(Addition = ifelse(Frozen_add==1, "Frozen",ifelse(AP_add_N>0, "Add", "Remove"))) %>% select(Plot, Addition)
) %>% mutate(Addition = ifelse(SeasonYear == "Spring16","Remove",Addition)) %>% mutate(Addition = as.factor(Addition)) %>% separate(SeasonYear, into=c("Season", "Year"), sep=-2, remove=F)

# Load experimental animal data
aboveanimals <- read_csv("Data/CH4_notwormanimals_Dec18.csv", 
                         col_types = cols(Date = col_date(format = "%m/%d/%y"), 
                                          ExpermentStart = col_date(format = "%m/%d/%y"), 
                                          Herb_SOAL_Percent = col_number(), 
                                          Hopper3rd = col_number(), Hopper4th = col_number(), 
                                          Hopper5th = col_number(), HopperWt = col_number(), 
                                          Plot = col_character()))

# Load experimental percent cover data
percentcover <- read_csv("Data/CH4_percentcover_Sept10.csv", 
                         col_types = cols(Date = col_date(format = "%m/%d/%y"), 
                                          ExperimentStart = col_date(format = "%m/%d/%y"), 
                                          MANE = col_integer(), Plot = col_character())) %>% mutate(ASsp = ASsp + UK6 + PointedLeaf) %>% select(-UK6, - PointedLeaf) %>% # add UK6 to ASVA because photos show it is the same
  mutate(GAAP = GAAP + UK7) %>% select(-UK7) %>% # add UK7 to GAAP because field notes identified it as GAAP
  select(-UKSEEDLING) %>% #looks like an oak, but not using the 2016 data were this point is
  mutate(PEPE = PEPE + UK1847) %>% select(-UK1847) %>% # lab ID from photos and later visits
  mutate(UK = UK + UK821 + UK1452 + Ukclusterdroop) %>% select(-UK821, -UK1452, -Ukclusterdroop) # a group of small, immature plants eatten/dead before ID could be obtained

# Load experimental Solidago growth data
SOAL2017 <- read_csv("Data/CH4_SOAL2017.csv", 
                     col_types = cols(Date = col_date(format = "%m/%d/%y"), 
                                      ExperimentStart = col_date(format = "%m/%d/%y"), 
                                      Plot = col_character()))

# .........Identify plots with significant caterpillar damage----

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# now created a set of robust vectors based on 2 criteria:
# 1. has a caterpillar
# 2. has greater than 25% herbivory on Solidago
# 3. because Solidago grows back well the percent cover estimate is probably not great

# 2017 => 10,18,20,47,48
# 2018 => 2, 10, 14,18,29,48,49
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

addtopb = aboveanimals %>% 
  separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>% 
  group_by(Plot, Year) %>% 
  summarize(TotWB = sum(WoolyBear, na.rm=T),
            TotC = sum(Caterpillar, na.rm=T), 
            SpittleBug = sum(SpittleBug, na.rm=T), 
            Moth = sum(Moth, na.rm=T),
            Kaytidid= sum(Kaytidid, na.rm=T),
            HERB = max(Herb_SOAL_Percent, na.rm=T)) %>% 
  mutate(TotCat = TotWB + TotC) %>% select(Plot, Year, TotCat, HERB)

list17 = as.character(c(10,18,20,47,48))
list18 = as.character(c(2, 10, 14,18,29,48,49))

plotdata = plotdata %>% mutate(Decimated17 =ifelse(Plot %in% list17, 1,0), Decimated18 = ifelse(Plot %in% list18, 1,0))

rm(list17, list18)

# ---- Create experimental dataframe to join together ----

pbdata = plotdata %>% select(Plot, PlantBiomass16,PlantBiomass17,PlantBiomass18,HopperN17,HopperWt17,HopperN18,HopperWt18,SoilV, Submerged17, Submerged18, SoilGroup, SoilAdded_Mar16, Decimated17, Decimated18) %>% gather(-Plot, -PlantBiomass16, - SoilV, -SoilGroup,-SoilAdded_Mar16, key=Key, value=Measure) %>% separate(Key, into=c("Type", "Year"), sep=-2) %>% spread(key=Type, value=Measure) %>% 
  mutate(Submerged = as.character(Submerged)) %>%# creates the plot data, with hoppers and plant biomass, plant biomass 2016 is the baseline
  left_join(wormdata2 %>% filter(Season=="Fall" & Year!="16") %>% select(Plot, Year, WORM_N, Worm_Wt, AP_N, LUM_N, AP_Wt, LUM_Wt, SoilTavg,VWCavg, LUTE,CloudCover, Addition)) %>% # Adding the worm data for 2017 and 2018
  left_join(aboveanimals %>% filter(HopperAdd >0) %>% select(Plot) %>% add_column(HopperAdd = "Add")) %>% mutate(HopperAdd = ifelse(is.na(HopperAdd), "Remove",ifelse(Addition=="Frozen", "Frozen", "Add"))) %>% # add hoppper treatments
  left_join(soildata2 %>% filter(Season=="Fall" & Year!="16") %>% select(-SeasonYear, -Season, -FrozenSIR, -pH)) %>% # Adding soil data for 2017 and 2018
  left_join(soildata2 %>% filter(SeasonYear=="Spring16") %>% select(-SeasonYear, -Season,-Year, -FrozenSIR, -pH) %>% gather(-Plot, key=Key, value=Measure) %>% mutate(Key2 = rep("16",600)) %>% unite(Key, Key, Key2, sep="") %>% spread(key=Key, value=Measure)) %>% # Adding baseline plot data
  full_join(addtopb) %>% # add TotCat and HERB
  mutate(Decimated = as.factor(Decimated)) %>% # make Decimated a factor
  mutate(BiomassCtrl = ifelse(Addition=="Frozen", "Yes", "No")) %>%
  distinct() #remove duplicated rows

pbdataWE = WE %>% select(Plot, Treatment, PlantBiomass16,PlantBiomass17,PlantBiomass18) %>% gather(-Plot,-Treatment, -PlantBiomass16, key=Key, value=Measure) %>% separate(Key, into=c("Type", "Year"), sep=-2) %>% spread(key=Type, value=Measure) %>% # creates the plot data, with hoppers and plant biomass, plant biomass 2016 is the baseline
  left_join(wormWE %>% separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>% filter(Season=="Fall" & Year!="16")%>% select(Plot, Year, WORM_N, Worm_Wt, AP_N, LUM_N, AP_Wt, LUM_Wt, SoilTavg)) %>% # Adding the worm data for 2017 and 2018
  left_join(soildata %>% filter(Project =="W") %>% select(-Project) %>% spread(key=Type, value=Measure) %>% mutate(FrozenSIR = ifelse(is.na(FrozenSIR), 0, FrozenSIR)) %>% filter(Season=="Fall" & Year!="16") %>% select(-SeasonYear, -Season, -FrozenSIR, -pH)) %>% # Adding soil data for 2017 and 2018
  left_join(soildata %>% filter(Project =="W") %>% select(-Project) %>% spread(key=Type, value=Measure) %>% mutate(FrozenSIR = ifelse(is.na(FrozenSIR), 0, FrozenSIR)) %>% filter(Year=="16") %>% select(-SeasonYear, -Season, -FrozenSIR, -pH) %>% select(Plot, Year, SIR, NTs, NTm) %>% gather(-Plot, -Year, key=Key, value=Value) %>% mutate(Key = paste0(Key, Year)) %>% select(-Year) %>% spread(key=Key, value=Value)) #Add 2016 soil data

# ---- Create data for RDA ----

# Experiment
RDAdata_1 = percentcover %>% separate(SeasonYear, into=c("Season", "Year"), sep=-2) %>% 
  select(-Date, -ExperimentStart, - DOE, -Season) %>%
  filter(DOY > 190 & Year !="16") %>% select(-DOY) %>%# remove spring surveys to make 2017 and 2018 consistent, and remove 2016 baseline data (will be predictors)
  gather(-Plot, -Year, key=Species, value=Cover) %>% group_by(Plot, Year, Species) %>% summarise(Cover = mean(Cover)) %>% spread(key=Species, value=Cover) %>% # create full data set of plant cover
  select(Plot:VICR) %>% gather(-Plot, -Year, key=Species, value=Cover)

RDAdata = RDAdata_1 %>% #right_join(RDAdata_1 %>% group_by(Species) %>% summarize(Cover = mean(Cover)) %>% filter(Cover >= 0.5) %>% ungroup() %>% select(Species)) %>% # subset that data set to only include plants with high cover
  spread(key=Species, value=Cover) %>%
  full_join(pbdata %>% select(Plot:SIR)) %>% ungroup() %>%
  mutate(wxh = paste0(ifelse(HopperAdd=="Add", "H",ifelse(HopperAdd=="Frozen","fH","0")),ifelse(Addition=="Add", "W", ifelse(Addition=="Frozen","fW","0")))) %>%
  mutate(BioCtrl = ifelse(Addition=="Frozen", "N", "Y")) %>%
  left_join(percentcover %>% filter(SeasonYear == "Fall16") %>%
              select(Plot, Grass, TRPR) %>% rename(Grass16 = Grass, TRPR16 = TRPR)) %>%  # add in baseline plant data
  rename(Earthworm = WORM_N, Grasshopper = HopperN)
rm(RDAdata_1)

# WE plots
RDAWdata = percentcoverWE %>% select(-Date, - Date2, -Treatment) %>% gather(-Plot, -Year, key=Species, value=Cover) %>% group_by(Plot, Year, Species) %>% summarize(Cover= mean(Cover)) %>% ungroup() %>% spread(key=Species, value=Cover) %>%
  right_join(pbdataWE)





# ---- Convert data to the model formats -----

if(createmodeldata){
  # Identify the correct treatments in model syntax
  CH4treatments = pbdata %>% 
    mutate(Treatment = ifelse(Addition=="Add", 
                              ifelse(HopperAdd=="Add", "HW","W"), 
                              ifelse(HopperAdd=="Add", "H","N"))) %>% 
    select(Plot, Treatment) %>% distinct()
  
  datatomodel = bind_rows(
    
    # Start with the earthworm data
    wormdata %>% select(Plot, SeasonYear, Worm_Wt) %>% 
      left_join(wormavgsample) %>% #add in model standardized dates
      mutate(W = Worm_Wt*0.1/0.212264) %>% # convert worm weight to g[N]~m^-2; worms 10% nitrogen
      select(-Worm_Wt, -SeasonYear) %>%
      left_join(CH4treatments) %>% # add the treatments
      gather(-Plot, -time, -Treatment, key=StateVar, value = Biomass),
    
    # Now add the grasshopper and plant data
    plotdata %>% 
      select(Plot, PlantBiomass16,PlantBiomass17,PlantBiomass18,HopperWt17,HopperWt18) %>% 
      gather(-Plot, key=Key, value=Measure) %>% 
      separate(Key, into=c("Type", "Year"), sep=-2) %>% 
      spread(key=Type, value=Measure) %>% 
      mutate(HopperWt = replace_na(HopperWt, 0)) %>% 
      left_join(
        data.frame(
          Year = c("16", "17", "18"),
          time = c(665, 1030, 1395)
        )
      ) %>%
      mutate(P = PlantBiomass*0.02/0.849056, # plants 2% nitrogen, standardized area to 1-m^2
             H = HopperWt*0.11/0.849056) %>% # grasshoppers 11% nitrogen, standardized area to 1-m^2
      select(-Year, -HopperWt, -PlantBiomass) %>%
      left_join(CH4treatments) %>% # add the treatments
      gather(-Plot, -time, -Treatment, key=StateVar, value = Biomass),
    
    # WORM EXTRACTION PLOTS  
    wormWE %>% select(Plot, Treatment, SeasonYear, Worm_Wt) %>% 
      filter(!is.na(Worm_Wt)) %>% # get rid of plots not sampled
      left_join(wormWEavgsample) %>% #add in model standardized dates
      mutate(Biomass = Worm_Wt*0.1/0.212264) %>% # convert worm weight to g[N]~m^-2
      select(-Worm_Wt, -SeasonYear) %>% 
      mutate(Treatment = ifelse(Treatment =="Remove", "RmW", "Rt")) %>%
      mutate(StateVar = "W") %>%
      select(Plot, time, Treatment, StateVar, Biomass),
    
    
    WE %>% 
      select(Plot, Treatment, PlantBiomass16,PlantBiomass17,PlantBiomass18) %>% 
      gather(-Plot,-Treatment, key=Key, value=Measure) %>% 
      separate(Key, into=c("Type", "Year"), sep=-2) %>% 
      spread(key=Type, value=Measure) %>% 
      left_join(
        data.frame(
          Year = c("16", "17", "18"),
          time = c(665, 1030, 1395)
        )
      ) %>%
      mutate(Biomass = PlantBiomass*0.02) %>% # plants 2% N, already 1-m^2
      select(-Year, -PlantBiomass) %>%
      mutate(Treatment = ifelse(Treatment =="Remove", "RmW", "Rt")) %>% # add treats
      mutate(StateVar = "P") %>%
      select(Plot, time, Treatment, StateVar, Biomass),
    
    # Soil data for both projects
    soildata %>% 
      select(-Season, -Year) %>% 
      filter(Type %in% c("NTe", "SIR")) %>%
      left_join(soilavgsample) %>% select(-SeasonYear) %>% #add in model standardized dates
      spread(key=Type, value=Measure)  %>% 
      mutate(N = NTe*49000*1e-6,
             M = ((((((((SIR/(12.01/44.01))/(1000*1000*44.01))*295*0.08205746)/0.998621))*100000)*40.04 + 0.37)/100)*49000*1000*1e-6
      ) %>% 
      select(Project, Plot, time, M, N) %>% 
      gather(-Project, -Plot, -time, key=StateVar, value=Biomass) %>% left_join(
        bind_rows(
          CH4treatments %>% mutate(Project = "CH4"),
          WE %>% select(Plot, Treatment) %>% 
            mutate(Treatment = ifelse(Treatment =="Remove", "RmW", "Rt")) %>% 
            mutate(Project ="W"))
      ) %>% # add in the treatment codes
      select(-Project) %>%
      select(Plot, time, Treatment, StateVar, Biomass)
    
  ) %>%
    
    # Modify plot names by project
    left_join(
      
      bind_rows(
        CH4treatments %>% mutate(Project = "E"),
        WE %>% select(Plot, Treatment) %>% 
          mutate(Treatment = ifelse(Treatment =="Remove", "RmW", "Rt")) %>% 
          mutate(Project ="W")) %>%
        mutate(Plot2 = paste0(Project, Plot)) %>%
        select(-Project)
    ) %>%
    select(-Plot) %>%
    rename(Plot = Plot2) %>%
    select(Plot, time, Treatment, StateVar, Biomass)
  
  if(verbose) datatomodel %>% write_csv("Data/datatomodel.csv")
  
  # ---- Load and analyze other people's data
  
  #Experimental Data
  
  datatomodel = datatomodel %>% 
    mutate(Experiment = "Rob_WE")
  
  # Add grasshopper data from experiments
  zack2 = read_csv("Data/zackmiller_data.csv") %>% 
    mutate(SOL = Sol_Biomass*Sol_perN_post/100,
           GR = Gr_Biomass*Gr_perN_post/100) %>% 
    mutate(PTOT = SOL + GR, PTOTb = Sol_Biomass + Gr_Biomass) %>% 
    mutate(SOLb = Sol_Biomass, GRb = Gr_Biomass) %>% 
    filter(Env =="C", Trophic <3) %>% select(ID, Trophic, SOL:GRb) %>% 
    mutate(Field="Xmas", Experiment="Zack", Trophic=as.factor(Trophic))
  
  adam2 = read_csv("Data/adamhouston_data.csv") %>% 
    filter(N.Treatment=="None" & Trophic.Treatment < 3 & !is.na(totmass)) %>% 
    mutate(Trophic.Treatment = as.factor(Trophic.Treatment)) %>% 
    mutate(SOL = smass*solidago.final.percentN/100,
           GR = gmass*grass.final.percentN/100) %>% 
    mutate(PTOT = SOL + GR, PTOTb = smass + gmass) %>% 
    mutate(Trophic = Trophic.Treatment, SOLb = smass, 
           GRb = gmass,
           HopperN = Grasshopper.weight*Grasshopper.percentN/100) %>% 
    mutate(ID = paste0("Adam",Cage.ID)) %>%
    select(ID, Field, SOL:HopperN) %>% filter(!is.na(PTOT)) %>% 
    mutate(Experiment = "Adam")
  
  hopperdata = bind_rows(zack2, adam2) %>%
    select(ID, Trophic, Experiment, PTOT, HopperN) %>%
    mutate(P = PTOT/ 3.14*0.26*0.26, # correct area
           H = HopperN/ 3.14*0.26*0.26) %>% # correct area
    select(-PTOT, -HopperN) %>%
    gather(-ID, -Trophic, - Experiment, key=StateVar, value = Biomass) %>%
    filter(!(is.na(Biomass))) %>% 
    mutate(Treatment = ifelse(Trophic==2, "RtH","RmH"), time=300) %>% 
    rename(Plot = ID) %>%
    select(Plot, time, Treatment, StateVar, Biomass, Experiment)
  
  # Add MESc data
  
  MESc_data <- read_csv("Data/MESc_data.csv") %>% 
    filter(Treatment %in% c("A", "B")) %>% # select only hopper treatments
    
    # Create plot ID
    left_join(
      data.frame(Site = rep(c("AO", "CP", "Ps", "Rr", "WR"), 2),
                 Treatment = rep(c("A","B"), each = 5),
                 Plot = as.character(paste0("MESc", seq(1,10))))
    ) %>%
    select(-Site) %>%
    mutate(Treatment = ifelse(Treatment =="A", "RmH", "RtH"), # rename treatments
           area = 3.14*(0.25/2)^2,
           time = 300) %>%
    mutate(H = MEFE*0.11/area, # convert hopper to g-N m^-2
           P = SOLb*0.02/area, #convert plants to gN m^-2
           N = TN_ug_gdw*1e-6*(0.1764*3.14*12.5*12.5*30)/area, # Bulk density g cm^-3 times area cm^3 ; convert inorganic N to gN m^-2
           M = ((((((((SIR/(12.01/44.01))/(1000*1000*44.01))*295*0.08205746)/0.998621))*100000)*40.04 + 0.37)/100)*49000*1000*1e-6
    ) %>% 
    select(Plot, Treatment,time, P, N, M) %>% 
    gather(-Plot, -Treatment, -time, key=StateVar, value=Biomass) %>%
    mutate(Experiment = "Rob_MESc") %>%
    select(Plot, time, Treatment, StateVar, Biomass, Experiment)
  
  Os06 = data.frame(PBio = c((118.609+118.2+118.2)/3/3.14*0.26*0.26, 80.164/3.14*0.26*0.26),
                    Treatment = c("RmH", "RtH"),
                    Experiment = c("Os06", "Os06"),
                    time = c(1395, 1395)) %>% 
    mutate(Biomass = PBio*0.02, StateVar="P") %>% 
    select(-PBio) %>% bind_rows(
      data.frame(Biomass = c(0, 0.0644),
                 Treatment = c("RmH", "RtH"),
                 Experiment = c("Os06", "Os06"),
                 time = c(1395, 1395)) %>% 
        mutate(StateVar="H")
    ) %>%
    mutate(Plot = c("Os06_1","Os06_2","Os06_1","Os06_2")) %>%
    select(Plot, time, Treatment, StateVar, Biomass, Experiment)
  
  # Final data
  
  datatomodel2 = bind_rows(datatomodel, hopperdata,MESc_data, Os06)
  
  rm(datatomodel, hopperdata,MESc_data, Os06, zack2, adam2)
  
  write_csv(datatomodel2, "Data/datatomodel2.csv")
  
  saveRDS(datatomodel2, "Data/datatomodel2.Rds")
  
}

# ---- Load climate data and make function ----

if(createmodeldata){
  
  climate = read_csv("Data/temp_records.csv") %>% gather(key=Month, value = Value,-Type, -Year, -Day) %>% spread(key=Type, value=Value) %>% filter(Temp != "-") %>% mutate(Precip = ifelse(Precip=="T",0, Precip)) %>% mutate(Temp = as.numeric(ifelse(Temp=="M", -9999,Temp)),Precip = as.numeric(ifelse(Precip=="M", -9999,Precip))) %>% mutate(Temp = (Temp-32)*5/9, Precip= Precip*2.54)
  
  climate = climate %>% mutate(Month2 = match(Month,month.abb))
  climate = climate %>% mutate(Date = as.Date(paste(Year, Month2,Day,sep="-"),format = "%Y-%m-%d")) 
  climate = climate%>% mutate(DOY =  as.numeric(strftime(climate$Date, format = "%j")))
  
  climate2 = climate %>% filter(Year==2015) %>% select(DOY, Temp) %>% mutate(DOY = replace_na(DOY, 366))
  
  climate2[climate2$DOY==366,2] = climate2[climate2$DOY==365,2]
  
  climate2["TempK"] = climate2$Temp + 273.15
  
  # Model with cosine wave
  m1 = nls(TempK ~ a*cos(2*3.14/365*DOY + c) + b, start=c(a=-12, b=281.15, c=0), data=climate2)
  
  summary(m1)
  
  # curve(coef(m1)[1]*cos(2*3.14/365*x + coef(m1)[3])+coef(m1)[2], 0, 365, col="green", add=T)
  
  LTtemp = function(doy){
    
    -12.8244*cos(2*3.14/365*doy-0.3666)+281.9846
    
  }
  
  pdf("temp_approx.pdf")
  plot(TempK~DOY, data=climate2[order(climate2$DOY),], type="l")
  points(seq(1,365), LTtemp(seq(1,365)), type="l", col="purple")
  dev.off()
  
}

# ---- Additional analysis of Experimental Data for model ----

if(verbose){
  
  # Summarize growth of Solidago in plots...
  SOAL2017 %>% filter(DOY %in% c(171, 192, 267)) %>%
    mutate(time = ifelse(DOY==171, "One", 
                         ifelse(DOY==192, "Two", "Three"))) %>%
    mutate(PN = (-6.789731 + 0.148969*HtTOT)*0.02) %>% # from Cerina
    select(Plot, time, PN) %>%
    spread(key=time, value=PN) %>%
    mutate(growtha = (Three-One)/(267-171),
           growthb = (Three-Two)/(267-192),
           growthc = (Two-One)/(192-171)
    ) %>%
    select(Plot,growtha, growthb, growthc) %>%
    summarize(mean(growthc))
  
  # Make it a percent change, so no transformation is required
  SOAL2017 %>% filter(DOY %in% c(171, 192, 267)) %>%
    mutate(time = ifelse(DOY==171, "One", 
                         ifelse(DOY==192, "Two", "Three"))) %>%
    mutate(PN = (-6.789731 + 0.148969*HtTOT)*0.02) %>% # from Cerina
    select(Plot, time, PN) %>%
    spread(key=time, value=PN) %>%
    mutate(growth = (Three-One)/(One)) %>%
    select(Plot,growth) %>%
    summarize(sd(growth), mean(growth))
  
  
  SOAL2017 %>% filter(DOY %in% c(171, 192, 267)) %>%
    mutate(time = ifelse(DOY==171, "One", 
                         ifelse(DOY==192, "Two", "Three"))) %>%
    mutate(PN = (-6.789731 + 0.148969*HtTOT)*0.02) %>% # from Cerina
    select(Plot, DOY, PN) %>%
    left_join(
      pbdata %>% 
        mutate(Treatment = ifelse(Addition=="Add", 
                                  ifelse(HopperAdd=="Add", "HW","W"), 
                                  ifelse(HopperAdd=="Add", "H","N"))) %>% 
        select(Plot, Treatment) %>% distinct()
    ) %>%
    ggplot(aes(x=DOY, y=PN, group = Plot)) + geom_line() +
    theme_classic() + facet_wrap(.~Treatment)
  
  
  SOAL2017 %>% filter(DOY %in% c(171, 192, 267)) %>%
    mutate(time = ifelse(DOY==171, "One", 
                         ifelse(DOY==192, "Two", "Three"))) %>%
    mutate(PN = (-6.789731 + 0.148969*HtTOT)*0.02) %>% # from Cerina
    select(Plot, time, PN) %>%
    left_join(
      pbdata %>% 
        mutate(Treatment = ifelse(Addition=="Add", 
                                  ifelse(HopperAdd=="Add", "HW","W"), 
                                  ifelse(HopperAdd=="Add", "H","N"))) %>% 
        select(Plot, Treatment) %>% distinct()
    ) %>%
    select(Plot,time, PN) %>%
    spread(key = time, value = PN) %>%
    select(One, Two,Three) %>%
    cor()
  
  # Species Turnover Rate ---------------------------------------------------
  
  percentcoverWE %>% select(-Date, - Date2, -Treatment) %>% gather(-Plot, -Year, key=Species, value=Cover) %>% group_by(Plot, Year, Species) %>% summarize(Cover= mean(Cover)) %>% ungroup() %>% spread(key=Year, value=Cover) %>% mutate(P67_dis = ifelse(`17`==0 & `16`!=0,1,0),P67_ap = ifelse(`16`==0 & `17`!=0,1,0),P78_dis = ifelse(`18`==0 & `17`!=0,1,0),P78_ap = ifelse(`17`==0 & `18`!=0,1,0)) %>% select(-`16`, -`17`, -`18`, -Species) %>% gather(-Plot, key=Type, value=Part) %>% separate(Type, into=c("Time", "Type"),sep="_") %>% group_by(Plot, Time, Type) %>% summarize(Count = sum(Part)) %>% ungroup() %>% group_by(Type, Plot) %>% arrange(Time) %>% mutate(cs = cumsum(Count)) %>% group_by(Type) %>% summarise(sd=sd(Count), Count=mean(Count))
  
  percentcoverWE %>% select(-Date, - Date2, -Treatment) %>% gather(-Plot, -Year, key=Species, value=Cover) %>% group_by(Plot, Year, Species) %>% summarize(Cover= mean(Cover)) %>% ungroup() %>% spread(key=Year, value=Cover)  
  
  # Figure out earthworm population reduction -------------------------------
  
  wormdata %>% filter(SeasonYear=="Spring16") %>% summarize(Nbar = mean(WORM_N),Wbar = mean(Worm_Wt))
  
  wormWE %>% filter(SeasonYear=="Spring17") %>% summarize(Nbar = mean(WORM_N, na.rm=T),Wbar = mean(Worm_Wt, na.rm=T))
  
  wormWE %>% filter(SeasonYear=="Fall16") %>% summarize(Nbar = mean(WORM_N, na.rm=T),Wbar = mean(Worm_Wt, na.rm=T))
  
  # Test for lag effect -----------------------------------------------------
  
  m1 = lm(PlantBiomass ~ WORM_N + HopperN , data=pbdata %>% filter(Year == "17") %>% select(Plot, WORM_N, HopperN, AP_N, LUM_N) %>% full_join(pbdata %>% filter(Year == "18") %>% select(Plot, PlantBiomass)))
  summary(m1)
  
  m1 = lm(NTs ~ WORM_N + HopperN , data=pbdata %>% filter(Year == "17") %>% select(Plot, WORM_N, HopperN, AP_N, LUM_N) %>% full_join(pbdata %>% filter(Year == "18") %>% select(Plot, NTs)))
  summary(m1)
  
  m1 = lm(SIR ~ WORM_N + HopperN , data=pbdata %>% filter(Year == "17") %>% select(Plot, WORM_N, HopperN, AP_N, LUM_N) %>% full_join(pbdata %>% filter(Year == "18") %>% select(Plot, SIR)))
  summary(m1)
  
  m1 = lm(NTm ~ WORM_N + HopperN , data=pbdata %>% filter(Year == "17") %>% select(Plot, WORM_N, HopperN, AP_N, LUM_N) %>% full_join(pbdata %>% filter(Year == "18") %>% select(Plot, NTm)))
  summary(m1)
  
}
