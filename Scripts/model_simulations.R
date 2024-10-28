# Model for the cluster

require(tidyverse)
require(deSolve) # version 1.21

# How many replicates do you want to run in this code file?
NTOT = 100

# How many years do you want to simulate the experiments?
simyear = 5

# Load in the required functions ----

# A temperature function that takes day of the year (doy) and returns the temperature in Kelvin (calculated in extract_data_model_analysis.R)

LTtemp = function(doy){
  
  -12.7084*cos(2*3.14/365*doy-0.3758)+281.8291
  
}

# The complex model
singlemodel <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    TEMP = LTtemp(t %% 365) # Calculate temperature via day of the year
    
    # Model of earthworm growth. From ASA Johnston
    A_W = exp(-Ea/Kappa*(1/TEMP-1/Tref_W)) 
    
    # Model of temperature sensitive plant growth from FENG 1990
    A_P = exp(-B/TEMP)/(1 + (B/(D-B))*exp(D*(1/(Tref_P) - 1/TEMP)))
    
    # Modeled microbial dynamics MIMICS
    tempC = (TEMP-273.15)
    Vlm = exp(Vslope*tempC + Vint)*Vlm_mod
    Vsm = exp(Vslope*tempC + Vint)*Vsm_mod
    Klm = exp(Kslope*tempC + Kint)*Klm_mod
    Ksm = exp(Kslope*tempC + Kint)*Ksm_mod
    
    dL = tp*P*P + (1-SUEh)*A_W*Vhp*H*P + th*H*H + tw*W*W - Vlm*L*M/(Klm + M) - A_W*Vlw*L*W - l*L
    
    dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - tm*M - SUEwm*A_W*Vsw*W*M
    
    dW = SUEwl*A_W*Vlw*L*W + SUEws*A_W*Vsw*S*W + SUEwm*A_W*W*Vsw*M - tw*W*W
    
    dR = SUErl*A_W*Vlr*L*R + SUErm*A_W*R*Vsr*M - tr*R*R
    
    dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - A_P*Vpf*N*P/(Kpf+N)
    
    dS = tm*M + (1-SUEwl)*A_W*Vlw*L*W - Vsm*S*M/(Ksm + M) - SUEws*A_W*Vsw*S*W + fi*N - fo*S
    
    dP = A_P*Vpf*N*P/(Kpf+N) - tp*P*P - A_W*Vhp*H*P
    
    dH = SUEh*A_W*Vhp*H*P - th*H*H
    
    dy = c(dP, dL, dM, dW, dR, dN, dS, dH)
    
    return(list(c(dy)))
    
  }
  )
}

source("Scripts/parameters.R")

# Add the isopods (R) to the yint vector:
yint = c(yint[1:4], R = 9.63, yint[5:7])

# Add the isopods (R) to the params vector:
params = c(params, c(Vlr =2.400000e-06, Vsr = 4.100000e-05, SUErl = 2e-02, SUErm = 3e-01, tr = 1e-05)) # THESE ARE THE SAME AS FOR EARTHWORMS --- NEED TO UPDATE FOR WOODLICE!

# Simulate to produce a stable equilibrium for the new temperature values in the model:

yts = 2000

stablerun = ode(y=yint,times = seq(1, 365*yts,1), func=singlemodel, parms=params)

# Check to make sure the equilibrium is stable!
if(dim(stablerun)[1] == 365*yts){
  ystable = stablerun[(365*yts),-1]
  
  ERRRRR = max(abs(stablerun[(365*yts),-1] - stablerun[(365*(yts-1)),-1])) > 1e-4
} ; print(ERRRRR) # Needs to be FALSE!

stablerun %>% data.frame() %>% tibble() %>%
  filter(time > 365*1990) %>%
  pivot_longer(!time) %>% ggplot(aes(x = time, y = value)) + geom_line() + facet_wrap(~name, scales = "free_y")

ystable

# Simulate the litter removal and addition treatments:

litter_removal <- data.frame(var = "L", # Affected pool is litter
                        time =  685, # Modify it at day of the year 320, which is Nov. 15th. Do this in the second year of simulation, so 320 + 365 = 685
                        value = 0, # Multiply by zero to remove it all.
                        method = "mult") # Multiply



litter_removal_sim <- ode(y=ystable, # stable equilibrium
                          times = 1:(365*3), # Three years 
                          func=singlemodel, # Model to simulate
                          parms=params,
                          events = list(data=litter_removal))

litter_addition <- data.frame(var = "L", # Affected pool is litter
                             time =  685, # Modify it at day of the year 320, which is Nov. 15th. Do this in the second year of simulation, so 320 + 365 = 685
                             value = 2, # Multiply by two to double the litter amount.
                             method = "mult") # Multiply


litter_addition_sim <- ode(y=ystable, # stable equilibrium
                          times = 1:(365*3), # Three years 
                          func=singlemodel, # Model to simulate
                          parms=params,
                          events = list(data=litter_addition))


litter_nothing <- data.frame(var = "L", # Affected pool is litter
                             time =  685, # Modify it at day of the year 320, which is Nov. 15th. Do this in the second year of simulation, so 320 + 365 = 685
                             value = 1, # Multiply by one to have no effect on the litter amount.
                             method = "mult") # Multiply

litter_nothing_sim <- ode(y=ystable, # stable equilibrium
                           times = 1:(365*3), # Three years 
                           func=singlemodel, # Model to simulate
                           parms=params,
                           events = list(data=litter_nothing))

# Join together the two simulations:
output = litter_removal_sim %>% data.frame() %>% tibble() %>%
  pivot_longer(!time) %>% 
  mutate(Treatment = "Removal") %>%
  bind_rows(
    litter_addition_sim %>% data.frame() %>% tibble() %>%
      pivot_longer(!time) %>% 
      mutate(Treatment = "Addition")
  ) %>%
  bind_rows(
    litter_nothing_sim %>% data.frame() %>% tibble() %>%
      pivot_longer(!time) %>% 
      mutate(Treatment = "Nothing")
  ) %>% 
  inner_join(
    tibble(name = c("H", "L", "M", "N", "P", "R", "S", "W"),
           namefull = c("Herbivore", "Litter", "Microbes", "Inorganic N", "Plant", "Woodlice", "Soil", "Earthworm")), by = join_by(name)
  ) 

# Plot the raw simulations:
output %>%
  ggplot(aes(x = time, y = value, color = Treatment)) + geom_line() + facet_wrap(~namefull, scales = "free_y")

# Plot the treatment percent change:
output %>%
  pivot_wider(names_from = Treatment, values_from = value) %>%
  mutate(`Removal change` = 100*(Removal-Nothing)/Nothing,
         `Addition change` = 100*(Addition-Nothing)/Nothing) %>% filter(time > 650) %>%
  select(!Removal & !Addition & !Nothing & !name) %>%
  pivot_longer(contains("change")) %>%
  ggplot(aes(x = time, y = value, color = name)) + geom_line() + facet_wrap(~namefull, scales = "free_y") + ylab("Change (%)") + xlab("Days of simulation") +
  geom_vline(xintercept = c(685, 685 + 365), linetype = 2)
  
