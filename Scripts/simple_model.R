# Analysis of the simple model and convenience code for plotting the results of the complex model in the same format. 

#The complex model output are available in the data provided, but were run separately using other scripts on a cluster.

library(FME) # version 1.3.5
library(tidyverse)

# How many replicates do you want?
REPS = 10000

# Show the equilibrium of the complex model to help with parameter conversion for these models:

Equilibrium = c(74.6263897,16.2370121,6.7150554,10.3387229,0.1480483,80.5771975,0.0134762)
names(Equilibrium) = c("P","L" ,"M","W","N","S","H")
rm(Equilibrium)

# Analysis of simple DC model ----

# A function that calculates the equilibrium of the simple model give input parameters.
# nvec = the list of state variables that matches the order of equilibrium expressions in the following models
# m1 = no animals in the model
# m2 = Detritivores only in the model
# m3 = Herbivores only in the model
# m4 = Herbivores and detritivores both in the model
test <-function(parms){
  with(as.list(c(parms)),{
    
    out = data.frame(
      nvec = c("P", "N", "L", "H", "W"),
      m1 = c((IN*(k + l)*Vpn)/(tp*(k*q + l*q + l*Vpn)),(IN*(k + l))/(k*q + l*q + l*Vpn),(IN*Vpn)/(k*q + l*q + l*Vpn),0,0),
      m2 = c((IN*Vpn*(k + l + Vwl - ew*Vwl))/(tp*(k*q + l*q + l*Vpn + q*Vwl - ew*q*Vwl)),
                            (IN*(k + l + Vwl - ew*Vwl))/(k*q + l*q + l*Vpn + q*Vwl - ew*q*Vwl),(IN*Vpn)/(k*q + l*q + l*Vpn + q*Vwl - ew*q*Vwl),
                            0,
                            (ew*IN*Vpn*Vwl)/(tw*(k*q + l*q + l*Vpn + q*Vwl - ew*q*Vwl))),
                     
                     m3 = c((IN*(k + l)*Vpn)/(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn),
                            (IN*(k + l)*(tp + Vhp))/(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn),
                            (IN*(tp*Vpn + eh*Vhp*Vpn))/(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn),
                            (eh*IN*(k + l)*Vhp*Vpn)/(th*(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn)),0),
                     
                     m4 = c((IN*Vpn*(k + l + Vwl - ew*Vwl))/
                              (k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn + q*tp*Vwl - ew*q*tp*Vwl + q*Vhp*Vwl - ew*q*Vhp*Vwl),
                            (IN*(tp + Vhp)*(k + l + Vwl - ew*Vwl))/
                              (k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn + q*tp*Vwl - ew*q*tp*Vwl + q*Vhp*Vwl - ew*q*Vhp*Vwl),
                            (IN*(tp + eh*Vhp)*Vpn)/(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn + q*tp*Vwl - ew*q*tp*Vwl + q*Vhp*Vwl - ew*q*Vhp*Vwl),
                            (eh*IN*Vhp*Vpn*(k + l + Vwl - ew*Vwl))/
                              (th*(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn + q*tp*Vwl - ew*q*tp*Vwl + q*Vhp*Vwl - ew*q*Vhp*Vwl)),
                            (ew*IN*(tp + eh*Vhp)*Vpn*Vwl)/
                              (tw*(k*q*tp + l*q*tp + k*q*Vhp + l*q*Vhp + l*tp*Vpn + eh*l*Vhp*Vpn + q*tp*Vwl - ew*q*tp*Vwl + q*Vhp*Vwl - ew*q*Vhp*Vwl)))
      )
    
    # Set the best fitting criteria based on reasonable relationships for an old-field..NOT FIT TO ACTUAL DATA!
    if(all(c(out[out$nvec == "N",-1] < out[out$nvec == "P",-1],
          out[out$nvec == "N",-1] < out[out$nvec == "L",-1],
          out[out$nvec == "H",-1] < 100*out[out$nvec == "P",-1],
          out[out$nvec == "W",-1] < 100*out[out$nvec == "L",-1]))){
      
      Best = 1
    }else{
      Best = 0
    }
    
    out[,"Best"] = Best
    rm(Best)
    
    return(out)
    
  })
  
}

# Generate lists to save the data for 10,000 parameter sets
outlist <- vector(mode = "list", length = REPS)
paramlist <- vector(mode = "list", length = REPS)

# Run 10,000 random parameter sets
t1 = Sys.time()
for(i in 1:REPS){
  
  # Randomly draw the new parameters: See notes for how they compare to the complex model
  params = c(Vpn = rlnorm(1,meanlog = log(0.3794274), sdlog = 0.3536), # Equals A_P*Vpf*P/(Kpf + N*) from the complex model
             tp = rlnorm(1,meanlog = log(0.000005/74.6263897), sdlog = 0.3536), # Equals tp/P* in the complex model
             th = rlnorm(1,meanlog = log(1484.098), sdlog = 0.3536), # Equals th/H* in the complex model
             Vhp = rlnorm(1,meanlog = log(0.0001888764), sdlog = 0.3536), # Equals A_W*Vhp*H calculated at 25C in the complex model
             IN = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536), # Same as the complex model
             q = rlnorm(1,meanlog = log(0.1), sdlog = 0.3536),# Same as the complex model
             k = rlnorm(1,meanlog = log(0.005713002), sdlog = 0.3536), # Equals: Vlm*M/(Klm + M) + A_W*Vlw*W, where Vlm and Klm are modified by temperature so they are calculated at 25C
             Vwl = rlnorm(1,meanlog = log(3.477668e-05), sdlog = 0.3536), # Equals A_W*Vwl*W calculated at 25C in the complex model
             tw = rlnorm(1,meanlog = log(9.672375e-07), sdlog = 0.3536), # Equals tw/W* in the complex model
             l = rlnorm(1,meanlog = log(1e-4), sdlog = 0.3536), # Same as the complex model
             eh = rlnorm(1,meanlog = log(0.7), sdlog = 0.3536), # Same as the complex model
             ew = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536)) # Same as SUEwl in the complex model
  
  # Save the new parameters
  paramlist[[i]] = c(params, N = i)
  
  # Calculate the equilibrium and save it
  outlist[[i]] = cbind(test(params), N = i)
 
  # Create a update note for the user to see the progress
  if(i %% 500 == 0){
    print(paste("Done", i, "in:"))
    print(round(Sys.time() - t1))
    }
}

# Bind the parameter sets
outlist2 = do.call("rbind", outlist)

# Delete animal pools
outlist2 = outlist2[!(outlist2$nvec %in% c("H", "W")),]

# Calculate effects
outlist2[,"WE"] = (outlist2$m2 - outlist2$m1) # Detritivore effect
outlist2[,"HE"] = (outlist2$m3 - outlist2$m1) # Herbivore effect
outlist2[,"IE"] = (outlist2$m4 - outlist2$m3 - outlist2$m2 + outlist2$m1) # Interaction effect
outlist2[,"IEpred"] = outlist2[,"WE"] + outlist2[,"HE"] # Linear combination of herbivore and detritivore effects
outlist2[,"IEacc"] = outlist2$m4 - outlist2$m1 # Combined effect

# Plotting preparation
outlist2[,"ncol"] = ifelse(outlist2$nvec == "P", "#009E73",
                           ifelse(outlist2$nvec == "L", "#E69F00",
                                  ifelse(outlist2$nvec == "N", "#56B4E9",
                                  "black")))

outlist2[,"npch"] = ifelse(outlist2$nvec == "P", 1,
                           ifelse(outlist2$nvec == "L", 2,
                                  ifelse(outlist2$nvec == "N", 3,
                                         4)))

outlist2[,"ncex"] = ifelse(outlist2$Best == 1, 1, 0.5)

write_csv(outlist2, "Data/simplemodel_DC_Nov2020.csv")

# Analysis of simple LV model ----

# A function that runs the simple model with Lotka-Volterra functional responses
# nvec = the list of state variables that matches the order of equilibrium expressions in the following models
# m1 = no animals in the model
# m2 = Detritivores only in the model
# m3 = Herbivores only in the model
# m4 = Herbivores and detritivores both in the model

# A function to run the models:
test2 <- function(Ni){
  k = 1 
  
  while(k == 1){ # Draw parameters until a stable version, as judged by the principle eigenvalue of the Jacobian matrix, is found
    
    # Randomly draw the new parameters: See notes for how they compare to the complex model
    params = c(Vpn = rlnorm(1,meanlog = log(0.005084359), sdlog = 0.3536), # Equals A_P*Vpf/(Kpf + N*) from the complex model
               tp = rlnorm(1,meanlog = log(0.000005/74.6263897), sdlog = 0.3536), # Equals tp/P* in the complex model
               th = rlnorm(1,meanlog = log(20), sdlog = 0.3536), # Equals th/H* in the complex model
               Vhp = rlnorm(1,meanlog = log(0.01401555), sdlog = 0.3536), # Equals A_W*Vhp calculated at 25C in the complex model
               IN = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536), # Same as the complex model
               q = rlnorm(1,meanlog = log(0.1), sdlog = 0.3536),# Same as the complex model
               k = rlnorm(1,meanlog = log(0.005713002), sdlog = 0.3536), # Equals: Vlm*M/(Klm + M) + A_W*Vlw*W, where Vlm and Klm are modified by temperature so they are calculated at 25C
               Vwl = rlnorm(1,meanlog = log(3.363731e-06), sdlog = 0.3536), # Equals A_W*Vwl calculated at 25C in the complex model
               tw = rlnorm(1,meanlog = log(9.672375e-07), sdlog = 0.3536), # Equals tw/W* in the complex model
               l = rlnorm(1,meanlog = log(1e-4), sdlog = 0.3536), # Same as the complex model
               eh = rlnorm(1,meanlog = log(0.7), sdlog = 0.3536), # Same as the complex model
               ew = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536)) # Same as SUEwl in the complex model
    
    # Model 1 
    m1 = with(as.list(c(params)),{
      
      cbind(c((-(k*q*tp) - l*q*tp - sqrt(k + l)*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(2.*l*tp*Vpn),
              (-((q*tp)/Vpn) - (k*q*tp)/(l*Vpn) - (sqrt(k + l)*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(l*Vpn))/(2.*Vpn),
              (((q^2)*tp)/Vpn + (k*(q^2)*tp)/(l*Vpn) + 2*IN*Vpn + 
                 (sqrt(k + l)*q*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(l*Vpn))/(2.*l*Vpn)),
            c((-(k*q*tp) - l*q*tp + sqrt(k + l)*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(2.*l*tp*Vpn),
              (-((q*tp)/Vpn) - (k*q*tp)/(l*Vpn) + (sqrt(k + l)*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(l*Vpn))/(2.*Vpn),
              (((q^2)*tp)/Vpn + (k*(q^2)*tp)/(l*Vpn) + 2*IN*Vpn - 
                 (sqrt(k + l)*q*sqrt(tp)*sqrt(k*(q^2)*tp + l*(q^2)*tp + 4*IN*l*(Vpn^2)))/(l*Vpn))/(2.*l*Vpn)))
      
    }
    )
    
    if(any(apply(m1 >= 0,2,all))){
      m1 = c(m1[,apply(m1 >= 0,2,all)],0,0)
    }else{
      m1 = c(0,0,0,0,0)
    }
    
    
    P = m1[1]
    Iorg = m1[2]
    L = m1[3]
    H = m1[4]
    W = m1[5]
    m1[6] = with(as.list(c(params)),{
      max(Re(eigen(cbind(c(-2*P*tp + Iorg*Vpn,P*Vpn,0),c(-(Iorg*Vpn),-q - P*Vpn,k),c(2*P*tp,0,-k - l)))$values))
    } 
    )
    
    # Model 2
    m2 = with(as.list(c(params)),{
      c(
        (Vpn*(-(l*tw) + ew*IN*Vwl))/(ew*q*tp*Vwl),
        
        (-(l*tw) + ew*IN*Vwl)/(ew*q*Vwl),
        
        tw/(ew*Vwl),
        
        0,
        
        (-((l^2)*(tw^2)*(Vpn^2)) + ew*k*(q^2)*tp*tw*Vwl + ew*l*(q^2)*tp*tw*Vwl + 2*ew*IN*l*tw*(Vpn^2)*Vwl - (ew^2)*(IN^2)*(Vpn^2)*(Vwl^2))/((-1 + ew)*ew*(q^2)*tp*tw*(Vwl^2)))
    })
    
    P = m2[1]
    Iorg = m2[2]
    L = m2[3]
    H = m2[4]
    W = m2[5]
    
    m2[6] = with(as.list(c(params)),{
      max(Re(eigen(rbind(c(-2*P*tp + Iorg*Vpn,P*Vpn,0,0),c(-(Iorg*Vpn),-q - P*Vpn,k + (1 - ew)*Vwl*W,(1 - ew)*L*Vwl),c(2*P*tp,0,-k - l - Vwl*W,tw - L*Vwl),
                         c(0,0,ew*Vwl*W,-tw + ew*L*Vwl)))$values))
    } 
    )
    
    # Model 3
    m3 = with(as.list(c(params)),{
      
      c(
        th/(eh*Vhp),
        
        (-(l*(th^2)*tp) + eh*l*(th^2)*tp + (eh^2)*IN*k*(Vhp^2) + (eh^2)*IN*l*(Vhp^2))/((eh^2)*Vhp*(k*q*Vhp + l*q*Vhp + l*th*Vpn)),
        
        -((th*(-(q*th*tp) + eh*q*th*tp - (eh^2)*IN*Vhp*Vpn))/((eh^2)*Vhp*(k*q*Vhp + l*q*Vhp + l*th*Vpn))),
        
        (-(eh*k*q*th*tp*Vhp) - eh*l*q*th*tp*Vhp - l*(th^2)*tp*Vpn + (eh^2)*IN*k*(Vhp^2)*Vpn + (eh^2)*IN*l*(Vhp^2)*Vpn)/((eh^2)*(Vhp^2)*(k*q*Vhp + l*q*Vhp + l*th*Vpn)),
        
        0
        
      )
    })
    
    P = m3[1]
    Iorg = m3[2]
    L = m3[3]
    H = m3[4]
    W = m3[5]
    
    m3[6] = with(as.list(c(params)),{
      max(Re(eigen(rbind(c(-2*P*tp - H*Vhp + Iorg*Vpn,P*Vpn,0,-(P*Vhp)),c((1 - eh)*H*Vhp - Iorg*Vpn,-q - P*Vpn,k,(1 - eh)*P*Vhp),c(2*P*tp,0,-k - l,th),
                         c(eh*H*Vhp,0,0,-th + eh*P*Vhp)))$values))
    } 
    )
    
    # Model 4
    m4 = with(as.list(c(params)),{
      
      c(
        th/(eh*Vhp),
        
        (-(l*tw) + ew*IN*Vwl)/(ew*q*Vwl),
        
        tw/(ew*Vwl),
        
        (-(eh*l*tw*Vhp*Vpn) - ew*q*th*tp*Vwl + eh*ew*IN*Vhp*Vpn*Vwl)/(eh*ew*q*(Vhp^2)*Vwl),
        
        ((eh^2)*k*q*tw*(Vhp^2) + (eh^2)*l*q*tw*(Vhp^2) + (eh^2)*l*th*tw*Vhp*Vpn - ew*q*(th^2)*tp*Vwl + eh*ew*q*(th^2)*tp*Vwl - (eh^2)*ew*IN*th*Vhp*Vpn*Vwl)/((eh^2)*(-1 + ew)*q*tw*(Vhp^2)*Vwl))
    })
    
    P = m4[1]
    Iorg = m4[2]
    L = m4[3]
    H = m4[4]
    W = m4[5]
    
    m4[6] = with(as.list(c(params)),{
      max(Re(eigen(rbind(c(-2*P*tp - H*Vhp + Iorg*Vpn,P*Vpn,0,-(P*Vhp),0),c((1 - eh)*H*Vhp - Iorg*Vpn,-q - P*Vpn,k + (1 - ew)*Vwl*W,(1 - eh)*P*Vhp,(1 - ew)*L*Vwl),c(2*P*tp,0,-k - l - Vwl*W,th,tw - L*Vwl),c(eh*H*Vhp,0,0,-th + eh*P*Vhp,0),c(0,0,ew*Vwl*W,0,-tw + ew*L*Vwl)))$values))
    } 
    )
    
    if(all(c(m1[c(1,2,3,6)]*c(1,1,1,-1) > 0,m2[c(1,2,3,5,6)]*c(1,1,1,1,-1) > 0,m3[c(1,2,3,4,6)]*c(1,1,1,1,-1) > 0,m4*c(1,1,1,1,1,-1) > 0))) k = 2
    
  }
  
  output = data.frame(nvec =c("P","N","L","H","W","rmax"),
                       m1 = m1,
                       m2 = m2,
                       m3 = m3,
                       m4 = m4,
                       N = Ni)
  
  if(all(c(output[output$nvec == "N",-1] < output[output$nvec == "P",-1],
           output[output$nvec == "N",-1] < output[output$nvec == "L",-1],
           output[output$nvec == "H",-1] < 100*output[output$nvec == "P",-1],
           output[output$nvec == "W",-1] < 100*output[output$nvec == "L",-1]))){
    
    Best = 1
  }else{
    Best = 0
  }
  
  output[,"Best"] = Best
  
  return(list(output, c(params, N = Ni)))
}

# Run the model
out2 = lapply(seq(1,REPS,1), FUN = test2)

# Reorganize the data
out3 = vector(mode = "list", length = REPS)

for(j in 1:REPS){
  out3[[j]] = out2[[j]][[1]]
}

out4 = do.call("rbind", out3)

out4 = out4[!(out4$nvec %in% c("H", "W", "rmax")),]

# Calculate effects
out4[,"WE"] = (out4$m2 - out4$m1) # Detritivore effect
out4[,"HE"] = (out4$m3 - out4$m1) # Herbivore effect
out4[,"IE"] = (out4$m4 - out4$m3 - out4$m2 + out4$m1) # Interaction effect
out4[,"IEpred"] = out4[,"WE"] + out4[,"HE"] # Linear combination of herbivore and detritivore effects
out4[,"IEacc"] = out4$m4 - out4$m1 # Combined effect


# Plotting preparation
out4[,"ncol"] = ifelse(out4$nvec == "P", "#009E73",
                           ifelse(out4$nvec == "L", "#E69F00",
                                  ifelse(out4$nvec == "N", "#56B4E9",
                                         "black")))

out4[,"npch"] = ifelse(out4$nvec == "P", 1,
                           ifelse(out4$nvec == "L", 2,
                                  ifelse(out4$nvec == "N", 3,
                                         4)))

out4[,"ncex"] = ifelse(out4$Best == 1, 1, 0.5)

write_csv(out4, "Data/simplemodel_LV_Nov2020.csv")

# Clean out the earlier analysis -----

rm(test, test2)

# Load in the simple model data if necessary -----

outlist2 = read_csv("Data/simplemodel_DC_Nov2020.csv")
out4 = read_csv("Data/simplemodel_LV_Nov2020.csv")

# Plot the most complex model ----

# Load in the equilibrium cluster data (This data is calculated using the script "complex_model_non-eqm_to_cluster.R")

if(F){ # Only necessary if loading data directly from cluster. Provided data is loaded below
  dirtoload = "Model_eqm_Jun2020/"
  
  ftoload = list.files(dirtoload)
  
  ftoload = ftoload[grepl("_eqm.csv", x = ftoload)]
  
  listOfDataFrames <- vector(mode = "list", length = length(ftoload))
  
  for(ii in 1:length(ftoload)){
    listOfDataFrames[[ii]] <- read.csv(paste0(dirtoload,ftoload[ii]))
  }
  
  data2 <- do.call("rbind", listOfDataFrames)
  
  
  ftoload = list.files(dirtoload)
  
  ftoload = ftoload[grepl("_param.csv", x = ftoload)]
  
  listOfDataFrames <- vector(mode = "list", length = length(ftoload))
  
  for(ii in 1:length(ftoload)){
    listOfDataFrames[[ii]] <- read.csv(paste0(dirtoload,ftoload[ii]))
  }
  
  paramsEQM <- do.call("rbind", listOfDataFrames)
  
  apply(paramsEQM,2, sd)
  
  ftoload = list.files(dirtoload)
  
  ftoload = ftoload[grepl("_noneqm.csv", x = ftoload)]
  
  listOfDataFrames <- vector(mode = "list", length = length(ftoload))
  
  for(ii in 1:length(ftoload)){
    listOfDataFrames[[ii]] <- read.csv(paste0(dirtoload,ftoload[ii]))
  }
  
  data2_noneqm <- do.call("rbind", listOfDataFrames)
  
  rm(listOfDataFrames,ftoload)
  
  data2aa = data2 %>%
    separate(X, into = c(NA, "X"), sep = 7) %>%
    filter(Stable == 1) %>% 
    rename(Treatment = X) %>%
    mutate(Type = "EQM")
  
  data2_noneqm_aa = data2_noneqm %>%
    left_join(
      data2aa %>% select(ID, Treatment, Stable) %>% distinct()
    ) %>%
    mutate(Type = "NON-EQM")
  
  colnames(data2aa)
  colnames(data2_noneqm_aa)
  
  data2 = data2aa %>%
    bind_rows(
      data2_noneqm_aa
    )
  
  write_rds(data2, "Data/complex_model_10000_eqm.rds")
  
}

data2 = read_rds("Data/complex_model_10000_eqm.rds") %>% distinct()

# Analyze the data
data2a = data2 %>% as_tibble() %>% 
  gather(-Treatment, -Stable, -ID, - Type, key = StateVar, value = biomass) %>%
  filter(Stable == 1) %>% select(-Stable)
data2b = data2a %>%
  group_by(ID, Type) %>% summarize(N = n()) %>% filter(N == 28) %>%
  select(ID, Type) %>%
  left_join(
    data2a
  ) %>%
  distinct() %>%
  spread(key = Treatment, value = biomass) %>%
  filter(!StateVar %in% c("H", "W"))

data2c = data2b %>%
  mutate(WE = W - N,
         HE = H - N,
         IE = HW - H - W + N) %>%
  mutate(IEacc = HW -N,
         IEpred = HE + WE) %>%
  left_join(
    data.frame(StateVar = c("P", "L", "N", "S", "M"),
               ncol = as.character(c("#009E73","#E69F00","#56B4E9","#F0E442", "#0072B2")),
               npch = c(1,2,3,4,5))
  ) %>%
  left_join(
    data2 %>% as_tibble() %>% 
      mutate(Best = ifelse(N < P & N < L & H < 100*P & W < 100*L, 1,0)) %>%
      group_by(ID, Type) %>% summarize(N = sum(Best)) %>%
      mutate(ncex = ifelse(N == 4, 1, 0.5)) %>%
      select(ID, Type, ncex)
  )

data2c$ncol = as.character(data2c$ncol)

# Split the two types
data2c_noneqm = data2c %>% filter(Type == "NON-EQM") %>% select(-Type)
data2c = data2c %>% filter(Type == "EQM") %>% select(-Type)

plot((abs(IE+1e-6))~(abs(WE+1e-6)), data = data2c, log = 'xy', col = ncol, pch = npch,
     xlab= "Detritivore Effect (log|DE|)", ylab = "Interaction Effect (log|IE|)", cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
legend("topleft", legend = "E", bty = "n")

plot((abs(IE+1e-6))~(abs(HE+1e-6)), data = data2c, log = 'xy', col = ncol, pch = npch,
     xlab= "Herbivore Effect (log|HE|)", ylab = "", cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
legend("topleft", legend = "F", bty = "n")

# Plot all four models together ----

# A. Simple model, donor-controlled
# B. Simple model, Lotka-Volterra
# C. Complex model, equilibrium
# D. Complex model, non-equilibrium

outfinal3 = read_rds("Data/TrueLinearComplex.rds") # Load in the non-equilibrium cluster data

outfinal3$ncol = as.character(outfinal3$ncol)

png("Plots/Bothmodel2.png", width = 12, height = 8, units = "in", res = 600)

par(oma=c(2,2,0,0))
par(mfrow=c(2,3), mar = c(3,3,1,1))
plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outlist2, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", type = "n", main = "Simple: Donor controlled",
     axes = F, ylim = c(1e-11, 1e36), xlim = c(1e-11, 1e36))
axis(1, at = c(1e-11,1e-3,1e5,1e13, 1e21, 1e31), labels = expression(10^-11,10^-3,10^5,10^13, 10^21, 10^31))
axis(2, at = c(1e-11,1e-3,1e5,1e13, 1e21, 1e31), labels = expression(10^-11,10^-3,10^5,10^13, 10^21, 10^31))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outlist2, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "A", bty = "n")
legend("bottomright", legend = c("Plant", "Litter", "Inorganic N", "Soil", "Microbe"),
       col = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), pch = 1:5, bty = "n", title = "State Variable")


plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = out4, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", type = "n", main = "Simple: Type I",
     axes = F, ylim = c(1e-11, 1e36), xlim = c(1e-11, 1e36))
axis(1, at = c(1e-11,1e-3,1e5,1e13, 1e21, 1e31), labels = expression(10^-11,10^-3,10^5,10^13, 10^21, 10^31))
axis(2, at = c(1e-11,1e-3,1e5,1e13, 1e21, 1e31), labels = expression(10^-11,10^-3,10^5,10^13, 10^21, 10^31))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = out4, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "B", bty = "n")

plot(1,1, type = "n",axes = F, ylab = "", xlab = "", main = "Complex Model Animal Biomass")
legend("topleft", legend = "C", bty = "n")

plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = data2c, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", type = "n", main = "Complex: Equilibrium",
     axes = F, ylim = c(1e-11, 1e5), xlim = c(1e-11, 1e5))
axis(1, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
axis(2, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = data2c, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "D", bty = "n")

plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = data2c_noneqm, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", type = "n", main = "Complex: Non-equilibrium",
     axes = F, ylim = c(1e-11, 1e5), xlim = c(1e-11, 1e5))
axis(1, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
axis(2, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = data2c_noneqm, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "E", bty = "n")

plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outfinal3, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", 
     type = "n", main = "Complex: Field simulation",
     axes = F, ylim = c(1e-11, 1e5), xlim = c(1e-11, 1e5))
axis(1, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
axis(2, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outfinal3, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "F", bty = "n")


mtext(text="Combined Effect (log|x|)",side=2,line=0,outer=TRUE,cex=1)
mtext(text="Herbivore + Detritivore Effect (log|x|)",side=1,line=0,outer=TRUE,cex=1)
dev.off()


# Another version of the both model plot -----

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

scientific3 <- function(x){
  x = signif(x, digits = 1)
  ifelse(x==0, "0", gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
}

selected_pool = "P"

IEplot = outlist2 %>% 
  filter(nvec == selected_pool) %>%
  mutate(IEscale = abs(IE)+1e-6,
         ID = "Simple: Donor-controlled") %>%
  select(IEscale, ID) %>%
  bind_rows(
    out4 %>% 
      filter(nvec == selected_pool) %>%
      mutate(IEscale = abs(IE)+1e-6,
             ID = "Simple: Lotka-Volterra") %>%
      select(IEscale, ID)
    ) %>% 
  bind_rows(
    data2c %>% 
      filter(StateVar == selected_pool) %>%
      ungroup() %>%
      mutate(IEscale = abs(IE)+1e-6,
             ID = "Complex: Equilibrium") %>%
      select(IEscale, ID)
  ) %>% 
  bind_rows(
    data2c_noneqm %>% 
      filter(StateVar == selected_pool) %>%
      ungroup() %>%
      mutate(IEscale = abs(IE)+1e-6,
             ID = "Complex: Non-equilibrium") %>%
      select(IEscale, ID)
  ) %>% 
  bind_rows(
    outfinal3 %>% 
      filter(StateVar == selected_pool) %>%
      ungroup() %>%
      mutate(IEscale = abs(IE)+1e-6,
             ID = "Complex: Field simulation") %>%
      select(IEscale, ID)
  )
  
IEtext = IEplot %>% group_by(ID) %>%
  summarise(X = median(IEscale)) %>%
  mutate(Y = c(0.3, 0.5, 0.4, 0.8, 0.5)+ 0.05) %>%
  mutate(t = X) %>%
  mutate(t = scientific3(t)) %>%
  mutate(X2 = ifelse(X > 1000, X*10, X))

png(paste0("Plots/Figure5_Nov2020_",selected_pool,".png"), width = 7, height = 4, units = "in", res = 600)
IEplot %>% ggplot(aes(x = IEscale, fill = ID)) + 
  geom_density(alpha = 0.7) + 
  geom_text(aes(x = X2, y = Y, label = t, col = ID),data = IEtext, parse = T) + 
  theme_classic() + 
  scale_x_log10(labels = scientific, name = "Interaction effect onto plants") + 
  scale_fill_manual(name = "Simulation Type", values = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), breaks = c("Simple: Donor-controlled", "Simple: Lotka-Volterra","Complex: Equilibrium","Complex: Non-equilibrium","Complex: Field simulation")) + 
  scale_color_manual(guide = F, values = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), breaks = c("Simple: Donor-controlled", "Simple: Lotka-Volterra","Complex: Equilibrium","Complex: Non-equilibrium","Complex: Field simulation")) +
  ylab("Density") +
  theme(legend.position = c(0.3, 0.55),
        legend.justification = c(1, 0),
        legend.box = "horizontal")
dev.off()

# Get the abundance data from the complex model and plot it -----

outlist2

head(out4)

datahist <- data2 %>%
  filter(Stable == 1) %>% select(Treatment, W, H, Type) %>%
  bind_rows(
    read_rds("Data/TrueLinearComplex_StateVar.rds") %>% select(Treatment, H, W) %>% mutate(Type = "4-year")
  ) %>%
  left_join(
    tibble(Type = c("EQM", "NON-EQM", "4-year"),
           Type2 = c("Equilibrium", "Non-equilibrium", "Field simulation"))
  ) %>%
  left_join(
    tibble(Treatment = c("N", "H", "W", "HW"),
           Treatment2 = c("Neither", "Herbivore", "Detritivore", "Both"))
  ) %>%
  select(Treatment2, Type2, H, W) %>%
  rename(Treatment = Treatment2, Type = Type2)

datahist$Type <- factor(datahist$Type, levels = c("Equilibrium", "Non-equilibrium", "Field simulation"))
datahist$Treatment <- factor(datahist$Treatment, levels = c("Neither", "Herbivore", "Detritivore", "Both"))

scientific2 <- function(x){
  ifelse(x==0, "0", parse(text=gsub("1 %*% ","",x = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))), fixed = T)))
}

png("Plots/complex_model_HW.png", width = 6, height = 5, units = "in", res = 600)
datahist %>%
  gather(-Treatment, -Type, key = StateVar, value = Biomass) %>%
  filter(Biomass > 0) %>%
  left_join(
    tibble(StateVar = c("H", "W"),
           StateVar2 = c("Herbivore (Grasshopper)", "Detritivore (Earthworm)"))
  ) %>%
  ggplot(aes(y = Biomass, x = Treatment, fill = Type)) + geom_boxplot() + theme_classic() + scale_y_log10(labels = scientific2) + facet_wrap(.~StateVar2) + scale_fill_viridis_d()  +
  theme(legend.position = c(0.4, 0.05),
        legend.justification = c(0, 0),
        legend.box = "horizontal")
dev.off()
