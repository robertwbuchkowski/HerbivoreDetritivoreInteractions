# Analysis of the simple model and convenience code for plotting the results of the complex model in the same format. 
# Robert W. Buchkowski
# Dec. 2020
# Note: The complex model output are available in the data provided, but are produced in a different script using a computing cluster.
# Note: On line 365, you can load in the model data if you do not want to simulate them. My machine runs the entire script in about 5 minutes.

library(FME) # version 1.3.5
library(tidyverse) # version 1.2.1

# The number of replicates for the two versions of the simple model:
REPS = 10000

# Show the equilibrium of the complex model to help with parameter conversion for these models ----
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
test <-function(
  parms # A function of the model parameters
  ){
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
    
    # Produce a list of model outputs that are reasonable for an old field, these are not used in the paper.
    # Set based on reasonable relationships for an old-field..NOT FIT TO ACTUAL DATA!
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

# Generate lists to save the data for 10,000 parameter sets and those parameter sets
outlist <- vector(mode = "list", length = REPS)
paramlist <- vector(mode = "list", length = REPS)

# Run 10,000 random parameter sets
t1 = Sys.time()
for(i in 1:REPS){
  
  # Randomly draw the new parameters: See notes for how they compare to the complex model
  params = c(Vpn = rlnorm(1,meanlog = log(0.3794274), sdlog = 0.3536), # Equals A_P*Vpf*P/(Kpf + N*) from the complex model
             tp = rlnorm(1,meanlog = log(0.000005/74.6263897), sdlog = 0.3536), # Equals tp/P* in the complex model
             th = rlnorm(1,meanlog = log(14.84098), sdlog = 0.3536), # Equals th/H* in the complex model : INCREASED TWO ORDERS OF MAGNITUDE (OOM) TO MATCH EQM
             Vhp = rlnorm(1,meanlog = log(0.01888764), sdlog = 0.3536), # Equals A_W*Vhp*H calculated at 25C in the complex model : INCREASED TWO OOM TO MATCH EQM
             IN = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536), # Same as the complex model
             q = rlnorm(1,meanlog = log(0.1), sdlog = 0.3536),# Same as the complex model
             k = rlnorm(1,meanlog = log(0.005713002), sdlog = 0.3536), # Equals: Vlm*M/(Klm + M) + A_W*Vlw*W, where Vlm and Klm are modified by temperature so they are calculated at 25C
             Vwl = rlnorm(1,meanlog = log(3.477668e-03), sdlog = 0.3536), # Equals A_W*Vwl*W calculated at 25C in the complex model : INCREASED TWO OOM TO MATCH EQM
             tw = rlnorm(1,meanlog = log(9.672375e-06), sdlog = 0.3536), # Equals tw/W* in the complex model : INCREASED TWO OOM TO MATCH EQM
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

# Save range in pool sizes for m4, the model with both herbivores and detritivores for later comparison.
DCeqm = outlist2 %>% select(nvec, m4) %>%
  group_by(nvec) %>%
  summarise(lq = quantile(m4, 0.25),
            med = median(m4),
            uq = quantile(m4, 0.75))

# Delete animal pools
outlist2 = outlist2[!(outlist2$nvec %in% c("H", "W")),]

# Calculate effects
outlist2[,"WE"] = (outlist2$m2 - outlist2$m1) # Detritivore effect
outlist2[,"HE"] = (outlist2$m3 - outlist2$m1) # Herbivore effect
outlist2[,"IE"] = (outlist2$m4 - outlist2$m3 - outlist2$m2 + outlist2$m1) # Interaction effect
outlist2[,"IEpred"] = outlist2[,"WE"] + outlist2[,"HE"] # Linear combination of herbivore and detritivore effects
outlist2[,"IEacc"] = outlist2$m4 - outlist2$m1 # Combined effect

# Plotting preparation colors and symbols: Not used in the final version of the plots, but is helpful if the user wants to make scatter plots
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

rm(params)

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
               th = rlnorm(1,meanlog = log(14.84098), sdlog = 0.3536), # Equals th/H* in the complex model : INCREASED TWO OOM TO MATCH EQM
               Vhp = rlnorm(1,meanlog = log(1.401555), sdlog = 0.3536), # Equals A_W*Vhp calculated at 25C in the complex model : INCREASED TWO OOM TO MATCH EQM
               IN = rlnorm(1,meanlog = log(0.02), sdlog = 0.3536), # Same as the complex model
               q = rlnorm(1,meanlog = log(0.1), sdlog = 0.3536),# Same as the complex model
               k = rlnorm(1,meanlog = log(0.005713002), sdlog = 0.3536), # Equals: Vlm*M/(Klm + M) + A_W*Vlw*W, where Vlm and Klm are modified by temperature so they are calculated at 25C
               Vwl = rlnorm(1,meanlog = log(3.363731e-04), sdlog = 0.3536), # Equals A_W*Vwl calculated at 25C in the complex model : INCREASED TWO OOM TO MATCH EQM
               tw = rlnorm(1,meanlog = log(9.672375e-06), sdlog = 0.3536), # Equals tw/W* in the complex model : INCREASED TWO OOM TO MATCH EQM
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
  
  # Create a dataframe. This data frame also has rmax to check for stability.
  output = data.frame(nvec =c("P","N","L","H","W","rmax"),
                       m1 = m1,
                       m2 = m2,
                       m3 = m3,
                       m4 = m4,
                       N = Ni)
  
  # Produce a list of model outputs that are reasonable for an old field, these are not used in the paper.
  # Set based on reasonable relationships for an old-field..NOT FIT TO ACTUAL DATA!
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

# Run the model: Note if you choose parameter values that do not produce a stable model for all four cases then this will run forever, just FYI...
out2 = lapply(seq(1,REPS,1), FUN = test2)

# Reorganize the data to extract the model results
out3 = vector(mode = "list", length = REPS)

# Flip the list:
for(j in 1:REPS){
  out3[[j]] = out2[[j]][[1]]
}

# Join the list into a data frame
out4 = do.call("rbind", out3)

# Save equilibrium for comparisons across models
LVeqm = out4 %>% select(nvec, m4) %>%
  group_by(nvec) %>%
  summarise(lq = quantile(m4, 0.25),
            med = median(m4),
            uq = quantile(m4, 0.75)) %>%
  filter(nvec != "rmax")

# Remove the animals and the stability test results
out4 = out4[!(out4$nvec %in% c("H", "W", "rmax")),]

# Calculate effects
out4[,"WE"] = (out4$m2 - out4$m1) # Detritivore effect
out4[,"HE"] = (out4$m3 - out4$m1) # Herbivore effect
out4[,"IE"] = (out4$m4 - out4$m3 - out4$m2 + out4$m1) # Interaction effect
out4[,"IEpred"] = out4[,"WE"] + out4[,"HE"] # Linear combination of herbivore and detritivore effects
out4[,"IEacc"] = out4$m4 - out4$m1 # Combined effect


# Plotting preparation colors and symbols: Not used in the final version of the plots, but is helpful if the user wants to make scatter plots
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

# Load in the simple model data if necessary -----

outlist2 = read_csv("Data/simplemodel_DC_Nov2020.csv")
out4 = read_csv("Data/simplemodel_LV_Nov2020.csv")

# Plot the most complex model ----

# Load in the equilibrium cluster data (This data is calculated using the script "complex_model_non-eqm_to_cluster.R")
# Running the code in the following 'if' statement is only necessary if loading data directly from cluster. Provided data is loaded below
if(F){ 
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

# Load in the provided data:
data2 = read_rds("Data/complex_model_10000_eqm.rds") %>% distinct()

# Analyze the data for effects
data2a = data2 %>% as_tibble() %>% 
  gather(-Treatment, -Stable, -ID, - Type, key = StateVar, value = biomass) %>%
  filter(Stable == 1) %>% select(-Stable)
data2b = data2a %>%
  group_by(ID, Type) %>% 
  summarize(N = n()) %>% 
  filter(N == 28) %>%
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

# Split the two types: Complex: Equilibrium and Complex: Non-equilibrium
data2c_noneqm = data2c %>% filter(Type == "NON-EQM") %>% select(-Type)
data2c = data2c %>% filter(Type == "EQM") %>% select(-Type)

# Load in the field simulation data: These are just a subset of the data used for the approximate Bayesian computaiton.
outfinal3 = read_rds("Data/TrueLinearComplex.rds")
outfinal3$ncol = as.character(outfinal3$ncol)

# Take out the 48 "best fitting" simulations that were added to this one.
outfinal3 %>% filter(Best == "No") %>%
  bind_rows(
    outfinal3 %>% 
      filter(Best == "Yes") %>%
      filter(Run %in% c("2332952347464863", "7303471844607670")) # Keep two best fitting simulations because there were 2 in the original sample.
  )

# Plot the equilibria values together: Only works if you ran the models, otherwise, load the data below
comeqm = DCeqm %>%
  mutate(ID = "Simple: Donor-controlled") %>%
  mutate(Model = "S-DC") %>%
  bind_rows(
    LVeqm %>%
      mutate(ID = "Simple: Lotka-Volterra") %>%
      mutate(Model = "S-LV")
  ) %>%
  bind_rows(
    data2a %>%
      group_by(ID, Type) %>% 
      summarize(N = n()) %>% 
      filter(N == 28) %>%
      select(ID, Type) %>%
      left_join(
        data2a
      ) %>%
      distinct() %>%
      spread(key = Treatment, value = biomass) %>% 
      filter(Type == "EQM") %>% 
      filter(!(StateVar %in% c("S", "M"))) %>%
      ungroup() %>%
      select(StateVar, HW) %>%
      group_by(StateVar) %>%
      summarise(lq = quantile(HW, 0.25),
                med = median(HW),
                uq = quantile(HW, 0.75)) %>%
      mutate(ID = "Complex: Equilibrium") %>%
      mutate(Model = "Complex") %>%
      rename(nvec = StateVar)
  ) %>%
  filter(nvec != "M" & nvec != "S") %>%
  left_join(
    tibble(nvec = c("H", "L", "P", "W", "N"),
           nvec2 = c("Herbivore", "Litter", "Plant", "Detritivore", "Inorganic N")
    )
  ) %>% 
  write_csv("Data/compeqm_Nov2020.csv")

comeqm = read_csv("Data/compeqm_Nov2020.csv")

cmsv = comeqm %>%
  ggplot(aes(x = Model)) + geom_pointrange(aes(y = med, ymin = lq, ymax = uq, color = ID)) + theme_classic() + facet_wrap(.~nvec2, scales = "free") +
  scale_color_manual(name = "Simulation Type", values = c("#009E73","#F0E442","#0072B2"), breaks = c("Simple: Donor-controlled", "Simple: Lotka-Volterra","Complex: Equilibrium"))  +
  theme(legend.position = c(0.95, 0.08),
        legend.justification = c(1, 0),
        legend.box = "horizontal")

png("Plots/model_compare_state_var.png", width = 8, height = 6, units = "in", res = 600)
cmsv
dev.off()

# Plot all four models together ----

# Three functions for scientific notation in the graphs:
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

scientific2 <- function(x){
  ifelse(x==0, "0", parse(text=gsub("1 %*% ","",x = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))), fixed = T)))
}

scientific3 <- function(x){
  x = signif(x, digits = 1)
  ifelse(x==0, "0", gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
}


# This function will plot the outcomes by state variable
toplot_int <- function(selected_pool, # Which state variable (P, N, L)?
                       xlabel, # What should it be called in the axis label?
                       Yvec, # Where should the text labels be vertically?
                       Xvec, # Multiplier to adjust text labels horizontally
                       legpos = c(0.3, 0.55) # Position of the legend
                       ){
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
    mutate(Y = Yvec) %>%
    mutate(t = X) %>%
    mutate(t = scientific3(t)) %>%
    mutate(X2 = X*Xvec)
  
  lalala = IEplot %>% ggplot(aes(x = IEscale, fill = ID)) + 
    geom_density(alpha = 0.7) + 
    geom_text(aes(x = X2, y = Y, label = t, col = ID),data = IEtext, parse = T) + 
    theme_classic() + 
    scale_x_log10(labels = scientific, name = xlabel) + 
    scale_fill_manual(name = "Simulation Type", values = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), breaks = c("Simple: Donor-controlled", "Simple: Lotka-Volterra","Complex: Equilibrium","Complex: Non-equilibrium","Complex: Field simulation")) +
    scale_color_manual(guide = F, values = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), breaks = c("Simple: Donor-controlled", "Simple: Lotka-Volterra","Complex: Equilibrium","Complex: Non-equilibrium","Complex: Field simulation")) +
    ylab("Density") +
    theme(legend.position = legpos,
          legend.justification = c(1, 0),
          legend.box = "horizontal")
  return(lalala)
}

# Run the plant plot (Figure 2 in the manuscript):
a = toplot_int(selected_pool = "P",
           xlabel = "Interaction effect onto plants",
           Yvec = c(0.3, 0.5, 0.4, 0.8, 0.5)+ 0.05,
           Xvec = c(1,1,1,10,0.1))

png(paste0("Plots/Figure5_Nov2020_","P",".png"), width = 7, height = 4, units = "in", res = 600)
a
dev.off()

# Run the inorganic nitrogen plot:
aN = toplot_int(selected_pool = "N",
               xlabel = "Interaction effect onto inorganic N",
               Yvec = c(0.95, 0.85, 0.85, 0.75, 0.85)+ 0.05,
               Xvec = c(1,2,1,1,10),
               legpos = c(0.6, 0.55))

png(paste0("Plots/Figure5_Nov2020_","N",".png"), width = 7, height = 4, units = "in", res = 600)
aN
dev.off()

# Run the litter plot:
aL = toplot_int(selected_pool = "L",
               xlabel = "Interaction effect onto litter",
               Yvec = c(1.05, 0.85, 0.95, 0.85, 0.95)+ 0.05,
               Xvec = c(1,1,1,1,10),
               legpos = c(0.4, 0.55))

png(paste0("Plots/Figure5_Nov2020_","L",".png"), width = 7, height = 4, units = "in", res = 600)
aL
dev.off()

# Get the abundance data from the Complex model: field simulation and plot it -----
# This is the plot showing grasshopper and earthworm abundance: Appendix 2: Figure 1D
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

cmsv2 = datahist %>%
  gather(-Treatment, -Type, key = StateVar, value = Biomass) %>%
  filter(Biomass > 0) %>%
  left_join(
    tibble(StateVar = c("H", "W"),
           StateVar2 = c("Herbivore (Grasshopper)", "Detritivore (Earthworm)"))
  ) %>%
  ggplot(aes(y = Biomass, x = Treatment, fill = Type)) + geom_boxplot() + theme_classic() + scale_y_log10(labels = scientific2) + facet_wrap(.~StateVar2) +
  theme(legend.position = c(0.4, 0.05),
        legend.justification = c(0, 0),
        legend.box = "horizontal") + 
  scale_fill_manual(name = "Simulation Type", values = c("#009E73","#56B4E9","#E69F00")) + ylab("Biomass (log scale)")

png("Plots/complex_model_HW.png", width = 6, height = 5, units = "in", res = 600)
cmsv2
dev.off()


# Combine the N, L, abundance equilibrium, and abundance in complex model plots for the supplemental ----
png("Plots/supplementary_figure_model_comparison.png", width = 14, height = 8, units = "in", res = 600)
cowplot::plot_grid(aN,aL, cmsv, cmsv2, labels = "AUTO") # Appendix S2: Figure S1
dev.off()
