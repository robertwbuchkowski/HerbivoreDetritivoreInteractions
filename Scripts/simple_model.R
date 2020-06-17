# Analysis of the simple model and convenience code for plotting the results of the complex model in the same format. 

#The complex model output are available in the data provided, but were run separately using other scripts on a cluster.

library(FME) # version 1.3.5
library(tidyverse)

# Analysis of simple DC model ----

params = c(Vpn = 1,
           tp = 1,
           th = 1,
           Vhp = 1,
           IN = 1,
           q = 0.5,
           k = 1,
           Vwl = 1,
           tw = 1,
           tp = 1,
           l = 0.5,
           eh = 0.8,
           ew = 0.5)

# A function that calculates the equilibrium of the simple model give input parameters.

# m1 = no animals
# m2 = Detritivores
# m3 = Herbivores
# m4 = Herbivores and detritivores

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

# Generate lists to save the data
outlist <- vector(mode = "list", length = 10000)
paramlist <- vector(mode = "list", length = 10000)
simlist <- vector(mode = "list", length = 10000)

# Create dynamic models 1-4 (as above) to simulate the non-equilibrium dynamics

simpleDCmodel1 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg - tp*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg
    dL = tp*P - k*L - l*L
    
    return(list(c(dP, dIorg, dL)))
  })
}

simpleDCmodel2 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg - tp*P;
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg + (1 - ew)*Vwl*L
    dL = tw*W + tp*P - k*L - l*L - Vwl*L
    dW = ew*Vwl*L - tw*W
    
    return(list(c(dP, dIorg, dL, dW)))
  })
}

simpleDCmodel3 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg - tp*P - Vhp*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg + (1 - eh)*Vhp*P
    dL = th*H + tp*P - k*L - l*L
    dH = eh*Vhp*P - th*H
    
    return(list(c(dP, dIorg, dL, dH)))
  })
}

simpleDCmodel4 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg - tp*P - Vhp*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg + (1 - eh)*Vhp*P + (1 - ew)*Vwl*L
    dL = tw*W + th*H + tp*P - k*L - l*L - Vwl*L
    dH = eh*Vhp*P - th*H
    dW = ew*Vwl*L - tw*W
    
    return(list(c(dP, dIorg, dL, dH, dW)))
  })
}


# Run 10,000 random parameter sets
t1 = Sys.time()
for(i in 1:10000){
  
  params = c(Vpn = rlnorm(1,meanlog = 1, sdlog = 10),
             tp = rlnorm(1,meanlog = 1, sdlog = 10),
             th = rlnorm(1,meanlog = 1, sdlog = 10),
             Vhp = rlnorm(1,meanlog = 1, sdlog = 10),
             IN = rlnorm(1,meanlog = 1, sdlog = 10),
             q = rlnorm(1,meanlog = 1, sdlog = 10),
             k = rlnorm(1,meanlog = 1, sdlog = 10),
             Vwl = rlnorm(1,meanlog = 1, sdlog = 10),
             tw = rlnorm(1,meanlog = 1, sdlog = 10),
             tp = rlnorm(1,meanlog = 1, sdlog = 10),
             l = rlnorm(1,meanlog = 1, sdlog = 10),
             eh = runif(1, 0.01, 1),
             ew = runif(1, 0.01, 1))
  
  paramlist[[i]] = c(params, N = i)
  
  outlist[[i]] = cbind(test(params), N = i)
  
  yint = test(params)$m4
  names(yint) = c("P", "Iorg", "L","H", "W")
  # yint = yint*c(0.05, 1, 0.05, 1, 0.05)
  if(F){
    mm1 = ode(y = yint[1:3], times = c(seq(1,1.9, 0.1),seq(2,(365*4), length = 10)), func = simpleDCmodel1, parms = params)
    mm2 = ode(y = yint[c(1,2,3,5)], times = c(seq(1,1.9, 0.1),seq(2,(365*4), length = 10)), func = simpleDCmodel2, parms = params)
    mm3 = ode(y = yint[1:4], times = c(seq(1,1.9, 0.1),seq(2,(365*4), length = 10)), func = simpleDCmodel3, parms = params)
    mm4 = ode(y = yint, times = c(seq(1,1.9, 0.1),seq(2,(365*4), length = 10)), func = simpleDCmodel4, parms = params)
    
    mm5 = (mm4[,c(1:4)] - mm3[,c(1:4)] - mm2[,c(1:4)] + mm1[,c(1:4)])
    mm5[,1] = mm4[,1]
    mm5 = data.frame(mm5)
    mm5[,"Effect"] = "IE"
    
    mm6 = (mm2[,c(1:4)] - mm1[,c(1:4)])
    mm6[,1] = mm4[,1]
    mm6 = data.frame(mm6)
    mm6[,"Effect"] = "WE"
    
    mm7 = (mm3[,c(1:4)] - mm1[,c(1:4)])
    mm7[,1] = mm4[,1]
    mm7 = data.frame(mm7)
    mm7[,"Effect"] = "HE"
    
    simlist[[i]] = cbind(rbind(mm5, mm6, mm7), N = i)
  }

  
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

png("Plots/DCmodel.png", width = 8, height = 5, units = "in", res = 600)
par(mfrow=c(1,2))
plot((abs(IE+1e-6))~(abs(WE+1e-6)), data = outlist2, log = 'xy', col = ncol, pch = npch,
     xlab= "Detritivore Effect (log|WE|)", ylab = "Interaction Effect (log|IE|)", cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
legend("topleft", legend = c("Plant", "Litter", "Inorganic N"),
       col = c("#009E73","#E69F00","#56B4E9"), pch = 1:3)
plot((abs(IE+1e-6))~(abs(HE+1e-6)), data = outlist2, log = 'xy', col = ncol, pch = npch,
     xlab= "Herbivore Effect (log|HE|)", ylab = "", cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
dev.off()

# Create plot of simulations over time

simlist = do.call("rbind", simlist)

simlist %>% gather(-time, -N, -Effect, key = StateVar, value = value) %>%
  ggplot(aes(x = time, y = value, color = Effect, group = paste0(N, Effect))) + facet_grid(Effect~StateVar, scales = "free") + geom_line() + theme_classic()

aggregate(cbind(P, Iorg, L) ~ time, data = simlist2, FUN = mean)

# Create the other plot
outlist3 = outlist2

outlist3[,"IEm1"] = signif(abs(outlist2[,"IE"]/ outlist2[,"m1"]), digits = 6)

outlist3 = outlist3[order(outlist3$IEm1),c("nvec", "IEm1")]
dim(outlist3)

tempt = tibble(nvec = rep("A", 30000),
           IEm1 = rep(1, 30000))

tempt$nvec = outlist3$nvec
tempt$IEm1 = outlist3$IEm1

tempt %>% arrange(IEm1) %>% group_by(nvec) %>%
  mutate(cumsum = cumsum(IEm1)) %>% 
  left_join(
    tempt %>% group_by(nvec) %>% summarize(tot = sum(IEm1))
  ) %>%
  mutate(prop = cumsum/tot) %>%
  select(nvec, IEm1, prop) %>%
  distinct() %>%
  ggplot(aes(x = IEm1, y = prop, color = nvec)) + geom_line() + theme_classic()


# Analysis of simple LV model ----

# A function that runs the simple model with Lotka-Volterra functional responses

# m1 = no animals
# m2 = Detritivores
# m3 = Herbivores
# m4 = Herbivores and detritivores

test2 <- function(Ni){
  k = 1 
  
  while(k == 1){ # Draw parameters until a stable version, as judged by the principle eigenvalue of the Jacobian matrix, is found
    
    params = c(Vpn = rlnorm(1,meanlog = 1, sdlog = 10),
               tp = rlnorm(1,meanlog = 1, sdlog = 10),
               th = rlnorm(1,meanlog = 1, sdlog = 10),
               Vhp = rlnorm(1,meanlog = 1, sdlog = 10),
               IN = rlnorm(1,meanlog = 1, sdlog = 10),
               q = rlnorm(1,meanlog = 1, sdlog = 10),
               k = rlnorm(1,meanlog = 1, sdlog = 10),
               Vwl = rlnorm(1,meanlog = 1, sdlog = 10),
               tw = rlnorm(1,meanlog = 1, sdlog = 10),
               tp = rlnorm(1,meanlog = 1, sdlog = 10),
               l = rlnorm(1,meanlog = 1, sdlog = 10),
               eh = runif(1, 0.01, 1),
               ew = runif(1, 0.01, 1))
    
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
out2 = lapply(seq(1,10000,1), FUN = test2)

# Reorganize the data
out3 = vector(mode = "list", length = 10000)

for(j in 1:10000){
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

png("Plots/LVmodel.png", width = 8, height = 5, units = "in", res = 600)
par(mfrow=c(1,2))
plot((abs(IE+1e-6))~(abs(WE+1e-6)), data = out4, log = 'xy', col = ncol, pch = npch,
     xlab= "Earthworm Effect (log|WE|)", ylab = "Interaction Effect (log|IE|)",cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
legend("topleft", legend = c("Plant", "Litter", "Inorganic N"),
       col = c("#009E73","#E69F00","#56B4E9"), pch = 1:3)
plot((abs(IE+1e-6))~(abs(HE+1e-6)), data = out4, log = 'xy', col = ncol, pch = npch,
     xlab= "Herbivore Effect (log|HE|)", ylab = "",cex = ncex)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
dev.off()

# Simulate the Lotka-Volterra simple model over time 

simpleLVmodel1 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg*P - tp*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg*P
    dL = tp*P - k*L - l*L
    
    return(list(c(dP, dIorg, dL)))
  })
}

simpleLVmodel2 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg*P - tp*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg*P + (1 - ew)*Vwl*L*W
    dL = tw*W + tp*P - k*L - l*L - Vwl*L*W
    dW = ew*Vwl*L*W - tw*W
    
    return(list(c(dP, dIorg, dL, dW)))
  })
}

simpleLVmodel3 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg*P - tp*P - Vhp*H*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg*P + (1 - eh)*Vhp*P*H
    dL = th*H + tp*P - k*L - l*L
    dH = eh*Vhp*P*H - th*H
    
    return(list(c(dP, dIorg, dL, dH)))
  })
}

simpleLVmodel4 <- function(t,y, pars){
  with(as.list(c(pars,y)),{
    dP = Vpn*Iorg*P - tp*P - Vhp*H*P
    dIorg = IN - q*Iorg + k*L - Vpn*Iorg*P + (1 - eh)*Vhp*P*H + (1 - ew)*Vwl*L*W
    dL = tw*W + th*H + tp*P - k*L - l*L - Vwl*L*W
    dH = eh*Vhp*P*H - th*H
    dW = ew*Vwl*L*W - tw*W
    
    return(list(c(dP, dIorg, dL, dH, dW)))
  })
}

run1 = test2(1)

params = run1[[2]][1:13]

yint = run1[[1]]$m4[1:5]
names(yint) = run1[[1]]$nvec[1:5]
names(yint)[2] = "Iorg"

mm1 = ode(yint[1:3], times = 1:(365*4), func = simpleLVmodel1, parms = params)
mm2 = ode(yint[c(1,2,3,5)], times = 1:(365*4), func = simpleLVmodel2, parms = params)
mm3 = ode(yint[1:4], times = 1:(365*4), func = simpleLVmodel3, parms = params)
mm4 = ode(yint, times = 1:(365*4), func = simpleLVmodel4, parms = params)

mm5 = mm4[,1:4] - mm3[,1:4] - mm2[,1:4] + mm1[,1:4]
mm5[,1] = mm4[,1]

mm6 = mm2[,1:4] - mm1[,1:4]
mm6[,1] = mm4[,1]

mm7 = mm3[,1:4] - mm1[,1:4]
mm7[,1] = mm4[,1]

run1[[1]][1:3, "m4"] - run1[[1]][1:3, "m3"] -run1[[1]][1:3, "m2"] + run1[[1]][1:3, "m1"]

par(mfrow=c(2,2))
plot(P~time, data = mm5, type = "l", lwd = 2, ylim = c(range(c(mm5[,"P"], mm6[,"P"], mm7[,"P"]))))
points(P~time, data = mm6, type = "l", col = "red", lwd = 2, lty = 2)
points(P~time, data = mm7, type = "l", col = "green", lwd = 2, lty = 3)

plot(Iorg~time, data = mm5, type = "l", lwd = 2, ylim = c(range(c(mm5[,"Iorg"], mm6[,"Iorg"], mm7[,"Iorg"]))))
points(Iorg~time, data = mm6, type = "l", col = "red", lwd = 2, lty = 2)
points(Iorg~time, data = mm7, type = "l", col = "green", lwd = 2, lty = 3)

plot(L~time, data = mm5, type = "l", lwd = 2, ylim = c(range(c(mm5[,"L"], mm6[,"L"], mm7[,"L"]))))
points(L~time, data = mm6, type = "l", col = "red", lwd = 2, lty = 2)
points(L~time, data = mm7, type = "l", col = "green", lwd = 2, lty = 3)

# Plot the most complex model ----

# Load in the equilibrium cluster data (run using the script "complex_model_non-eqm_to_cluster.R")

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

data2c %>% ggplot(aes(x = abs(IE+ 1e-6), fill = Type)) + geom_histogram() + facet_wrap(.~StateVar) + scale_x_log10() + theme_classic()

data2c %>% select(ID, Type, StateVar, IE) %>%
  spread(key = Type, value = IE) %>%
  mutate(DIFF = `NON-EQM`-EQM) %>%
  filter(!is.na(DIFF)) %>%
  ungroup() %>%
  summarize(mean(DIFF), min(DIFF), quantile(DIFF, 0.25),quantile(DIFF, 0.5), quantile(DIFF, 0.75), max(DIFF))

data2c %>% select(ID, Type, StateVar, IE) %>%
  spread(key = Type, value = IE) %>%
  mutate(DIFF = `NON-EQM`-EQM) %>%
  ggplot(aes(x = abs(DIFF)+1e-6)) + geom_histogram() + facet_wrap(.~StateVar) + theme_classic() + scale_x_log10()

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


data.frame(StateVar = c("P", "L", "N", "S", "M"),
           ncol = c("#009E73","#E69F00","#56B4E9","#F0E442", "#0072B2"),
           npch = c(1,2,3,4,5))

# Plot all four models together ----

# A. Simple model, donor-controlled
# B. Simple model, Lotka-Volterra
# C. Complex model, equilibrium
# D. Complex model, non-equilibrium

png("Plots/Bothmodel2.png", width = 8, height = 8, units = "in", res = 600)

par(oma=c(2,2,0,0))
par(mfrow=c(2,2), mar = c(3,3,1,1))
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
arrows(1e13, 1e13, 1e21, 1e2, length = 0.1)
text(1e21, 1e-1, "Interaction effect \n increases")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outlist2, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "A", bty = "n")

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
legend("topleft", legend = "C", bty = "n")

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
legend("topleft", legend = "D", bty = "n")

outfinal3 = read_rds("Data/TrueLinearComplex.rds") # Load in the non-equilibrium cluster data

outfinal3$ncol = as.character(outfinal3$ncol)

plot((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outfinal3, log = 'xy', col = ncol, pch = npch, cex = ncex,
     xlab= "", 
     ylab = "", 
     type = "n", main = "Complex: 4-year treatment",
     axes = F, ylim = c(1e-11, 1e5), xlim = c(1e-11, 1e5))
axis(1, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
axis(2, at = c(1e-11,1e-3,1e5), labels = expression(10^-11,10^-3,10^5))
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
points((abs(IEacc+1e-6))~(abs(IEpred+1e-6)), data = outfinal3, col = ncol, pch = npch, cex = ncex)
legend("topleft", legend = "D", bty = "n")
legend("bottomright", legend = c("Plant", "Litter", "Inorganic N", "Soil", "Microbe"),
       col = c("#009E73","#E69F00","#56B4E9", "#F0E442","#0072B2"), pch = 1:5)

mtext(text="Combined Effect (log|x|)",side=2,line=0,outer=TRUE,cex=1)
mtext(text="Herbivore + Detritivore Effect (log|x|)",side=1,line=0,outer=TRUE,cex=1)
dev.off()