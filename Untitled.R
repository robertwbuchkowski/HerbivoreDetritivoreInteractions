# Random attempt #2

library(FME)
library(tidyverse)

verbose = F

datatomodel2 = read_csv("Data/datatomodel2.csv")

if(verbose){
  datatomodel2 %>%
    ggplot(aes(x=time, y=Biomass)) + geom_point() + 
    facet_grid(StateVar~Treatment, scales="free")
  
  # Use the average biomass in my field treatments to calibrate the model equilibrium
  
  datatomodel2 %>% 
    filter(Treatment == "Rt" & Experiment == "Rob_WE") %>%
    group_by(StateVar) %>%
    summarize(Biomass = mean(Biomass))
}

# Model formulation ----

Buchkowski2019 <-function(t, y,pars){
  
  with(as.list(c(pars,y)),{
    
    if(AAT == 1){
      Temp = 281.9906
    }else{
      Temp = -12.8244*cos(2*3.14/365*(t %% 365)-0.3666)+281.9846
    }
    
    # Temperature-dependent functions
    A_W = exp(-0.25/8.62e-05*(1/Temp-1/288.15)) # Based on Johnston et al. 2018 
    A_P = exp(-2493/Temp)/(1 + (2493/(26712-2493))*exp(26712*(1/(297.65) - 1/Temp))) # Based on Feng
    
    V0 = exp(0.063*(Temp - 273.15) + 5.47) # Based on Wieder et al. 2013
    K0 = exp(0.007*(Temp - 273.15) + 3.19) # Based on Wieder et al. 2013
    
    mfr = (V0*Vlm*L*M/(K0*Klm + M))/((V0*Vsm*S*M/(K0*Ksm + M)) + V0*Vlm*L*M/(Klm + M))
    
    dL = tp*P + tw*W*W + th*H*H + # litter inputs from dead bodies
      (1-a_hp)*A_W*Vhp*H*P - # sloppy herbivore feeding/ urine loss
      V0*Vlm*L*M/(K0*Klm + M) - A_W*Vlw*L*W - # litter consumed by worms and microbes
    q2*L
      
    dM = a_m*(V0*Vlm*L*M/(K0*Klm + M) + V0*Vsm*S*M/(K0*Ksm + M)) - # microbes eat litter and soil
      A_W*Vlw*mfr*M*W - A_W*Vsw*(1-mfr)*M*W - # worms eat microbes in litter and soil
      tm*M # microbes die
    
    dW = a_wl*(A_W*Vlw*L*W) + # worms eat litter
      A_W*Vsw*S*W + # worms eat soil (unassimlated soil stays put)
      A_W*Vlw*mfr*M*W + A_W*Vsw*(1-mfr)*M*W - # worm eat microbes (unassimulated microbes survive gut passage)
      tw*W*W # worms die
    
    dN = IN - q*N - # system gains and losses
      fi*N + fo*S + # exchange with SOM and IORG
      (1-a_m)*(V0*Vlm*L*M/(K0*Klm + M) + V0*Vsm*S*M/(K0*Ksm + M)) - # microbial mineralization
      A_P*Vnp*N*P # plant uptake of nitrogen 
    
    dS = tm*M - # input from dead microbes
      V0*Vsm*S*M/(K0*Ksm + M) - # loss to microbes
      A_W*Vsw*S*W + # loss to worms
      (1-a_wl)*(Vlw*L*W) + # unassimulated worm faeces
      fi*N - fo*S # exchange with IORG
    
    dP = A_P*Vnp*N*P - # nitrogen gain
      tp*P - # density-dependent death
      A_W*Vhp*H*P # herbivory
    
    dP2 = A_P*Vnp*1.5*N*P2 - # nitrogen gain
      tp*P2 - # density-dependent death
      A_W*Vhp*2*H*P2 # herbivory
    
    dH = a_hp*A_W*Vhp*H*P + # herbivory normal plant
      a_hp*A_W*Vhp*H*P2 - # herbivory second plant
      th*H*H # death
    
    return(list(c(dL, dM, dW, dN, dS, dP, dP2, dH)))
    
  }
  )
}

params<- c(Vlm = 0.025/414.4523,
           Klm = 0.005/25.83898,
           Vsm = 0.011/414.4523,
           Ksm = 2500/25.83898,
           Vlw = 0.0013/0.8026427,
           Vsw = 0.01/0.8026427,
           Vnp = 1244.671,
           Vhp = 0.01/0.8026427,
           a_hp = 0.5,
           a_m = 0.50,
           a_wl = 0.3,
           q   = 0.2,
           q2  = 0.2,
           IN = 0.1,
           tm = 0.01,
           tw = 0.001,
           tp = 0.01,
           th = 0.1,
           fi = 0.6,
           fo = 0.002,
           AAT = 1)

yint_base1 = c(L=26, # WE plots with C:N ratio from files
             M=6.73, # WE plots
             W=7.18, # WE plots
             N=0.192, # WE plots
             S=135, # historical data
             P=10.7, # historical data
             P2 = 0, # only use if second plant exists
             H=0.0149) # plot averages

out0_1 <- ode(y = yint_base1, times = 1:1000, parms = params, func = Buchkowski2019); plot(out0_1)

(runout0_1 <- runsteady(y = yint_base1, parms = params, func = Buchkowski2019))

yint_base2 = yint_base1

yint_base2["P2"] = yint_base2["P"] = yint_base2["P"]/2

out0_2 <- ode(y = yint_base2, times = 1:1000, parms = params, func = Buchkowski2019); plot(out0_2)

(runout0_2 <- runsteady(y = yint_base2, parms = params, func = Buchkowski2019))

length(params)
lparams = log(params[1:20])
nparams = params[21]

fitparams <- function(ptf, pntf, input){
  
  ptf2 = c(exp(ptf), pntf)
  
  # run steady works better than stode!
  a <- runsteady(y = input, parms = ptf2, 
                 func = Buchkowski2019)$y
  
  # hmax = 1 was not enough.
  
  Wt = (input-a)/input
  
  return(sum(Wt^2, na.rm= T))
}

# ... Functions for replicating the treatments ----

hopsamp = seq(0,3,1)*365 + 266
soilsamp = c(seq(1,3,1)*365 + 170, seq(0,3,1)*365 + 300)
wormsamp = c(114, 311, c(114, 309, 314)+365, c(99, 280, 323, 114)+365*2)

samp = unique(c(hopsamp, soilsamp, wormsamp))

timestosample = samp[order(samp)]

rm(hopsamp, soilsamp, wormsamp,samp)

# Creates dataframes that can be used as events in ODE simulations

timelist2_expt = sort(c(311,674,1047,114,474,844)) + 365

tmax2 = 1200 + 365

timelist2 = c(timelist2_expt, rep(seq(4,23, by=1), each=2)*365 + rep(c(121,305), 20))

eshock_WE <- data.frame(var = rep("W", length(timelist2)),
                        time =  timelist2,
                        value = rep(0.65, length(timelist2)),
                        method = rep("mult", length(timelist2)))


eshock <- data.frame(var = c("P", "L", "W", rep("W", length(timelist2))),
                     time =  c(rep(290, 3), timelist2),
                     value = c(0.1, 0.1, 0.05, rep(0.2, length(timelist2))),
                     method = rep("mult", length(timelist2)+3))

eadd <- data.frame(var = c("P", "L", "W", rep("W", length(timelist2))),
                   time =  c(rep(290, 3), timelist2),
                   value = c(0.1, 0.1, 0.05, c(0,7.75, 5.73, 9.87, 6.9,0)*0.1, rep(0, 40)),
                   method = c(rep("mult",3),rep("add", length(timelist2))))

rm(timelist2_expt,timelist2)

# Fit the EQM model with 1 plant species -----

(fit1 <- optim(par = lparams, fn = fitparams, 
               pntf = nparams, input = yint_base1, 
               control = list(maxit = 5000)))


if(fit1$convergence!=0){print("WARNING!!!!!! DID NOT CONVERGE")}

(params1 = c(exp(fit1$par), nparams))

params1["a_wl"] = 1

(runout1 <- runsteady(y = yint_base1, parms = params1, func = Buchkowski2019))

out1 <- ode(y = runout1$y, times = 1:10000, parms = params1, func = Buchkowski2019); plot(out1)

params1_1 = params1

params1_1["AAT"] = 0

out1_1 <- ode(y = runout1$y, times = 1:5000, parms = params1_1, func = Buchkowski2019); plot(out1_1)

# Test the difference between models with and without animals ----

testNp = ode(y = runout1$y*c(1,1,0,1,1,1,1,0), times=1:10000, parms = params1_1, func = Buchkowski2019)
if(verbose) plot(testNp)

testHp = ode(y = runout1$y*c(1,1,1,1,1,1,1,0), times=1:10000, parms = params1_1, func = Buchkowski2019)
if(verbose) plot(testHp)

testWp = ode(y = runout1$y*c(1,1,0,1,1,1,1,1), times=1:10000, parms = params1_1, func = Buchkowski2019)
if(verbose) plot(testWp)

testHWp = ode(y = runout1$y*c(1,1,1,1,1,1,1,1), times=1:10000, parms = params1_1, func = Buchkowski2019)
if(verbose) plot(testHWp)

testN = stode(y = runout1$y*c(1,1,0,1,1,1,1,0), parms = params1, func = Buchkowski2019)
testH = stode(y = runout1$y*c(1,1,0,1,1,1,1,1), parms = params1, func = Buchkowski2019)
testW = stode(y = runout1$y*c(1,1,1,1,1,1,1,0), parms = params1, func = Buchkowski2019)
testHW = stode(y = runout1$y*c(1,1,1,1,1,1,1,1), parms = params1, func = Buchkowski2019)

eqbm = data.frame(
  EOWeq = (testW$y - testN$y)/testN$y,
  EOHeq = (testH$y - testN$y)/testN$y,
  IEeq = (testHW$y - testW$y - testH$y + testN$y)/testN$y
)

eqbm = eqbm[c("L", "M", "N", "S", "P"),]


EOWeqp = data.frame((testWp - testNp)/testNp)
EOHeqp = data.frame((testHp - testNp)/testNp)
IEeqp = data.frame((testHWp - testWp - testHp + testNp)/testNp)

oneperyear = T
if(oneperyear){
  sel = seq(180, dim(EOWeqp)[1], 365)
  EOWeqp = EOWeqp[sel,] 
  EOHeqp = EOHeqp[sel,] 
  IEeqp = IEeqp[sel,] 
}

dev.off()

xrange = range(c(EOWeqp$P,eqbm["P", "EOWeq"]))
yrange = range(c(IEeqp$P,eqbm["P", "IEeq"]))

toplot = c("L", "M", "N", "S", "P")

par(mfrow=c(2,5))
for(i in 1:length(toplot)){
  sel = toplot[i]
  plot(IEeqp[,sel]~EOWeqp[,sel], type="l", col="red", main=sel)
  points(IEeqp[1,sel]~EOWeqp[1,sel], type="p", col="red", pch=19)
  
  plot(IEeqp[,sel]~EOHeqp[,sel], type="l", col="red", main=sel)
  points(IEeqp[1,sel]~EOHeqp[1,sel], type="p", col="red", pch=19)
  
}

par(mfrow=c(2,5))
for(i in 1:length(toplot)){
  sel = toplot[i]
  
  xrange = range(c(EOWeqp[,sel],eqbm[sel, "EOWeq"]))
  yrange = range(c(IEeqp[,sel],eqbm[sel, "IEeq"]))
  
  plot(IEeqp[,sel]~EOWeqp[,sel], type="l", col="red", main=sel, ylim = yrange, xlim = xrange)
  points(IEeqp[1,sel]~EOWeqp[1,sel], type="p", col="red", pch=19)
  points(IEeq~EOWeq, data=eqbm[c(sel),], type="p", col="red", pch=1)
  
  xrange = range(c(EOHeqp[,sel],eqbm[sel, "EOHeq"]))
  yrange = range(c(IEeqp[,sel],eqbm[sel, "IEeq"]))
  
  plot(IEeqp[,sel]~EOHeqp[,sel], type="l", col="red", main=sel, ylim = yrange, xlim = xrange)
  points(IEeqp[1,sel]~EOHeqp[1,sel], type="p", col="red", pch=19)
  points(IEeq~EOHeq, data=eqbm[c(sel),], type="p", col="red", pch=1)
  
}


# Fit the 1 plant species model to actual data ----

 rootfun <- function(Time, State, Pars) {
      dstate <- unlist(Buchkowski2019(Time, State, Pars))
   sum(abs(dstate)) - 1e-4
  }

custom_modCost <- function(paramscur1, datatofit, otpt = "Cost"){
  
  paramscur = c(paramscur1, AAT = 1)
  
  # Run the model steady-state
  runout_full_0 <- ode(y = yint_base1, times = 1:10000, parms = paramscur, 
                              func = Buchkowski2019, 
                       rootfun = rootfun)
  
  if(any(is.na(runout_full_0))) return(5000)
  
  
  if(dim(runout_full_0)[1] == 10000) print("Steady run not there at 10000")
  
  if(min(runout_full_0[dim(runout_full_0)[1],-1]) < 0) return(5000)
  
  # Simulate a few years to make sure it is stable
  paramscur["AAT"] = 0 # activate temperature
  
  out_full_0 <- ode(y = runout_full_0[dim(runout_full_0)[1],-1], times = 1:(365*5),
                    parms = paramscur, func = Buchkowski2019)
  # plot(out_full_0)
  
  if(is.na(out_full_0[(365*4 + 180),2])) return(5000)

  error = sum(abs(out_full_0[(365*4 + 180),-1] - out_full_0[(365*3 + 180),-1]))

  if(error > 1e-5) print(paste("Annual run not super stable:", error)) #stop("Not a stable cycle")
  
  ystable = out_full_0[(365*5),-1]
  
  output2_WE_Return = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur)
  
  output2_WE_Remove = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur,
                          events = list(data=eshock_WE))
  
  output2_HW = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur,
                   events = list(data=eadd))
  
  output2_H = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur,
                  events = list(data=eshock))
  
  ystable["H"] = 0
  
  output2_W = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur,
                  events = list(data=eadd))
  
  output2_0 = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur,
                  events = list(data=eshock))
  
  output2_WE_Hopper = ode(y=ystable,times = 1:tmax2, func=Buchkowski2019, parms=paramscur)
  
  out = rbind(output2_WE_Return, output2_WE_Hopper,
              output2_WE_Return,output2_WE_Remove,
              output2_HW, output2_H, 
              output2_W, output2_0)
  
  out = data.frame(out, Treatment = rep(c("RtH","RmH","Rt","RmW","HW","H","W","N"), each=tmax2))
  
  if(otpt == "Cost"){
    datatofit %>%
      left_join(
        datatofit %>%
          group_by(StateVar) %>%
          summarize(Bmad = sd(Biomass)) %>%
          ungroup()
      ) %>%
      left_join(
        out %>% 
          as_tibble() %>%
          gather(-time, -Treatment, key = StateVar, value = Model)
      ) %>% 
      mutate(Diff = abs((Model - Biomass)/Bmad)) %>%
      summarize(Cost = sum(Diff)) %>% unlist() %>% return()
  }else{
    datatofit %>%
      left_join(
        datatofit %>%
          group_by(StateVar) %>%
          summarize(Bmad = sd(Biomass)) %>%
          ungroup()
      ) %>%
      left_join(
        out %>% 
          as_tibble() %>%
          gather(-time, -Treatment, key = StateVar, value = Model)
      ) %>% return()
  }
  
  
}

custom_modCost(params[-21], datatomodel2)

fit_test = optim(par = paramsfit, fn = custom_modCost,
                 datatofit = datatomodel2,
                 lower = rep(0, 20),
                 upper = c(1,1,1,2000,1,
                            1,2000, 1,1,1,
                            1,1,1,10,1,
                            1,1,1,1,1))

write.csv(fit_test$par, "par_fit.csv")

paramsfit = params[1:20]
paramsfit[1:20] = read.csv("par_fit.csv")$x[1:20]

datafit = custom_modCost(paramsfit, datatomodel2, otpt = "Other")

debugonce(custom_modCost)


if(verbose){

  out1 %>% as_tibble() %>% 
    gather(-time, -Treatment, key= StateVar, value = Size) %>% 
    ggplot(aes(x=time, y=Size, color = Treatment)) + geom_line() + 
    facet_wrap(.~StateVar, scales="free") + 
    geom_point(
      data = datatomodel2 %>%
        select(-Experiment) %>%
        group_by(time, StateVar, Treatment) %>%
        summarize(Bio = mean(Biomass)) %>%
        filter(time != 114), 
      aes(x=time, y = Bio)
    )
  
  out1 %>% as_tibble() %>% 
    gather(-time, -Treatment, key= StateVar, value = Size)
  
  datatomodel2 %>%
    select(-Experiment) %>%
    group_by(time, StateVar, Treatment) %>%
    summarize(Bio = mean(Biomass)) %>%
    filter(time != 114) %>%
    left_join(
      out1 %>% as_tibble() %>% 
        gather(-time, -Treatment, key= StateVar, value = Size)
    )  %>% 
    ggplot(aes(x=time, y=Size, color = Treatment)) + geom_line() + 
    facet_wrap(.~StateVar, scales="free") + geom_point(aes(y = Bio))
  
  datatomodel2 %>%
    select(-Experiment) %>%
    group_by(time, StateVar, Treatment) %>%
    summarize(Bio = mean(Biomass)) %>%
    filter(time != 114) %>%
    left_join(
      out1 %>% as_tibble() %>% 
        gather(-time, -Treatment, key= StateVar, value = Size)
    )  %>% 
    ggplot(aes(x=Size, y=Bio, color = Treatment)) + geom_point() + 
    facet_wrap(.~StateVar, scales="free")
  
}

# Fit the model with 2 plant species (NOT WORKING NOW) -----

out2_1 <- ode(y = runout1$y + c(0,0,0,0,0,0, 5, 0), times = 1:5000, parms = params1_1, func = Buchkowski2019); plot(out2_1)


fitparams <- function(ptf, pntf, input){
  
  ptf2 = c(exp(ptf), pntf)
  
  # run steady works better than stode!
  a <- runsteady(y = input, parms = ptf2, func = Buchkowski2019)$y
  
  Wt = (input-a)/input
  
  return(sum(Wt^2, na.rm= T))
}



lparams = log(params[1:24])
nparams = params[25]

(fit2 <- optim(par = lparams, fn = fitparams, pntf = nparams, input = yint_base2, control = list(maxit = 5000)))

if(fit1$convergence!=0){print("WARNING!!!!!! DID NOT CONVERGE")}










params2 = params1

(runout2 <- runsteady(y = yint_base2, parms = params, func = Buchkowski2019))

out2 <- ode(y = yint_base2, times = 1:10000, parms = params2, func = Buchkowski2019); plot(out2)

params3 = params2

params3["AAT"] = 0

out3 <- ode(y = runout2$y, times = 1:5000, parms = params3, func = Buchkowski2019); plot(out3)



