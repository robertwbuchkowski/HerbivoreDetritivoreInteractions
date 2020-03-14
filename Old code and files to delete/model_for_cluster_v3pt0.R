# Model for cluster: Version 3.0

source("treatment_apply.R")
require(deSolve)
NTOT = 15
singlerun <- function(idx){
  ID = round(runif(1)*1e8,0) # set round ID
  ID2 = round(runif(1)*1e8,0) # set backup round ID
  paramscur = params # import params
  spinyears = c(500, 1500, 3000) # spin length
  ERRRRR = T # set loop

  while(ERRRRR){
    kkkk = 1
    while(kkkk == 1){
      for(i in 1:31){
        paramscur[i] = rlnorm(1, meanlog = log(params[i]), sdlog = 0.3536)
      }
      
      if(max(paramscur[22:31]) <=1) break
    }
    rm(kkkk)
    
    ystable = yint
    
    for(i in 1:3){
      # Try first quick spin
      stablerun = ode(y=ystable,times = 1:(365*spinyears[i]), 
                      func=plantmodel, 
                      parms=paramscur)
      
      stableannual <- as.data.frame(stablerun)
      stableannual<- stableannual[stableannual$time %in% seq(180, 365*spinyears[i], by=365),]
      
      stabletest <- stableannual[dim(stableannual)[1],c(-1, -11, -12, -13)] -
        stableannual[(dim(stableannual)[1]-1),c(-1, -11, -12, -13)]
      
      ystable = stablerun[365*(spinyears[i]-1),c(-1, -11, -12, -13)]
      
      if(!all(ystable > 1e-4)) break
      if(max(abs(stabletest)) < 1e-3) break
    }
    
    ERRRRR = !all(c(max(abs(stabletest)) < 1e-3, ystable > 1e-4))
  }
  
  ystable2 = ystable
  
  output2_WE_Return = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur)
  
  output2_WE_Remove = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                          events = list(data=eshock_WE))
  
  output2_HW = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                   events = list(data=eadd))
  
  output2_H = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eshock))
  
  ystable["H"] = 0
  
  output2_W = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eadd))
  
  output2_0 = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur,
                  events = list(data=eshock))
  
  output2_WE_Hopper = ode(y=ystable,times = 1:tmax2, func=plantmodel, parms=paramscur)
  
  out = rbind(output2_WE_Return, output2_WE_Hopper,
              output2_WE_Return,output2_WE_Remove,
              output2_HW, output2_H, 
              output2_W, output2_0)
  
  out1 = cbind(out, Treatment = rep(c(7,6,5,4,3,2,1,0), each=tmax2),
               Expt = rep(c(1,1,1,1,2,2,2,2), each=tmax2),
               Run = rep(ID, dim(out)[1]),
               Run2 = rep(ID2, dim(out)[1]))
  
  out2 = out1[out1[,"time"] %in% timestosample,]
  
  lg <- dim(out2)[1]
  
  out3 <- cbind(out2,
                PARS = c(unname(paramscur),i, rep(-2, lg -1 - length(paramscur))),
                YSTABLE = c(unname(ystable), rep(-2, lg - length(ystable2))))
  
  return(out3)

}

repseq = seq(1, NTOT, 1)

out1 = lapply(repseq, FUN=singlerun)

outf <- do.call("rbind", out1)

fname = paste0("modelcluster_", round(runif(1), 7)*10000000,round(runif(1), 7)*10000000, ".csv")

write.csv(outf, fname, row.names = F)