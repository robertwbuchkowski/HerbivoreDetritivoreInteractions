# Analysis of simple DC model

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
    
    return(out)
    
  })
  
}


outlist <- vector(mode = "list", length = 10000)
paramlist <- vector(mode = "list", length = 10000)


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
  
  
}

outlist2 = do.call("rbind", outlist)

outlist2 = outlist2[!(outlist2$nvec %in% c("H", "W")),]

outlist2[,"WE"] = (outlist2$m2 - outlist2$m1)#/outlist2$m1
outlist2[,"HE"] = (outlist2$m3 - outlist2$m1)#/outlist2$m1
outlist2[,"IE"] = (outlist2$m4 - outlist2$m3 - outlist2$m2 + outlist2$m1)#/outlist2$m1

outlist2[,"ncol"] = ifelse(outlist2$nvec == "P", "#009E73",
                           ifelse(outlist2$nvec == "L", "#E69F00",
                                  ifelse(outlist2$nvec == "N", "#56B4E9",
                                  "black")))

outlist2[,"npch"] = ifelse(outlist2$nvec == "P", 1,
                           ifelse(outlist2$nvec == "L", 2,
                                  ifelse(outlist2$nvec == "N", 3,
                                         4)))
png("Plots/DCmodel.png", width = 8, height = 5, units = "in", res = 600)
par(mfrow=c(1,2))
plot((abs(IE+1e-6))~(abs(WE+1e-6)), data = outlist2, log = 'xy', col = ncol, pch = npch,
     xlab= "Earthworm Effect (log|WE|)", ylab = "Interaction Effect (log|IE|)")
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
legend("topleft", legend = c("Plant", "Litter", "Inorganic N"),
       col = c("#009E73","#E69F00","#56B4E9"), pch = 1:3)
plot((abs(IE+1e-6))~(abs(HE+1e-6)), data = outlist2, log = 'xy', col = ncol, pch = npch,
     xlab= "Herbivore Effect (log|HE|)", ylab = "")
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(a = -log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = log(10), b = 1, lty = 2, lwd = 1.5, col = "grey")
abline(a = -log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
abline(a = log(100), b = 1, lty = 3, lwd = 1.5, col = "grey")
dev.off()

# Analysis of simple LV model

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


k = 1 

while(k == 1){
  
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
