
LF = list.files("Data/")

LF = LF[-8]

CSV = grepl("csv",LF)

dd <- vector("list", length(LF))

for(i in 1:length(LF)){
  if(CSV[i]){
    dd[[i]] = data.frame(file = LF[i],
               column = colnames(read.csv(paste0("Data/", LF[i])))
    )
  }else{
    dd[[i]] = data.frame(file = LF[i],
               column = colnames(readRDS(paste0("Data/", LF[i])))
    )
  }
}

dd
dd <- do.call("rbind", dd)
write.csv(dd, "Data/schema.csv")
