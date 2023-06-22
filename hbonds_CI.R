hb_abierta <- read.delim("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds-details_CI_abierta.dat")
hb_abierta[] <- lapply(hb_abierta[], gsub, pattern = "%", replacement = "", fixed = TRUE)
hb_abierta[] <- lapply(hb_abierta[], gsub, pattern = " ", replacement = "", fixed = TRUE)
hb_abierta[] <- lapply(hb_abierta[], gsub, pattern = "-Main", replacement = "", fixed = TRUE)
hb_abierta[] <- lapply(hb_abierta[], gsub, pattern = "-Side", replacement = "", fixed = TRUE)
hb_abierta$occupancy <- as.numeric(hb_abierta$occupancy)
hb_cerrada <- read.delim("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-details_CI_cerrada.dat")
hb_cerrada[] <- lapply(hb_cerrada[], gsub, pattern = "%", replacement = "", fixed = TRUE)
hb_cerrada[] <- lapply(hb_cerrada[], gsub, pattern = " ", replacement = "", fixed = TRUE)
hb_cerrada[] <- lapply(hb_cerrada[], gsub, pattern = "-Main", replacement = "", fixed = TRUE)
hb_cerrada[] <- lapply(hb_cerrada[], gsub, pattern = "-Side", replacement = "", fixed = TRUE)
hb_cerrada$occupancy <- as.numeric(hb_cerrada$occupancy)
index <- read.delim("~/Desktop/Sysbio/cSTING/corridas restringidas/index.txt")
index$index <- as.numeric(index$index)
index$X6nt6 <- as.character(index$X6nt6)
index$X6nt7_holo <- as.character(index$X6nt7_holo)

for(i in 1:181){
  hb_abierta$donor[i] <- index$index[which(index$X6nt6 == hb_abierta$donor[i])] 
  hb_abierta$acceptor[i] <- index$index[which(index$X6nt6 == hb_abierta$acceptor[i])] 
}
for(i in 1:216){
  hb_cerrada$donor[i] <- index$index[which(index$X6nt7_holo == hb_cerrada$donor[i])]
  hb_cerrada$acceptor[i] <- index$index[which(index$X6nt7_holo == hb_cerrada$acceptor[i])]
}

hb <- merge(x = hb_abierta, y = hb_cerrada, by = c("donor", "acceptor"), all = TRUE)
names(hb)[names(hb)=="occupancy.x"] <- "Abierta"
names(hb)[names(hb)=="occupancy.y"] <- "Cerrada"

hb$Abierta[which(is.na(hb$Abierta))] = 0
hb$Cerrada[which(is.na(hb$Cerrada))] = 0
hb$dif <- hb$Abierta - hb$Cerrada

dif_hb_abierta <- hb[which(hb$dif > 0),]
dif_hb_cerrada <- hb[which(hb$dif < 0),]

for(i in 1:96){
  dif_hb_abierta$donor[i] <- index$X6nt6[as.numeric(dif_hb_abierta$donor[i])]
  dif_hb_abierta$acceptor[i] <- index$X6nt6[as.numeric(dif_hb_abierta$acceptor[i])]
}

for(i in 1:129){
  dif_hb_cerrada$donor[i] <- index$X6nt7_holo[as.numeric(dif_hb_cerrada$donor[i])]
  dif_hb_cerrada$acceptor[i] <- index$X6nt7_holo[as.numeric(dif_hb_cerrada$acceptor[i])]
}

missing <- c("ARG189", "THR190", "ASN191", "PRO192", "MET193", "LEU194", "ARG195", 
             "ALA196", "HIS197", "ARG198", "ASP199", "THR200", "GLU232", "THR233", 
             "ILE234", "LEU235", "THR236", "ARG237", "ALA238", "GLY239", "ILE240", 
             "LYS241", "ARG242", "ARG243", "VAL244", "ASP254", "LYS255", "ASP256", 
             "ASN257", "LYS258", "PRO322", "ALA323", "GLU324", "PRO325", "GLU326", 
             "SER327", "HIS328", "ARG386", "THR387", "ASN388", "PRO389", "MET390", 
             "LEU391", "ARG392", "ALA393", "HIS394", "ARG395", "ASP396", "THR397", 
             "GLU429", "THR430", "ILE431", "LEU432", "THR433", "ARG434", "ALA435", 
             "GLY436", "ILE437", "LYS438", "ARG439", "ARG440", "VAL441", "ASP451", 
             "LYS452", "ASP453", "ASN454", "LYS455", "PRO519", "ALA520", "GLU521", 
             "PRO522", "GLU523", "SER524", "HIS525")
### Dato de color --> Saque todos los puentes de hidrogeno que involucran residuos faltantes

## Puentes de hidrogeno que se pierden despues de la rotacion 
dif_hb_abierta$mis <- "false"
for(i in 1:96){
  if(dif_hb_abierta$donor[i] %in% missing){
    dif_hb_abierta$mis[i] <- "true"
  }
  if(dif_hb_abierta$acceptor[i] %in% missing){
    dif_hb_abierta$mis[i] <- "true"
  }
}
dif_hb_abierta <- dif_hb_abierta[-which(dif_hb_abierta$mis == "true"),]

## Puentes de hidrogeno que se ganan despues de la rotacion 
dif_hb_cerrada$mis <- "false"
for(i in 1:129){
  if(dif_hb_cerrada$donor[i] %in% missing){
    dif_hb_cerrada$mis[i] <- "true"
  }
  if(dif_hb_cerrada$acceptor[i] %in% missing){
    dif_hb_cerrada$mis[i] <- "true"
  }
}

dif_hb_cerrada <- dif_hb_cerrada[-which(dif_hb_cerrada$mis == "true"),]

alter <- c(1:71)
alter <- as.character(alter)

for(i in 1:71){
  if(nrow(dif_hb_cerrada[which(dif_hb_abierta$donor[i] == dif_hb_cerrada[which(dif_hb_abierta$acceptor[i] == dif_hb_cerrada$donor),]$acceptor),]) > 0){
    alter[i] = "true"
  }
}
alter
alter2 <- c(1:92)
alter2 <- as.character(alter2)
for(i in 1:92){
  if(nrow(dif_hb_abierta[which(dif_hb_cerrada$donor[i] == dif_hb_abierta[which(dif_hb_cerrada$acceptor[i] == dif_hb_abierta$donor),]$acceptor),]) > 0){
    alter[i] = "true"
  }
}
alter2
    